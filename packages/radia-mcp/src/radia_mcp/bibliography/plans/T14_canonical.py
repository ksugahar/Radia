"""Canonical bibliography access and manuscript-specific ``.bbl`` export.

The bundled ``references.bib`` is the single source of truth. Manuscript
directories do not receive private copies of it: they cite canonical keys and
ship the generated ``.bbl`` required by the publisher.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import shutil
import subprocess
import tempfile

import difflib
from .._bibparse import parse_bib, read_bib_file

CANONICAL = pathlib.Path(__file__).resolve().parents[1] / "data" / "references.bib"

_CITE = re.compile(r"\\(?P<command>(?:no)?cite[a-zA-Z]*|(?:auto|paren|text|foot|smart|super)cite)"
                   r"\*?\s*(?:\[[^\]]*\]\s*)*\{(?P<keys>[^}]*)\}")


def _canonical_snapshot():
    raw = CANONICAL.read_bytes()
    entries = [entry for entry in parse_bib(raw.decode("utf-8", errors="strict")) if entry.key]
    if not entries:
        raise ValueError("canonical bibliography contains no entries")
    keys = [entry.key for entry in entries]
    if len(set(keys)) != len(keys):
        raise ValueError("canonical bibliography contains duplicate citation keys")
    return raw, entries


def _citation_source_sha256(keys: list[str], source: bytes) -> str:
    """Fingerprint cited records/dependencies, not unrelated parent entries.

    Preserve raw value expressions so string-macro and literal values cannot
    collide. Global string/preamble directives are conservatively included.
    """
    from .._source_edit import literal_value

    text = source.decode("utf-8").replace("\r\n", "\n")
    entries = parse_bib(text)
    by_key = {entry.key: entry for entry in entries if entry.key}
    if len(by_key) != sum(bool(entry.key) for entry in entries):
        raise ValueError("duplicate canonical citation keys")
    selected = {}
    visiting = set()

    def visit(key):
        if key in visiting:
            raise ValueError(f"cyclic bibliography dependency: {key}")
        if key in selected:
            return
        if key not in by_key:
            raise ValueError(f"missing bibliography dependency: {key}")
        visiting.add(key)
        entry = by_key[key]
        fields = {name: text[start:end] for name, (start, end) in entry.field_spans.items()}
        for name in ("crossref", "xref", "xdata", "related"):
            if name in fields:
                target = literal_value(fields[name])
                if target is None:
                    raise ValueError(f"bibliography dependency requires literal keys: {key}.{name}")
                for dependency in re.split(r"[,\s]+", target.strip()):
                    if dependency:
                        visit(dependency)
        visiting.remove(key)
        selected[key] = {"kind": entry.kind, "fields": fields}

    for key in keys:
        visit(key)
    payload = {
        "keys": keys, "entries": selected,
        "directives": [(e.kind, e.raw_body) for e in entries if e.kind in ("@string", "@preamble")],
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")).hexdigest()


def bibliography_canonical_path() -> str:
    """Return the bundled canonical bibliography path and entry count."""
    if not CANONICAL.is_file():
        return f"canonical bibliography not found at {CANONICAL}"
    try:
        snapshot, entries = _canonical_snapshot()
    except (OSError, UnicodeError, ValueError) as exc:
        return f"Error: cannot read canonical bibliography: {exc}"
    return (
        "bibliography_canonical_path\n"
        f"  {CANONICAL}\n"
        f"  {len(entries)} entries, {len(snapshot) / 1024:.0f} KB\n"
        "  Cite this file as the single source of truth; ship the generated .bbl."
    )


def bibliography_get_entries(keys: str) -> str:
    """Return canonical BibTeX records as machine-readable JSON.

    ``keys`` is a comma-, semicolon-, or newline-separated list.  Order is
    preserved and duplicates are rejected: a slide or other non-LaTeX
    generator can therefore use symbolic citation keys without keeping a
    private bibliography or hand-maintained rendered references.
    """
    requested = [
        key.strip()
        for key in re.split(r"[,;\n]+", keys)
        if key.strip()
    ]
    if not requested:
        return json.dumps(
            {"ok": False, "error": "give at least one canonical citation key"},
            ensure_ascii=False,
        )
    duplicates = sorted({key for key in requested if requested.count(key) > 1})
    if duplicates:
        return json.dumps(
            {"ok": False, "error": "duplicate citation keys", "keys": duplicates},
            ensure_ascii=False,
        )
    if not CANONICAL.is_file():
        return json.dumps(
            {"ok": False, "error": f"canonical bibliography missing: {CANONICAL}"},
            ensure_ascii=False,
        )

    try:
        snapshot, entries = _canonical_snapshot()
    except (OSError, UnicodeError, ValueError) as exc:
        return json.dumps({"ok": False, "error": f"canonical bibliography unavailable: {exc}"})
    available = {entry.key: entry for entry in entries}
    missing = [key for key in requested if key not in available]
    if missing:
        return json.dumps(
            {
                "ok": False,
                "error": "citation keys are absent from canonical references.bib",
                "missing": missing,
            },
            ensure_ascii=False,
        )

    records = [
        {
            "key": key,
            "kind": available[key].kind,
            "fields": available[key].fields,
        }
        for key in requested
    ]
    return json.dumps(
        {
            "ok": True,
            "canonical_path": str(CANONICAL),
            "canonical_sha256": hashlib.sha256(snapshot).hexdigest(),
            "entries": records,
        },
        ensure_ascii=False,
        indent=2,
    )


def _citation_text(tex: str) -> str:
    """Exclude common literal-code environments and TeX line comments."""
    from ...paper_writing._tex_lex import mask_tex_noncode
    return mask_tex_noncode(tex)


def _keys_in_order(tex: str, include_wildcard: bool = False) -> list[str]:
    """Return citation keys in first-appearance order without duplicates."""
    seen: set[str] = set()
    keys: list[str] = []
    tex = _citation_text(tex)
    if re.search(r"\\(?:auto|paren|text|foot|smart|super)?cites\b|\\(?:InputIfFileExists|includeonly)\b", tex):
        raise ValueError("conditional inputs or multi-cite commands require explicit citation resolution")
    for match in _CITE.finditer(tex):
        before = tex[:match.start()]
        if (len(before) - len(before.rstrip("\\"))) % 2:
            continue
        group = match["keys"]
        for key in (part.strip() for part in group.split(",")):
            if key == "*" and match["command"] != "nocite":
                raise ValueError("wildcard is only supported with nocite")
            if key == "*" and match["command"] == "nocite" and not include_wildcard:
                continue
            if key and key not in seen:
                seen.add(key)
                keys.append(key)
    return keys


def _decode_process_output(data: bytes | str | None) -> str:
    if data is None:
        return ""
    if isinstance(data, str):
        return data
    for encoding in ("utf-8", "cp932"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            pass
    return data.decode("utf-8", errors="replace")


def _generated_keys(data: bytes) -> list[str]:
    """Inspect standard BibTeX bibitems; unknown output syntax fails closed."""
    # Citation keys are ASCII by the input contract; reference text need not be.
    # Latin-1 is a reversible byte mapping; reference text is not reencoded.
    active = _citation_text(data.decode("latin-1")).encode("latin-1")
    heads = list(re.finditer(rb"(?<!\\)\\bibitem\b", active))
    items = list(re.finditer(
        rb"(?<!\\)\\bibitem\s*(?:\[[^\]]*\]\s*)?\{([A-Za-z0-9_:./+-]+)\}", active,
    ))
    if len(heads) != len(items):
        raise ValueError("unsupported or malformed bibitem syntax")
    return [item[1].decode("ascii") for item in items]


def _resolve_installed_style(style: str, directory: pathlib.Path) -> pathlib.Path:
    """Locate the exact installed bst so generation never uses an untracked style."""
    command = shutil.which("kpsewhich")
    if command is None:
        raise ValueError("kpsewhich is required to locate installed styles; provide a local bst instead")
    try:
        result = subprocess.run(
            [command, "--format=bst", "--", f"{style}.bst"], cwd=directory,
            capture_output=True, timeout=15, check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ValueError(f"cannot locate bibliography style: {exc}") from exc
    output = _decode_process_output(result.stdout).strip()
    if result.returncode != 0 or not output or len(output.splitlines()) != 1:
        raise ValueError(f"installed bibliography style could not be resolved: {style}")
    path = pathlib.Path(output)
    if not path.is_absolute():
        path = directory / path
    path = path.resolve()
    if not path.is_file() or path.suffix.casefold() != ".bst":
        raise ValueError(f"resolved bibliography style is not a bst file: {path}")
    return path


def bibliography_make_bbl(
    tex_path: str,
    style: str = "",
    out_path: str | None = None,
    aux_path: str | None = None,
) -> str:
    """Generate canonical bbl from TeX or explicit notebook bibliography keys.

    The caller owns aux freshness and must regenerate it after source changes.
    No TeX or aux commands are executed here. Cooperating writers use an OS lock.
    Notebooks declare metadata.radia.bibliography.keys and optional style; their
    code cells are never executed or scanned for implicit citation identities.
    """
    from .._write_lock import target_lock
    destination = pathlib.Path(out_path).resolve() if out_path else pathlib.Path(tex_path).resolve().with_suffix(".bbl")
    if not destination.parent.is_dir() or destination.suffix.casefold() != ".bbl":
        return "Error: output requires an existing directory and .bbl suffix"
    try:
        with target_lock(destination):
            return _make_bbl_unlocked(tex_path, style, out_path, aux_path)
    except OSError as exc:
        return f"Error: bibliography generation unavailable: {exc}"


def _make_bbl_unlocked(
    tex_path: str,
    style: str = "",
    out_path: str | None = None,
    aux_path: str | None = None,
) -> str:
    """Generate one manuscript's ``.bbl`` from canonical ``references.bib``.

    The manuscript is not compiled. Its citation keys are written to a small
    temporary ``.aux`` file and processed by BibTeX, so missing figures or
    unrelated LaTeX errors do not prevent bibliography generation. If any
    citation key is absent from the canonical file, no partial ``.bbl`` is
    written.
    """
    source = pathlib.Path(tex_path).resolve()
    if not source.is_file():
        return f"Error: no such manuscript: {source}"
    if not CANONICAL.is_file():
        return f"Error: canonical bibliography missing: {CANONICAL}"
    try:
        source_bytes = source.read_bytes()
        tex = source_bytes.decode("utf-8")
    except (UnicodeError, OSError) as exc:
        return f"Error: manuscript is unreadable or not valid UTF-8: {source} ({exc})"
    from ...paper_writing._tex_resolver import resolve_input_chain

    compiled_style = ""
    if source.suffix.casefold() == ".ipynb":
        if aux_path:
            return "Error: compiled aux is not supported for notebook citation keys"
        try:
            notebook = json.loads(tex)
            bibliography = notebook["metadata"]["radia"]["bibliography"]
            keys = bibliography["keys"]
            if (not isinstance(keys, list) or not keys
                    or any(not isinstance(key, str) or not re.fullmatch(r"[A-Za-z0-9_:./+-]+", key) for key in keys)
                    or len(set(keys)) != len(keys)):
                raise ValueError("keys must be a nonempty list of unique safe canonical identifiers")
            if not style:
                style = bibliography.get("style", "IEEEtran")
            if not isinstance(style, str):
                raise ValueError("style must be a string")
        except (KeyError, TypeError, ValueError) as exc:
            return f"Error: invalid notebook bibliography metadata: {exc}"
        resolved = {
            "ok": True,
            "merged_tex": "".join(f"\\cite{{{key}}}" for key in keys),
            "files_resolved": [{"path": str(source), "sha256": hashlib.sha256(source_bytes).hexdigest()}],
        }
    elif aux_path:
        from .._compiled_aux import read_compiled_aux
        try:
            keys, compiled_style, snapshots = read_compiled_aux(pathlib.Path(aux_path))
            snapshots.append({"path": str(source), "sha256": hashlib.sha256(source.read_bytes()).hexdigest()})
        except (OSError, UnicodeError, ValueError) as exc:
            return f"Error: invalid compiled aux: {exc}"
        resolved = {"ok": True, "merged_tex": "", "files_resolved": snapshots}
    else:
        resolved = resolve_input_chain(str(source))
    if not resolved.get("ok"):
        return f"Error: failed to resolve TeX inputs: {resolved.get('error')}"
    try:
        if not aux_path:
            keys = _keys_in_order(resolved["merged_tex"], include_wildcard=True)
    except ValueError as exc:
        return f"Error: unsupported citation syntax: {exc}"
    if not keys:
        return f"Error: no \\cite keys found in {source.name}"

    try:
        canonical_bytes, entries = _canonical_snapshot()
    except (OSError, UnicodeError, ValueError) as exc:
        return f"Error: cannot read canonical bibliography: {exc}"
    known = {entry.key for entry in entries}
    if "*" in keys:
        keys = [key for key in keys if key != "*"]
        keys.extend(entry.key for entry in entries if entry.key not in keys)
    missing = [key for key in keys if key not in known]
    if missing:
        return (
            "Error: citation keys are absent from canonical references.bib: "
            + ", ".join(missing)
            + ". Add and verify them there; do not create a manuscript-local .bib."
        )
    if any(not re.fullmatch(r"[A-Za-z0-9_:./+-]+", key) for key in keys):
        return "Error: citation keys require a safe ASCII BibTeX identifier"

    notebook_fingerprint = ""
    if source.suffix.casefold() == ".ipynb":
        try:
            notebook_fingerprint = _citation_source_sha256(keys, canonical_bytes)
        except ValueError as exc:
            return f"Error: cannot resolve notebook bibliography dependencies: {exc}"

    if not style:
        match = re.search(r"\\bibliographystyle\s*\{([^}]*)\}", _citation_text(resolved["merged_tex"]))
        style = compiled_style or (match.group(1).strip() if match else "IEEEtran")
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+", style):
        return f"Error: unsafe BibTeX style name: {style!r}"

    destination = (
        pathlib.Path(out_path).resolve() if out_path else source.with_suffix(".bbl")
    )
    if destination.suffix.casefold() != ".bbl":
        return f"Error: output path must end in .bbl: {destination}"
    if not destination.parent.is_dir():
        return f"Error: output directory does not exist: {destination.parent}"
    try:
        previous_output = destination.read_bytes() if destination.exists() else None
        local_style = source.parent / f"{style}.bst"
        style_bytes = local_style.read_bytes() if local_style.exists() else None
    except OSError as exc:
        return f"Error: cannot read existing bbl or local bibliography style: {exc}"

    # Validate every input that does not require an external executable first.
    # In particular, an unknown citation key must fail closed even on a host
    # without a TeX installation, and must never replace an existing .bbl.
    bibtex = shutil.which("bibtex")
    if bibtex is None:
        return "Error: bibtex is not on PATH; a TeX installation is required"

    try:
        selected_style = local_style if style_bytes is not None else _resolve_installed_style(style, source.parent)
        selected_style_bytes = style_bytes if style_bytes is not None else selected_style.read_bytes()
    except (OSError, ValueError) as exc:
        return f"Error: cannot snapshot bibliography style: {exc}"

    with tempfile.TemporaryDirectory(prefix="radia-bbl-") as temp_name:
        work = pathlib.Path(temp_name)
        (work / "references.bib").write_bytes(canonical_bytes)

        # A publisher-provided style may live beside the manuscript rather than
        # in the TeX installation. Copy only the requested style.
        (work / local_style.name).write_bytes(selected_style_bytes)

        aux = "\\relax\n" + "".join(f"\\citation{{{key}}}\n" for key in keys)
        aux += f"\\bibstyle{{{style}}}\n\\bibdata{{references}}\n"
        (work / "manuscript.aux").write_text(aux, encoding="ascii")

        try:
            result = subprocess.run(
                [bibtex, "manuscript"],
                cwd=work,
                capture_output=True,
                timeout=180,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            return f"Error: BibTeX could not run: {exc}"

        log = _decode_process_output(result.stdout) + _decode_process_output(
            result.stderr
        )
        generated = work / "manuscript.bbl"
        if result.returncode != 0 or not generated.is_file():
            return "Error: BibTeX did not generate a .bbl:\n" + log[-1600:]

        try:
            data = generated.read_bytes()
            generated_keys = _generated_keys(data)
        except (OSError, ValueError) as exc:
            return f"Error: cannot validate generated bibliography: {exc}"
        bibitem_count = len(generated_keys)
        missing_items = set(keys) - set(generated_keys)
        unknown_items = set(generated_keys) - known
        duplicate_items = len(set(generated_keys)) != bibitem_count
        if missing_items or unknown_items or duplicate_items:
            return (
                "Error: BibTeX generated an incomplete or inconsistent bibliography "
                f"({bibitem_count} bibitems for {len(keys)} cited keys; "
                f"missing={sorted(missing_items)}, unknown={sorted(unknown_items)}, "
                f"duplicate_keys={duplicate_items}).\n"
                + log[-1200:]
            )

        temporary_output = None
        try:
            with tempfile.NamedTemporaryFile(
                dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp", delete=False,
            ) as output_handle:
                temporary_output = pathlib.Path(output_handle.name)
                output_handle.write(data)
            if CANONICAL.read_bytes() != canonical_bytes:
                return "Error: canonical bibliography changed during generation; retry"
            current_style = local_style.read_bytes() if local_style.exists() else None
            if current_style != style_bytes:
                return "Error: local bibliography style changed during generation; retry"
            if selected_style.read_bytes() != selected_style_bytes:
                return "Error: selected bibliography style changed during generation; retry"
            for item in resolved["files_resolved"]:
                current_hash = hashlib.sha256(pathlib.Path(item["path"]).read_bytes()).hexdigest()
                if current_hash != item.get("sha256"):
                    return "Error: TeX input changed during bibliography generation; retry"
            current_output = destination.read_bytes() if destination.exists() else None
            if current_output != previous_output:
                return "Error: destination bbl changed during generation; retry"
            temporary_output.replace(destination)
            temporary_output = None
        except OSError as exc:
            return f"Error: cannot publish generated bbl: {exc}"
        finally:
            if temporary_output is not None:
                temporary_output.unlink(missing_ok=True)

    citation_source = (
        "notebook explicit keys" if source.suffix.casefold() == ".ipynb"
        else "caller-supplied compiled aux (freshness caller-owned)" if aux_path
        else "static TeX scan"
    )
    return (
        f"bibliography_make_bbl: {source.name} -> {destination}\n"
        f"  cited {len(keys)} canonical keys; wrote {bibitem_count} bibitems; "
        f"style {style}\n"
        f"  canonical_sha256: {hashlib.sha256(canonical_bytes).hexdigest()}\n"
        + (f"  selected_source_sha256: {notebook_fingerprint}\n" if notebook_fingerprint else "")
        +
        f"  citation_source: {citation_source}\n"
        f"  style_source: {selected_style}\n"
        f"  style_sha256: {hashlib.sha256(selected_style_bytes).hexdigest()}\n"
        + (f"  local_style_sha256: {hashlib.sha256(style_bytes).hexdigest()}\n"
           if style_bytes is not None else "  installed bibliography style snapshot used\n")
        +
        "  references.bib remained canonical and was not copied into the manuscript."
    )


# Preserve additional read-only bibliography diagnostics.
def bibliography_verify_dois(limit: int = 0, delay: float = 0.3) -> str:
    """Check that every DOI in the canonical file resolves, and to the right work.

    A closed bibliography only stops invented citations if the bibliography
    itself is clean. Two failures are caught here:

      unresolvable -- Crossref has no such DOI. Either a typo or a fabrication;
        this repository carried ``10.1109/TMAG.2018.2848702`` on the Kameari
        2018 CLN paper, which does not exist. The real one is
        ``10.1109/TMAG.2017.2743224``.
      mismatched -- the DOI resolves, but to a different paper than the entry
        claims. The nastier case, because the reference list still compiles and
        looks plausible.
    """
    import time

    from .T1_doi_to_bibtex import _fetch_crossref

    def plain(s: str) -> str:
        s = re.sub(r"\\[a-zA-Z]+", "", s or "")
        return re.sub(r"[^a-z0-9]", "", s.lower())

    entries = [e for e in read_bib_file(CANONICAL) if e.fields.get("doi")]
    if limit:
        entries = entries[:limit]

    dead, wrong, ok, skipped = [], [], 0, 0
    for e in entries:
        doi = e.fields["doi"].strip().rstrip("}").replace("https://doi.org/", "")
        # Registries other than Crossref: arXiv mints through DataCite, and
        # JSIAM and its neighbours through JaLC. An absent record proves
        # nothing about these, so they are not counted as dead.
        if doi.lower().startswith(("10.48550/", "10.11540/", "10.14947/")):
            skipped += 1
            continue
        try:
            msg = _fetch_crossref(doi)
        except Exception:
            msg = None
        time.sleep(delay)
        if not msg:
            dead.append((e.key, doi))
            continue
        theirs = msg.get("title") or [""]
        theirs = plain(theirs[0] if isinstance(theirs, list) else theirs)
        ours = plain(e.fields.get("title", ""))
        # a book listed without its subtitle is the same book, so containment
        # counts as agreement before falling back to a similarity ratio
        # A Japanese title and its English translation share no characters, so
        # a ratio of zero here means the languages differ, not the works.
        ja = re.search(r"[\u3040-\u30ff\u4e00-\u9fff]",
                       e.fields.get("title", "") or "")
        same = (ours and theirs
                and (bool(ja) and theirs.isascii()
                     or ours.startswith(theirs) or theirs.startswith(ours)
                     or difflib.SequenceMatcher(
                         None, ours[:80], theirs[:80]).ratio() >= 0.70))
        # A conference digest and the journal paper that grew out of it share
        # almost the same title, so the words alone never separate them. The
        # container does: Henrotte's axisymmetric entry named IEEE Trans. Magn.
        # while its DOI resolved into the CEFC 1992 digest.
        # "Proceedings of the Royal Society" is a journal, so the word alone
        # proves nothing; what matters is that the container is not the
        # publication the entry names.
        container = " ".join(msg.get("container-title") or [])
        proceedings = re.search(
            r"(?i)digest|proceedings|conference|symposium|workshop", container)
        journal = e.fields.get("journal", "")
        # 年次大会 and "The Proceedings of Mechanical Engineering Congress,
        # Japan" are the same meeting under its two names; comparing scripts
        # that share no characters says nothing about the venue.
        same_venue = bool(journal) and (
            (not journal.isascii() and container.isascii())
            or difflib.SequenceMatcher(
                None, plain(journal)[:40], plain(container)[:40]).ratio() >= 0.55)
        if journal and proceedings and not same_venue \
                and not e.fields.get("booktitle"):
            wrong.append((e.key, doi,
                          "journal: " + e.fields["journal"][:40],
                          "conference: " + container[:40]))
        elif ours and theirs and not same:
            wrong.append((e.key, doi, e.fields.get("title", "")[:52],
                          (msg.get("title") or [""])[0][:52]))
        else:
            ok += 1

    out = [f"bibliography_verify_dois: {len(entries)} entries with a DOI",
           f"  resolved and matching : {ok}",
           f"  DOI does not resolve  : {len(dead)}",
           f"  resolves to another work: {len(wrong)}",
           f"  other registry (arXiv/JaLC, not in Crossref): {skipped}"]
    for k, d in dead[:20]:
        out.append(f"  DEAD  {k}: {d}")
    for k, d, o, t in wrong[:20]:
        out.append(f"  WRONG {k}: {d}")
        out.append(f"          ours: {o}")
        out.append(f"          xref: {t}")
    return "\n".join(out)


def bibliography_refresh_unpublished(min_year: str = "2024",
                                     threshold: float = 0.80) -> str:
    """Report canonical entries that Crossref now knows about.

    Reports only. Review each match before editing the entry: a short or
    generic title can match a different paper.
    """
    from .T3_search_crossref import bibliography_search_crossref

    def plain(s: str) -> str:
        s = re.sub(r"\\[a-zA-Z]+", "", s or "")
        return re.sub(r"[^a-z0-9 ]", "", s.lower()).strip()

    todo = [e for e in read_bib_file(CANONICAL)
            if not e.fields.get("doi") and e.fields.get("title")
            and e.fields.get("year", "0") >= min_year
            and not re.search(r"(?i)web ?site|github", e.fields["title"])]

    out = [f"bibliography_refresh_unpublished: {len(todo)} DOI-less entries "
           f"from {min_year}"]
    hits = 0
    for e in todo:
        title = plain(e.fields["title"])
        if not title:
            continue
        try:
            res = bibliography_search_crossref(title[:110], 3)
        except Exception as exc:                      # network / API trouble
            out.append(f"  ? {e.key}: lookup failed ({exc})")
            continue
        best, line = 0.0, ""
        for cand_line in res.splitlines():
            m = re.match(r"\s*\d+\.\s*\[(\d{4})\]\s*(.+?)\s+—\s+(.+)", cand_line)
            if not m:
                continue
            r = difflib.SequenceMatcher(None, title[:90],
                                        plain(m.group(3))[:90]).ratio()
            if r > best:
                best, line = r, cand_line.strip()
        if best >= threshold:
            hits += 1
            out.append(f"  * {best:.2f} {e.key}")
            out.append(f"      ours: {e.fields['title'][:70]}")
            out.append(f"      xref: {line[:100]}")
    out.append(f"  {hits} candidate(s); verify before editing the entry.")
    return "\n".join(out)


def bibliography_find_stray_bibs(root: str = r"W:\02_学会資料",
                                 max_files: int = 200) -> str:
    """List .bib files outside the canonical one, and whether they can go.

    The lab keeps one bibliography. A local copy in a paper folder is the thing
    that lets a DOI be right in one place and wrong in another, so this reports
    them and, for each, how many entries the canonical file already covers:

      absorbed   -- every entry is already central; the file can be deleted once
                    its .tex points at the canonical bibliography.
      to merge   -- entries not present centrally; merge those first, or they
                    are lost with the file.

    Reports only; it deletes nothing.
    """
    base = pathlib.Path(root)
    if not base.exists():
        return f"no such directory: {base}"

    def ident(e):
        t = re.sub(r"\\[a-zA-Z]+", "", (e.fields.get("title") or "").lower())
        t = re.sub("[^a-z0-9" + CJK_RANGE + "]", "", t)[:80]
        return t or (e.fields.get("doi") or e.key).lower()

    canon = {ident(e) for e in read_bib_file(CANONICAL)
             if not e.kind.startswith("@")}
    rows, unread, total_missing = [], [], 0
    for f in sorted(base.rglob("*.bib")):
        # _retired_bib holds files already taken out of paper folders;
        # counting them would make the migration look stalled
        # _retired_bib holds files already taken out of paper folders;
        # counting them would make the migration look stalled
        if (f.resolve() == CANONICAL.resolve() or ".git" in f.parts
                or "_retired_bib" in f.parts):
            continue
        try:
            entries = [e for e in read_bib_file(f) if not e.kind.startswith("@")]
        except Exception as exc:
            unread.append((f, exc))
            continue
        missing = [e for e in entries if ident(e) not in canon]
        total_missing += len(missing)
        rows.append((f, len(entries), missing))

    out = [f"bibliography_find_stray_bibs: {base}",
           f"  {len(rows)} 個のローカル .bib（正典以外）",
           f"  正典に無いエントリの総数: {total_missing}"]
    ready = [r for r in rows if not r[2]]
    out.append(f"  そのまま削除できる（全件が正典にある）: {len(ready)} ファイル")
    out.append("  ※ 削除前に、その .tex が正典を参照しているか確認すること。")

    for f, n, missing in rows[:max_files]:
        try:
            rel = f.relative_to(base)
        except ValueError:
            rel = f
        if missing:
            out.append(f"  [要マージ {len(missing):3d}/{n:3d}] {rel}")
            for e in missing[:3]:
                t = re.sub(r"[{}]", "", e.fields.get("title", ""))[:56]
                out.append(f"        - {e.key}: {t}")
            if len(missing) > 3:
                out.append(f"        ... 他 {len(missing)-3} 件")
        else:
            out.append(f"  [削除可   {n:3d}件] {rel}")
    if len(rows) > max_files:
        out.append(f"  ... 他 {len(rows)-max_files} ファイル")
    for f, exc in unread:
        out.append(f"  [読めない] {f.name}: {exc}")
    return "\n".join(out)


def bibliography_check_keys(bib_path: str = "") -> str:
    """Report cite keys that cannot survive a real BibTeX run.

    Two failures are silent, which is what makes them worth a check:

      non-ASCII -- classic BibTeX drops the entry with no error, so the
        reference is missing from the printed list while everything looks fine.
      LaTeX-special characters -- & % $ # ~ ^ backslash and braces are acted on
        while \\cite{} is scanned, so the citation never resolves.

    Also lists keys that carry no information about what they cite (POD, ref3),
    which are not broken but collide as soon as a second such paper arrives.
    """
    path = pathlib.Path(bib_path) if bib_path else CANONICAL
    entries = [e for e in read_bib_file(path) if not e.kind.startswith("@")]

    special = set("&%$#{}~^" + chr(92))
    nonascii = [e.key for e in entries if not e.key.isascii()]
    hostile = [(e.key, "".join(sorted(set(e.key) & special)))
               for e in entries if set(e.key) & special]
    vague = [e.key for e in entries
             if re.fullmatch(r"[A-Za-z]{2,6}\d*|ref\d+|bib\d+", e.key)]

    out = [f"bibliography_check_keys: {path.name}, {len(entries)} entries",
           f"  non-ASCII (BibTeX drops these silently): {len(nonascii)}",
           f"  LaTeX-special characters               : {len(hostile)}",
           f"  uninformative (collide on the next one): {len(vague)}"]
    for k in nonascii[:20]:
        out.append(f"  DROPPED  {k}")
    for k, ch in hostile[:20]:
        out.append(f"  UNCITABLE {k}  [{ch}]")
    for k in vague[:20]:
        out.append(f"  vague    {k}")
    if not (nonascii or hostile):
        out.append("  every key survives BibTeX.")
    return "\n".join(out)


def bibliography_search(query: str, limit: int = 10,
                        bib_path: str = "") -> str:
    """Find entries in the lab's bibliography by author, title, year or venue.

    This is the step before citing: it turns "Henrotte's axisymmetric paper"
    into the cite key that make_bbl and the slide tools need. All the terms must
    appear somewhere in the entry, so adding a word narrows rather than widens.

    Searches the lab's own holdings. bibliography_search_crossref looks outward,
    for works not here yet.
    """
    terms = [t for t in re.split(r"[\s,]+", query.lower()) if t]
    if not terms:
        return "bibliography_search: give at least one term"
    path = pathlib.Path(bib_path) if bib_path else CANONICAL
    entries = [e for e in read_bib_file(path) if not e.kind.startswith("@")]

    from .._bibparse import _strip_latex

    def blob(e):
        # An accented author is stored as N\'ed\'elec and typed as "nedelec".
        # Matching the raw field only worked here because the cite key happens
        # to be ASCII; a co-author's name would still have been unreachable.
        parts = [e.key] + [str(v) for v in e.fields.values()]
        return _strip_latex(" ".join(parts)).lower()

    # "the Italian's play model paper" says nothing the entry contains, so the
    # descriptions people actually use are written down alongside it
    from .T15_landmarks import landmark_keys_for
    by_description = landmark_keys_for(terms)

    hits = []
    for e in entries:
        b = blob(e)
        if e.key in by_description or all(t in b for t in terms):
            # a term in the title counts for more than one in a note
            title = (e.fields.get("title") or "").lower()
            score = sum(2 if t in title else 1 for t in terms)
            if e.key in by_description:
                score += 5      # a named landmark is what was being asked for
            hits.append((score, e))
    hits.sort(key=lambda x: (-x[0], x[1].fields.get("year", "")), reverse=False)
    hits.sort(key=lambda x: -x[0])

    def clean(x):
        x = re.sub(r"\\[a-zA-Z]+\s*", "", x or "")
        return " ".join(x.replace("{", "").replace("}", "").split())

    out = [f"bibliography_search: {len(hits)} hit(s) for {query!r}"]
    for _, e in hits[:limit]:
        out.append(f"  {e.key}")
        out.append(f"      {clean(e.fields.get('title', ''))[:78]}")
        out.append(f"      {clean(e.fields.get('author', ''))[:56]}"
                   f"  ({e.fields.get('year', '年不明')})")
        venue = clean(e.fields.get("journal") or e.fields.get("booktitle")
                      or e.fields.get("publisher") or "")
        if venue:
            out.append(f"      {venue[:66]}")
        if e.fields.get("doi"):
            out.append(f"      doi {clean(e.fields['doi'])}")
        if "unverified" in (e.fields.get("note") or ""):
            out.append("      [未検証: \\bibitem からの転記。使う前に確認すること]")
    if len(hits) > limit:
        out.append(f"  ... 他 {len(hits)-limit} 件。語を足すと絞り込める。")
    if not hits:
        out.append("  該当なし。bibliography_search_crossref で外部を探し、"
                   "見つかったら正典に追加すること。")
    return "\n".join(out)
