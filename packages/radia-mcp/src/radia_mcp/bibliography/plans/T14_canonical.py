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

from .._bibparse import parse_bib

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
    """Generate canonical bbl; optional fresh compiled aux resolves macros/conditionals.

    The caller owns aux freshness and must regenerate it after source changes.
    No TeX or aux commands are executed here. Cooperating writers use an OS lock.
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
        tex = source.read_text(encoding="utf-8")
    except (UnicodeError, OSError) as exc:
        return f"Error: manuscript is unreadable or not valid UTF-8: {source} ({exc})"
    from ...paper_writing._tex_resolver import resolve_input_chain

    compiled_style = ""
    if aux_path:
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

    return (
        f"bibliography_make_bbl: {source.name} -> {destination}\n"
        f"  cited {len(keys)} canonical keys; wrote {bibitem_count} bibitems; "
        f"style {style}\n"
        f"  canonical_sha256: {hashlib.sha256(canonical_bytes).hexdigest()}\n"
        f"  citation_source: {'caller-supplied compiled aux (freshness caller-owned)' if aux_path else 'static TeX scan'}\n"
        f"  style_source: {selected_style}\n"
        f"  style_sha256: {hashlib.sha256(selected_style_bytes).hexdigest()}\n"
        + (f"  local_style_sha256: {hashlib.sha256(style_bytes).hexdigest()}\n"
           if style_bytes is not None else "  installed bibliography style snapshot used\n")
        +
        "  references.bib remained canonical and was not copied into the manuscript."
    )
