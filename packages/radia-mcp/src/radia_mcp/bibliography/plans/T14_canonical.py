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

from .._bibparse import read_bib_file

CANONICAL = pathlib.Path(__file__).resolve().parents[1] / "data" / "references.bib"

_CITE = re.compile(r"\\(?P<command>(?:no)?cite[a-zA-Z]*|(?:auto|paren|text|foot|smart|super)cite)"
                   r"\*?\s*(?:\[[^\]]*\]\s*)*\{(?P<keys>[^}]*)\}")


def bibliography_canonical_path() -> str:
    """Return the bundled canonical bibliography path and entry count."""
    if not CANONICAL.is_file():
        return f"canonical bibliography not found at {CANONICAL}"
    entries = [entry for entry in read_bib_file(CANONICAL) if entry.key]
    return (
        "bibliography_canonical_path\n"
        f"  {CANONICAL}\n"
        f"  {len(entries)} entries, {CANONICAL.stat().st_size / 1024:.0f} KB\n"
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

    available = {
        entry.key: entry
        for entry in read_bib_file(CANONICAL)
        if entry.key
    }
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
            "canonical_sha256": hashlib.sha256(CANONICAL.read_bytes()).hexdigest(),
            "entries": records,
        },
        ensure_ascii=False,
        indent=2,
    )


def _citation_text(tex: str) -> str:
    """Exclude common literal-code environments and TeX line comments."""
    tex = re.sub(r"\\begin\{(verbatim\*?|Verbatim|lstlisting|minted)\}.*?\\end\{\1\}",
                 "", tex, flags=re.DOTALL)
    tex = re.sub(r"\\verb\*?([^\w\s]).*?\1", "", tex)
    lines = []
    for line in tex.splitlines(keepends=True):
        for index, char in enumerate(line):
            if char != "%":
                continue
            before = line[:index]
            slashes = len(before) - len(before.rstrip("\\"))
            if slashes % 2 == 0:
                line = before + "\n"
                break
        lines.append(line)
    return "".join(lines)


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


def bibliography_make_bbl(
    tex_path: str,
    style: str = "",
    out_path: str | None = None,
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
    except UnicodeDecodeError as exc:
        return f"Error: manuscript is not valid UTF-8: {source} ({exc})"
    from ...paper_writing._tex_resolver import resolve_input_chain

    resolved = resolve_input_chain(str(source))
    if not resolved.get("ok"):
        return f"Error: failed to resolve TeX inputs: {resolved.get('error')}"
    try:
        keys = _keys_in_order(resolved["merged_tex"], include_wildcard=True)
    except ValueError as exc:
        return f"Error: unsupported citation syntax: {exc}"
    if not keys:
        return f"Error: no \\cite keys found in {source.name}"

    entries = [entry for entry in read_bib_file(CANONICAL) if entry.key]
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
        style = match.group(1).strip() if match else "IEEEtran"
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+", style):
        return f"Error: unsafe BibTeX style name: {style!r}"

    destination = (
        pathlib.Path(out_path).resolve() if out_path else source.with_suffix(".bbl")
    )
    if destination.suffix.casefold() != ".bbl":
        return f"Error: output path must end in .bbl: {destination}"
    if not destination.parent.is_dir():
        return f"Error: output directory does not exist: {destination.parent}"

    # Validate every input that does not require an external executable first.
    # In particular, an unknown citation key must fail closed even on a host
    # without a TeX installation, and must never replace an existing .bbl.
    bibtex = shutil.which("bibtex")
    if bibtex is None:
        return "Error: bibtex is not on PATH; a TeX installation is required"

    with tempfile.TemporaryDirectory(prefix="radia-bbl-") as temp_name:
        work = pathlib.Path(temp_name)
        shutil.copyfile(CANONICAL, work / "references.bib")

        # A publisher-provided style may live beside the manuscript rather than
        # in the TeX installation. Copy only the requested style.
        local_style = source.parent / f"{style}.bst"
        if local_style.is_file():
            shutil.copyfile(local_style, work / local_style.name)

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

        data = generated.read_bytes()
        bibitem_count = data.count(b"\\bibitem")
        if bibitem_count < len(keys):
            return (
                "Error: BibTeX generated an incomplete bibliography "
                f"({bibitem_count} bibitems for {len(keys)} cited keys).\n"
                + log[-1200:]
            )

        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as output_handle:
            output_handle.write(data)
            temporary_output = pathlib.Path(output_handle.name)
        temporary_output.replace(destination)

    return (
        f"bibliography_make_bbl: {source.name} -> {destination}\n"
        f"  cited {len(keys)} canonical keys; wrote {bibitem_count} bibitems; "
        f"style {style}\n"
        "  references.bib remained canonical and was not copied into the manuscript."
    )
