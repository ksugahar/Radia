"""Canonical bibliography access and manuscript-specific ``.bbl`` export.

The bundled ``references.bib`` is the single source of truth. Manuscript
directories do not receive private copies of it: they cite canonical keys and
ship the generated ``.bbl`` required by the publisher.
"""
from __future__ import annotations

import pathlib
import re
import shutil
import subprocess
import tempfile

from .._bibparse import read_bib_file

CANONICAL = pathlib.Path(__file__).resolve().parents[1] / "data" / "references.bib"

_CITE = re.compile(r"\\(?:no)?cite[a-zA-Z]*\s*(?:\[[^\]]*\]\s*)*\{([^}]*)\}")


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


def _keys_in_order(tex: str) -> list[str]:
    """Return citation keys in first-appearance order without duplicates."""
    seen: set[str] = set()
    keys: list[str] = []
    for group in _CITE.findall(tex):
        for key in (part.strip() for part in group.split(",")):
            if key and key != "*" and key not in seen:
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
    bibtex = shutil.which("bibtex")
    if bibtex is None:
        return "Error: bibtex is not on PATH; a TeX installation is required"

    try:
        tex = source.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        return f"Error: manuscript is not valid UTF-8: {source} ({exc})"
    from ...paper_writing._tex_resolver import resolve_input_chain

    resolved = resolve_input_chain(str(source))
    if not resolved.get("ok"):
        return f"Error: failed to resolve TeX inputs: {resolved.get('error')}"
    keys = _keys_in_order(resolved["merged_tex"])
    if not keys:
        return f"Error: no \\cite keys found in {source.name}"

    entries = [entry for entry in read_bib_file(CANONICAL) if entry.key]
    known = {entry.key for entry in entries}
    missing = [key for key in keys if key not in known]
    if missing:
        return (
            "Error: citation keys are absent from canonical references.bib: "
            + ", ".join(missing)
            + ". Add and verify them there; do not create a manuscript-local .bib."
        )

    if not style:
        match = re.search(r"\\bibliographystyle\s*\{([^}]*)\}", tex)
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
