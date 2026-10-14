"""Tier 3 — bibliography_cite_validation.

Cross-check ``\\cite{key}`` references in .tex against keys defined in
.bib. The 2 failure modes:
    - cite-key used in .tex but not in .bib (compilation error)
    - bib key defined but never cited (bloat / forgotten reference)
"""
from __future__ import annotations

import pathlib

from .._bibparse import read_bib_file


def _collect_cited_keys(tex_path: pathlib.Path) -> set[str]:
    """Use the same input resolver and citation scanner as canonical bbl export."""
    from ...paper_writing._tex_resolver import resolve_input_chain
    from .T14_canonical import _keys_in_order

    resolved = resolve_input_chain(str(tex_path))
    if not resolved.get("ok"):
        raise ValueError(f"cannot resolve TeX inputs: {resolved.get('error')}")
    return set(_keys_in_order(resolved["merged_tex"], include_wildcard=True))


def bibliography_cite_validation(tex_path: str, bib_path: str) -> str:
    """Cross-check ``\\cite{}`` keys in .tex against entries in .bib.

    Args:
        tex_path: Root .tex file (follows ``\\input{}`` / ``\\include{}``).
        bib_path: Path to the .bib source.

    Returns:
        - Keys cited but missing from .bib (HIGH — compile error)
        - Keys in .bib but never cited (LOW — bloat)
    """
    tp = pathlib.Path(tex_path)
    if not tp.is_absolute():
        tp = pathlib.Path.cwd() / tp
    bp = pathlib.Path(bib_path)
    if not bp.is_absolute():
        bp = pathlib.Path.cwd() / bp
    if not tp.exists():
        return f"Error: tex file not found: {tp}"
    if not bp.exists():
        return f"Error: bib file not found: {bp}"

    try:
        cited = _collect_cited_keys(tp)
        entries = [e for e in read_bib_file(bp) if not e.kind.startswith("@")]
    except (OSError, UnicodeError, ValueError) as exc:
        return f"Error: citation validation incomplete: {exc}"
    bib_keys = {e.key for e in entries}
    if len(bib_keys) != len(entries):
        return "Error: duplicate bibliography keys make citation validation ambiguous"
    if "*" in cited:
        cited.remove("*")
        cited.update(bib_keys)

    missing = cited - bib_keys
    uncited = bib_keys - cited

    lines = [f"bibliography_cite_validation: tex={tp}, bib={bp}",
             f"  cited keys: {len(cited)}  bib keys: {len(bib_keys)}"]
    if missing:
        lines.append(f"[HIGH] {len(missing)} key(s) cited in .tex but missing from .bib:")
        lines.extend(f"  - {k}" for k in sorted(missing))
    if uncited:
        lines.append(f"[LOW] {len(uncited)} key(s) defined in .bib but never cited:")
        lines.extend(f"  - {k}" for k in sorted(uncited))
    if not missing and not uncited:
        lines.append("PASS — cite-key set matches bib-key set")
    return "\n".join(lines)
