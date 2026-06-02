"""Tier 4 — bibliography_normalize_journal_names.

Rewrite ``journal`` and ``booktitle`` fields to a target style: IEEE
abbreviation, ISO 4 abbreviation, IEICE, or full (de-abbreviated).
"""
from __future__ import annotations

import pathlib

from .._bibparse import read_bib_file, write_bib
from ..data.journal_abbreviations import lookup, reverse_lookup


def bibliography_normalize_journal_names(bib_path: str,
                                         style: str = "ieee",
                                         dry_run: bool = True) -> str:
    """Rewrite journal/booktitle fields to a target style.

    Args:
        bib_path: Path to .bib source.
        style: ``ieee`` (default), ``ieice``, ``iso4``, or ``full``
               (de-abbreviate, where possible).
        dry_run: If True (default), report only. If False, rewrite the
                 .bib in place.

    Returns:
        Per-entry rename summary.
    """
    p = pathlib.Path(bib_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"

    entries = read_bib_file(p)
    changes: list[tuple[str, str, str, str]] = []
    for e in entries:
        if e.kind.startswith("@"):
            continue
        for f in ("journal", "booktitle"):
            old = e.fields.get(f)
            if not old:
                continue
            if style == "full":
                new = reverse_lookup(old) or old
            else:
                new = lookup(old, style) or old
            if new != old:
                e.fields[f] = new
                changes.append((e.key, f, old, new))

    lines = [f"bibliography_normalize_journal_names: {p}",
             f"  style={style}  dry_run={dry_run}",
             f"  proposed changes: {len(changes)}"]
    for key, fld, old, new in changes:
        lines.append(f"  {key}.{fld}: {old!r} → {new!r}")

    if not dry_run and changes:
        p.write_text(write_bib(entries), encoding="utf-8")
        lines.append(f"  wrote: {p}")
    elif not changes:
        lines.append("PASS — no changes needed (already canonical for this style)")
    return "\n".join(lines)
