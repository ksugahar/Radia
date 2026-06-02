"""Tier 3 — bibliography_self_citation_ratio.

Compute the self-citation ratio for a target author. Common journal
reviewer flag: papers with >25% self-citation are usually scrutinized.
"""
from __future__ import annotations

import pathlib
import re

from .._bibparse import read_bib_file


def _normalize_name(s: str) -> str:
    s = s.lower().strip()
    s = re.sub(r"[^a-z]+", "", s)
    return s


def _entry_authors(author_field: str) -> list[str]:
    """Return list of normalized lastnames."""
    out = []
    for a in (author_field or "").split(" and "):
        a = a.strip()
        if not a:
            continue
        if "," in a:
            last = a.split(",", 1)[0]
        else:
            last = a.split()[-1] if a.split() else ""
        out.append(_normalize_name(last))
    return [x for x in out if x]


def bibliography_self_citation_ratio(bib_path: str,
                                     author_lastname: str) -> str:
    """Compute the fraction of bib entries authored by ``author_lastname``.

    Args:
        bib_path: Path to the .bib source.
        author_lastname: Target last name (case-insensitive, e.g.
                         ``"Sugahara"``).

    Returns:
        Count, ratio, list of matched cite-keys + severity hint.
    """
    p = pathlib.Path(bib_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    needle = _normalize_name(author_lastname)
    entries = [e for e in read_bib_file(p) if not e.kind.startswith("@")]
    if not entries:
        return f"bibliography_self_citation_ratio: {p}\n  no normal entries"

    self_keys = []
    for e in entries:
        if needle in _entry_authors(e.fields.get("author", "")):
            self_keys.append(e.key)

    n = len(entries)
    s = len(self_keys)
    ratio = s / n if n else 0
    lines = [f"bibliography_self_citation_ratio: {p}",
             f"  target: {author_lastname!r}",
             f"  matches: {s}/{n} ({ratio:.0%})"]
    for k in self_keys:
        lines.append(f"    - {k}")
    if ratio > 0.25:
        lines.append(f"[MEDIUM] {ratio:.0%} > 25% — most journals flag this as reviewer concern")
    elif ratio > 0.40:
        lines.append(f"[HIGH] {ratio:.0%} > 40% — rejection risk on this axis alone")
    else:
        lines.append("PASS — self-citation within typical range")
    return "\n".join(lines)
