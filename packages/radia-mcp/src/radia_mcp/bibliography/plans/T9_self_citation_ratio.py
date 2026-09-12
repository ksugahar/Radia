"""Tier 3 — bibliography_self_citation_ratio.

Compute a bibliography-level surname-match ratio using local screening thresholds.
Surname matching does not establish author identity or a publication verdict.
"""
from __future__ import annotations

import pathlib
import re
import unicodedata

from .._bibparse import read_bib_file, first_author_family, _split_name_parts


def _normalize_name(s: str) -> str:
    s = re.sub(r"\\[A-Za-z]+\*?", "", s)
    s = unicodedata.normalize("NFKD", s).casefold()
    return "".join(c for c in s if c.isalnum())


def _entry_authors(author_field: str) -> list[str]:
    """Return list of normalized lastnames."""
    out = []
    for a in _split_name_parts(author_field or "", r"\s+and\s+"):
        a = a.strip()
        if not a:
            continue
        last = first_author_family(a)
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
    if not needle:
        return "Error: author_lastname must contain name characters"
    try:
        entries = [e for e in read_bib_file(p) if not e.kind.startswith("@")]
    except (OSError, UnicodeError, ValueError) as exc:
        return f"Error: cannot read bibliography: {exc}"
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
    lines.append("  Surname-match screening only; confirm author identity and citation relevance.")
    for k in self_keys:
        lines.append(f"    - {k}")
    if ratio > 0.40:
        lines.append(f"[HIGH] {ratio:.0%} > 40% — above local screening threshold; review relevance")
    elif ratio > 0.25:
        lines.append(f"[MEDIUM] {ratio:.0%} > 25% — above local screening threshold; review relevance")
    else:
        lines.append("PASS — surname-match ratio within local screening threshold")
    return "\n".join(lines)
