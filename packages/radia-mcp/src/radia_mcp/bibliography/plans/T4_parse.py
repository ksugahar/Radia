"""Tier 2 — bibliography_parse.

Read a .bib file and return a structured summary. Useful as a quick
inventory tool and as input to other lints.
"""
from __future__ import annotations

import pathlib

from .._bibparse import read_bib_file


def bibliography_parse(bib_path: str) -> str:
    """Parse a .bib file and return per-entry inventory.

    Args:
        bib_path: Path to the ``.bib`` source.

    Returns:
        Tabular summary: ``key | kind | year | first_author | title``.
    """
    p = pathlib.Path(bib_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"
    entries = read_bib_file(p)
    lines = [f"bibliography_parse: {p}",
             f"  total entries: {len(entries)}",
             "  -" * 30]
    for e in entries:
        if e.kind.startswith("@"):
            lines.append(f"  {e.kind}: {e.raw_body[:60]}")
            continue
        author = e.fields.get("author", "")
        first_author = author.split(" and ")[0].split(",")[0].strip() if author else ""
        title = e.fields.get("title", "")[:60]
        year = e.fields.get("year", "?")
        lines.append(f"  {e.key:30s} {e.kind:14s} {year}  {first_author:20s}  {title}")
    return "\n".join(lines)
