"""Tier 2 — bibliography_canonicalize_keys.

Rename cite-keys to the lab convention ``<author><year><word>`` while
preserving existing semantic keywords when they look hand-picked.
"""
from __future__ import annotations

import pathlib
import re

from .._bibparse import (BibEntry, read_bib_file, write_bib,
                          make_cite_key, is_lab_style_key)


def _extract_keyword(key: str, author: str, year: str) -> str | None:
    """Pull the semantic word out of an existing lab-style key.

    For ``bobbio1997play`` with author=bobbio, year=1997 → returns ``play``.
    Returns None if the key doesn't have the expected prefix.
    """
    expected_prefix = f"{author}{year}"
    if key.startswith(expected_prefix):
        return key[len(expected_prefix):] or None
    return None


def bibliography_canonicalize_keys(bib_path: str,
                                   dry_run: bool = True) -> str:
    """Rewrite cite-keys to the lab style, preserving hand-picked words.

    Algorithm per entry:
        1. Compute ``new_prefix = author + year``.
        2. If old key already starts with the new prefix, **keep its
           semantic suffix** (it was hand-picked).
        3. Otherwise, fall back to ``first_title_word`` to generate one.

    Args:
        bib_path: Path to .bib source.
        dry_run: If True (default), only report the proposed renames.
                 If False, write the renamed .bib back to disk.

    Returns:
        Per-entry rename summary.
    """
    p = pathlib.Path(bib_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"

    entries = read_bib_file(p)
    renames: list[tuple[str, str]] = []
    new_entries: list[BibEntry] = []
    for e in entries:
        if e.kind.startswith("@"):
            new_entries.append(e)
            continue
        from .._bibparse import first_author_lastname
        author = first_author_lastname(e.fields.get("author", "")
                                        or e.fields.get("editor", ""))
        year = (e.fields.get("year", "") or "").strip()[:4] or "nodate"
        kw = _extract_keyword(e.key, author, year)
        new_key = make_cite_key(e, keyword_override=kw)
        if new_key != e.key:
            renames.append((e.key, new_key))
        new_e = BibEntry(kind=e.kind, key=new_key, fields=e.fields,
                          raw_body=e.raw_body)
        new_entries.append(new_e)

    lines = [f"bibliography_canonicalize_keys: {p}", f"  total entries: {len(entries)}",
             f"  proposed renames: {len(renames)}",
             f"  dry_run: {dry_run}"]
    for old, new in renames:
        lab_ok = "✓" if is_lab_style_key(new) else "?"
        lines.append(f"  {old!r:35s} -> {new!r:35s} [{lab_ok}]")

    if not dry_run and renames:
        out_path = p
        out_path.write_text(write_bib(new_entries), encoding="utf-8")
        lines.append(f"  wrote: {out_path}")
    return "\n".join(lines)
