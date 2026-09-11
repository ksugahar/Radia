"""Tier 2 — bibliography_canonicalize_keys.

Rename cite-keys to the lab convention ``<author><year><word>`` while
preserving existing semantic keywords when they look hand-picked.
"""
from __future__ import annotations

import pathlib
import re
from collections import Counter

from .._bibparse import (BibEntry,
                          make_cite_key, is_lab_style_key)
from .._source_edit import literal_value, read_source, write_source_edits


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

    try:
        original, source, entries = read_source(p)
    except (OSError, UnicodeError, ValueError) as exc:
        return f"Error: cannot read bibliography: {exc}"
    existing_keys = Counter(e.key for e in entries if not e.kind.startswith("@"))
    if any(count > 1 for count in existing_keys.values()):
        return "Error: duplicate input citation keys; no changes made"
    renames: list[tuple[str, str]] = []
    new_entries: list[BibEntry] = []
    for e in entries:
        if e.kind.startswith("@"):
            new_entries.append(e)
            continue
        for field in ("author", "editor", "title"):
            if field in e.field_spans and literal_value(source[slice(*e.field_spans[field])]) is None:
                return (f"Error: {e.key}.{field} requires explicit macro/concatenation "
                        "resolution before generating a citation key; no changes made")
        from .._bibparse import first_author_lastname
        author = first_author_lastname(e.fields.get("author", "")
                                        or e.fields.get("editor", ""))
        year = (e.fields.get("year", "") or "").strip()
        if year and (not re.fullmatch(r"[0-9]{4}", year) or year == "0000"):
            return f"Error: invalid publication year for {e.key!r}; no changes made"
        year = year or "nodate"
        kw = _extract_keyword(e.key, author, year)
        new_key = make_cite_key(e, keyword_override=kw)
        if new_key != e.key:
            renames.append((e.key, new_key))
        new_e = BibEntry(kind=e.kind, key=new_key, fields=e.fields,
                          raw_body=e.raw_body)
        new_entries.append(new_e)

    proposed_keys = Counter(e.key for e in new_entries if not e.kind.startswith("@"))
    collisions = sorted(key for key, count in proposed_keys.items() if count > 1)
    if collisions:
        return f"Error: proposed citation-key collisions {collisions}; no changes made"

    if renames:
        unresolved = []
        for entry in entries:
            for field in ("crossref", "xref", "xdata", "related"):
                if field in entry.field_spans:
                    expression = source[slice(*entry.field_spans[field])]
                    if literal_value(expression) is None:
                        unresolved.append(f"{entry.key}.{field}")
        if unresolved:
            return ("Error: reference macros/concatenations require explicit resolution "
                    "before key renaming: " + ", ".join(unresolved) + "; no changes made")

    lines = [f"bibliography_canonicalize_keys: {p}", f"  total entries: {len(entries)}",
             f"  proposed renames: {len(renames)}",
             f"  dry_run: {dry_run}"]
    for old, new in renames:
        lab_ok = "✓" if is_lab_style_key(new) else "?"
        lines.append(f"  {old!r:35s} -> {new!r:35s} [{lab_ok}]")

    if not dry_run and renames:
        mapping = dict(renames)
        edits = []
        for entry in entries:
            if entry.key in mapping and entry.key_span is not None:
                edits.append((*entry.key_span, mapping[entry.key]))
            for field in ("crossref", "xref", "xdata", "related"):
                value = entry.fields.get(field)
                if not value:
                    continue
                replaced = re.sub(r"[^,\s]+", lambda match: mapping.get(match[0], match[0]), value)
                if replaced != value:
                    edits.append((*entry.field_spans[field], "{" + replaced + "}"))
        try:
            write_source_edits(p, original, source, edits)
        except (OSError, UnicodeError, ValueError) as exc:
            return f"Error: cannot safely write bibliography: {exc}"
        lines.append(f"  wrote: {p}")
        lines.append("  Update external manuscript citations using the rename mapping above.")
    return "\n".join(lines)
