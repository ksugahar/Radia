"""Tier 2 — bibliography_dedupe.

Find duplicate entries by DOI (strict) and by fuzzy title (loose).
Useful when merging .bib files from collaborators.
"""
from __future__ import annotations

import pathlib
import re
import unicodedata
from collections import Counter
from difflib import SequenceMatcher

from .._bibparse import read_bib_file
from .._doi import normalize_doi


_FUZZY_THRESHOLD = 0.85


def _norm_title(t: str) -> str:
    """Normalize title for fuzzy comparison: lowercase, strip non-alphanum."""
    t = unicodedata.normalize("NFKC", t).casefold()
    return " ".join("".join(c if c.isalnum() else " " for c in t).split())


def bibliography_dedupe(bib_path: str) -> str:
    """Find duplicate entries in a .bib file.

    Two pass:
        1. Exact DOI match → HIGH severity (almost certainly same paper)
        2. Fuzzy title similarity ≥0.85 → MEDIUM severity (probably same)

    Args:
        bib_path: Path to .bib source.

    Returns:
        Per-duplicate pair report.
    """
    p = pathlib.Path(bib_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"

    try:
        entries = [e for e in read_bib_file(p) if not e.kind.startswith("@")]
    except (OSError, UnicodeError, ValueError) as exc:
        return f"Error: cannot read bibliography: {exc}"
    lines = [f"bibliography_dedupe: {p}", f"  entries scanned: {len(entries)}"]
    repeated_keys = {key: count for key, count in Counter(e.key for e in entries).items() if count > 1}
    for key, count in repeated_keys.items():
        lines.append(f"[HIGH] Citation-key collision: {key!r} occurs {count} times")

    doi_groups: dict[str, list] = {}
    for e in entries:
        doi = normalize_doi(e.fields.get("doi") or "").casefold()
        if doi:
            doi_groups.setdefault(doi, []).append(e)

    doi_dups = [(d, gs) for d, gs in doi_groups.items() if len(gs) > 1]
    if doi_dups:
        lines.append("[HIGH] DOI duplicates:")
        for doi, gs in doi_dups:
            keys = [e.key for e in gs]
            lines.append(f"  - {doi}: keys = {keys}")

    fuzzy_dups = []
    for i, e1 in enumerate(entries):
        t1 = _norm_title(e1.fields.get("title", ""))
        if not t1:
            continue
        for e2 in entries[i + 1:]:
            if e1.key == e2.key:
                continue
            doi1 = normalize_doi(e1.fields.get("doi") or "").casefold()
            doi2 = normalize_doi(e2.fields.get("doi") or "").casefold()
            if doi1 and doi2 and doi1 == doi2:
                continue  # already in doi_dups
            t2 = _norm_title(e2.fields.get("title", ""))
            if not t2:
                continue
            ratio = SequenceMatcher(None, t1, t2).ratio()
            if ratio >= _FUZZY_THRESHOLD:
                fuzzy_dups.append((e1.key, e2.key, ratio,
                                    e1.fields.get("title", "")[:40]))

    if fuzzy_dups:
        lines.append(f"[MEDIUM] Fuzzy-title near-duplicates (≥{_FUZZY_THRESHOLD}):")
        for k1, k2, r, t in fuzzy_dups:
            lines.append(f"  - {k1!r} ↔ {k2!r}: {r:.2f} — {t}")

    if not doi_dups and not fuzzy_dups and not repeated_keys:
        lines.append("PASS — no duplicates detected")
    return "\n".join(lines)
