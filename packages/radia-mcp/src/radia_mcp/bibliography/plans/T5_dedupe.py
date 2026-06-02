"""Tier 2 — bibliography_dedupe.

Find duplicate entries by DOI (strict) and by fuzzy title (loose).
Useful when merging .bib files from collaborators.
"""
from __future__ import annotations

import pathlib
import re
from difflib import SequenceMatcher

from .._bibparse import read_bib_file


_FUZZY_THRESHOLD = 0.85


def _norm_title(t: str) -> str:
    """Normalize title for fuzzy comparison: lowercase, strip non-alphanum."""
    return re.sub(r"[^a-z0-9]+", " ", t.lower()).strip()


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

    entries = [e for e in read_bib_file(p) if not e.kind.startswith("@")]
    lines = [f"bibliography_dedupe: {p}", f"  entries scanned: {len(entries)}"]

    doi_groups: dict[str, list] = {}
    for e in entries:
        doi = (e.fields.get("doi") or "").strip().lower()
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
            doi1 = (e1.fields.get("doi") or "").lower()
            doi2 = (e2.fields.get("doi") or "").lower()
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

    if not doi_dups and not fuzzy_dups:
        lines.append("PASS — no duplicates detected")
    return "\n".join(lines)
