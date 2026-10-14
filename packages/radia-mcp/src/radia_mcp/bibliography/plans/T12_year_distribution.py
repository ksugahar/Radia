"""Tier 3 — bibliography_year_distribution.

Show the year histogram of citations. Reviewer #2 cares: too-old refs
suggest the lit review is stale, too-young refs suggest the author
skipped foundational work.
"""
from __future__ import annotations

import datetime
import pathlib
import re

from .._bibparse import read_bib_file


def bibliography_year_distribution(bib_path: str) -> str:
    """Print a year histogram of citation entries.

    Args:
        bib_path: Path to the .bib source.

    Returns:
        Year bins + a recency assessment.
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
    now = datetime.date.today().year
    counts: dict[int, int] = {}
    unknown = future = 0
    for e in entries:
        y = (e.fields.get("year", "") or "").strip()
        if not re.fullmatch(r"[0-9]{4}", y) or int(y) == 0:
            unknown += 1
        elif int(y) > now:
            future += 1
        else:
            counts[int(y)] = counts.get(int(y), 0) + 1

    if not counts:
        return (f"bibliography_year_distribution: {p}\n  no valid non-future publication years\n"
                f"  unknown/invalid years: {unknown}; future years: {future}")

    total = sum(counts.values())
    recent_5 = sum(c for y, c in counts.items() if now - y <= 5)
    recent_10 = sum(c for y, c in counts.items() if now - y <= 10)

    lines = [f"bibliography_year_distribution: {p}",
             f"  total dated entries: {total}",
             f"  unknown/invalid years: {unknown}; future years: {future}",
             "  Percentages use only valid non-future dated entries.",
             f"  in last 5 years:  {recent_5}  ({recent_5/total:.0%})",
             f"  in last 10 years: {recent_10} ({recent_10/total:.0%})",
             "  histogram (year: count, bar):"]
    bin_size = 5
    bins: dict[int, int] = {}
    for y, c in counts.items():
        b = (y // bin_size) * bin_size
        bins[b] = bins.get(b, 0) + c
    max_bin = max(bins.values()) if bins else 1
    for b in sorted(bins):
        bar = "#" * round(bins[b] / max_bin * 20)
        lines.append(f"    {b}–{b+4}: {bins[b]:3d}  {bar}")

    # Recency advisory
    if unknown or future:
        lines.append("INCOMPLETE — verify excluded years before assessing literature recency")
    elif recent_5 / total < 0.20:
        lines.append("[MEDIUM] <20% in last 5 years — lit review may look stale")
    if not (unknown or future) and total >= 10 and recent_5 / total > 0.90 and recent_10 / total > 0.95:
        lines.append("[MEDIUM] >90% in last 5y — foundational refs may be missing")
    return "\n".join(lines)
