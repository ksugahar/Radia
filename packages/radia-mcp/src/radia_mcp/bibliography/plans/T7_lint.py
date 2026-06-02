"""Tier 3 — bibliography_lint.

Check for missing required fields per entry type and report shape
errors. Required-field tables follow BibTeX standard semantics with
lab-specific additions (note field encouraged for narrative).
"""
from __future__ import annotations

import pathlib
import re

from .._bibparse import read_bib_file, is_lab_style_key


# Standard BibTeX required-field tables (minus optional ones).
_REQUIRED = {
    "article":       {"author", "title", "journal", "year"},
    "inproceedings": {"author", "title", "booktitle", "year"},
    "conference":    {"author", "title", "booktitle", "year"},
    "book":          {"author", "title", "publisher", "year"},  # editor OR author
    "incollection":  {"author", "title", "booktitle", "publisher", "year"},
    "techreport":    {"author", "title", "institution", "year"},
    "phdthesis":     {"author", "title", "school", "year"},
    "mastersthesis": {"author", "title", "school", "year"},
    "manual":        {"title"},
    "unpublished":   {"author", "title", "note"},
    "misc":          set(),  # nothing required, anything goes
}
_RECOMMENDED = {
    "article":       {"volume", "pages", "doi"},
    "inproceedings": {"pages", "address"},
    "book":          {"address"},
}


def _missing_fields(entry, required):
    present = set(entry.fields)
    # ``author`` is sometimes provided as ``editor`` for books.
    if "author" in required and "editor" in present:
        required = required - {"author"}
    return required - present


def bibliography_lint(bib_path: str) -> str:
    """Lint a .bib for missing required fields, key shape, year sanity.

    Severity:
        - HIGH: required field missing
        - MEDIUM: cite-key doesn't match lab style ``<author><year><word>``
        - MEDIUM: year is in the future or missing
        - LOW: recommended field missing (for the entry type)
        - LOW: pages use a single hyphen instead of ``--``

    Args:
        bib_path: Path to the .bib source.

    Returns:
        Multi-line lint report.
    """
    p = pathlib.Path(bib_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"

    import datetime
    next_year = datetime.date.today().year + 1

    entries = [e for e in read_bib_file(p) if not e.kind.startswith("@")]
    issues: list[tuple[str, str, str]] = []  # severity, key, message
    for e in entries:
        req = _REQUIRED.get(e.kind, set())
        rec = _RECOMMENDED.get(e.kind, set())
        miss = _missing_fields(e, req)
        for f in miss:
            issues.append(("HIGH", e.key, f"missing required field {f!r} for kind={e.kind}"))
        miss_rec = rec - set(e.fields)
        for f in miss_rec:
            issues.append(("LOW", e.key, f"missing recommended field {f!r}"))

        if not is_lab_style_key(e.key):
            issues.append(("MEDIUM", e.key, "cite-key not in lab style (author+year+word, all lowercase)"))

        year = (e.fields.get("year", "") or "").strip()
        if year:
            try:
                y = int(year[:4])
                if y > next_year:
                    issues.append(("MEDIUM", e.key, f"year {y} > current+1"))
                if y < 1800:
                    issues.append(("LOW", e.key, f"year {y} suspiciously old"))
            except ValueError:
                issues.append(("LOW", e.key, f"year {year!r} not numeric"))

        pages = e.fields.get("pages", "")
        if re.match(r"^\d+\s*-\s*\d+$", pages):
            issues.append(("LOW", e.key, f"pages {pages!r} should use '--' not '-'"))

    sev_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
    issues.sort(key=lambda x: (sev_order[x[0]], x[1]))

    lines = [f"bibliography_lint: {p}", f"  entries: {len(entries)}"]
    if not issues:
        lines.append("PASS — no shape errors")
        return "\n".join(lines)
    lines.append(f"FAIL — {len(issues)} issue(s):")
    for sev, key, msg in issues:
        lines.append(f"  [{sev}] {key}: {msg}")
    return "\n".join(lines)
