"""Tier 5 — bibliography_health_report.

Aggregate all Tier 3 lints into a single weighted score. Same pattern
as poster_health_report.
"""
from __future__ import annotations

import re

from .T5_dedupe import bibliography_dedupe
from .T7_lint import bibliography_lint


_SEVERITY_PENALTY = {"HIGH": 15, "MEDIUM": 6, "LOW": 2}


def _count_severities(report: str) -> dict[str, int]:
    return {sev: len(re.findall(rf"\[{sev}\]", report))
            for sev in ("HIGH", "MEDIUM", "LOW")}


def bibliography_health_report(bib_path: str) -> str:
    """Composite health score for a .bib file.

    Weighted penalty: HIGH=15, MEDIUM=6, LOW=2 (from 100).

    Args:
        bib_path: Path to the .bib source.

    Returns:
        Per-axis severity counts + composite score + verdict.
    """
    reports = {
        "lint": bibliography_lint(bib_path),
        "dedupe": bibliography_dedupe(bib_path),
    }
    score = 100
    axis_sev = {}
    for axis, rep in reports.items():
        sev = _count_severities(rep)
        axis_sev[axis] = sev
        for level, n in sev.items():
            score -= _SEVERITY_PENALTY[level] * n
    score = max(0, score)

    lines = [f"bibliography_health_report: {bib_path}",
             f"  composite score: {score}/100"]
    for axis, sev in axis_sev.items():
        lines.append(f"  - {axis}: H={sev['HIGH']} M={sev['MEDIUM']} L={sev['LOW']}")
    if score >= 80:
        lines.append("VERDICT: GOOD")
    elif score >= 60:
        lines.append("VERDICT: ACCEPTABLE — address HIGH first")
    else:
        lines.append("VERDICT: NEEDS WORK")
    return "\n".join(lines)
