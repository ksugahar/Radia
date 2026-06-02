"""Tier 6 — poster_health_report.

Aggregates all Tier 1-5 lints into a single weighted health score
with severity breakdown. The score is a coarse 0-100 number — its main
purpose is to rank top issues, not to be calibrated to anything external.
"""
from __future__ import annotations

import re

from .._textutil import read_tex
from .T1_word_count import poster_word_count
from .T2_line_length import poster_line_length
from .T3_caption_self_contained import poster_caption_self_contained
from .T4_jp_font_check import poster_jp_font_check
from .T5_fontsize_by_distance import poster_fontsize_by_distance
from .T6_typography_lints import poster_typography_lints
from .T7_color_contrast_wcag import poster_color_contrast_wcag
from .T8_colorblind_hint import poster_colorblind_hint
from .T9_color_count_321 import poster_color_count_321


_SEVERITY_PENALTY = {"HIGH": 15, "MEDIUM": 6, "LOW": 2}


_REPORT_FNS = [
    ("word_count", poster_word_count),
    ("line_length", poster_line_length),
    ("caption_self_contained", poster_caption_self_contained),
    ("jp_font", poster_jp_font_check),
    ("fontsize_by_distance", poster_fontsize_by_distance),
    ("typography", poster_typography_lints),
    ("color_contrast", poster_color_contrast_wcag),
    ("colorblind", poster_colorblind_hint),
    ("color_count", poster_color_count_321),
]


def _count_severities(report: str) -> dict[str, int]:
    return {
        sev: len(re.findall(rf"\[{sev}\]", report))
        for sev in ("HIGH", "MEDIUM", "LOW")
    }


def poster_health_report(tex_path: str) -> str:
    """Run Tier 1-2 lints and produce a weighted health score.

    Severity weights: HIGH=15, MEDIUM=6, LOW=2 (deducted from 100).
    Score < 60 = needs work, 60-80 = acceptable, ≥80 = good.

    Returns the per-axis severity counts + the consolidated score +
    pointers to the per-axis tools for drill-down.
    """
    try:
        read_tex(tex_path)
    except FileNotFoundError as e:
        return f"Error: file not found: {e}"

    score = 100
    axis_severities = {}
    for axis, fn in _REPORT_FNS:
        try:
            rep = fn(tex_path)
        except Exception as e:
            axis_severities[axis] = {"error": str(e)}
            continue
        sev = _count_severities(rep)
        axis_severities[axis] = sev
        for level, n in sev.items():
            score -= _SEVERITY_PENALTY[level] * n
    score = max(0, score)

    lines = [f"poster_health_report: {tex_path}",
             f"  composite score: {score}/100"]
    for axis, sev in axis_severities.items():
        if "error" in sev:
            lines.append(f"  - {axis}: ERROR ({sev['error']})")
            continue
        h, m, l = sev["HIGH"], sev["MEDIUM"], sev["LOW"]
        lines.append(f"  - {axis}: H={h} M={m} L={l}")
    if score >= 80:
        lines.append("VERDICT: GOOD — poster is in submission shape")
    elif score >= 60:
        lines.append("VERDICT: ACCEPTABLE — address HIGH/MEDIUM before printing")
    else:
        lines.append("VERDICT: NEEDS WORK — fix HIGH severity items first")
    return "\n".join(lines)
