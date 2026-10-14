"""Tier 5 — doc_convert_health_report.

Dispatches to the right lint based on whether the pptx is a *figure*
(embedded into paper/poster) or a *presentation* (shown live). Auto-
detects the kind via ``doc_convert_classify``; can be overridden by the
``kind`` parameter.
"""
from __future__ import annotations

import re

from .T6_slide_lint import doc_convert_slide_lint
from .T7_hidden_audit import doc_convert_hidden_audit
from .T11_figure_readability import doc_convert_figure_readability_check
from .T12_classify import classify_impl


def _count_severities(report: str) -> dict[str, int]:
    return {sev: len(re.findall(rf"\[{sev}\]", report))
            for sev in ("HIGH", "MEDIUM", "LOW")}


def doc_convert_health_report(pptx_path: str,
                          kind: str | None = None,
                          target_width_cm: float = 8.0) -> str:
    """Composite readiness score with kind dispatch.

    Args:
        pptx_path: Source pptx.
        kind: Force ``"figure"`` or ``"presentation"``. If None (default),
              auto-classify via ``doc_convert_classify``.
        target_width_cm: Used only for figure mode; default 8 cm (single
                         column in 2-column journals).

    Returns:
        Per-axis severity counts + composite score + verdict.
    """
    if kind is None:
        kind = classify_impl(pptx_path)
    if kind not in ("figure", "presentation"):
        return f"Error: unknown kind {kind!r}. Use 'figure' or 'presentation'."

    score = 100
    axis_sev: dict[str, dict[str, int]] = {}
    lines = [f"doc_convert_health_report: {pptx_path}", f"  kind: {kind}"]

    if kind == "figure":
        rep = doc_convert_figure_readability_check(pptx_path,
                                                target_width_cm=target_width_cm)
        sev = _count_severities(rep)
        axis_sev["figure_readability"] = sev
        score -= sev["HIGH"] * 15 + sev["MEDIUM"] * 6 + sev["LOW"] * 2
        if "FAIL" in rep:
            score -= 10
    else:
        rep1 = doc_convert_slide_lint(pptx_path)
        rep2 = doc_convert_hidden_audit(pptx_path)
        sev1 = _count_severities(rep1)
        sev2 = _count_severities(rep2)
        axis_sev["slide_lint"] = sev1
        axis_sev["hidden_audit"] = sev2
        score -= sev1["HIGH"] * 15 + sev1["MEDIUM"] * 6 + sev1["LOW"] * 2
        score -= sev2["HIGH"] * 15 + sev2["MEDIUM"] * 6 + sev2["LOW"] * 2
    score = max(0, score)

    lines.append(f"  composite score: {score}/100")
    for axis, sev in axis_sev.items():
        lines.append(f"  - {axis}: H={sev['HIGH']} M={sev['MEDIUM']} L={sev['LOW']}")
    if score >= 80:
        lines.append("VERDICT: GOOD")
    elif score >= 60:
        lines.append("VERDICT: ACCEPTABLE")
    else:
        lines.append("VERDICT: NEEDS WORK")
    return "\n".join(lines)
