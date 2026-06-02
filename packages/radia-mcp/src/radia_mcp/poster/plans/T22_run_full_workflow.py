"""Tier 6 — poster_run_full_workflow.

One-call orchestration: intake → adaptive health → root cause →
next-5-actions → rewrite candidates for the top action target.
"""
from __future__ import annotations

from .T17_health_report import poster_health_report
from .T18_root_cause_diagnosis import poster_root_cause_diagnosis
from .T19_next_5_actions import poster_next_5_actions
from .T20_rewrite_suggest import poster_rewrite_suggest
from .T21_adaptive_health_report import poster_adaptive_health_report


def poster_run_full_workflow(tex_path: str,
                             conference: str = "IEEJ",
                             rewrite_target: str = "title") -> str:
    """Chain the Intelligence Layer phases into one call.

    Phase 0: read source.
    Phase 1: composite health score.
    Phase 2: adaptive health report (conference preset).
    Phase 3: root cause diagnosis (5 patterns).
    Phase 4: next-5 actions ranked by impact / effort.
    Phase 5: rewrite candidates for ``rewrite_target``.

    Args:
        tex_path: Path to the poster ``.tex`` source.
        conference: Preset key for adaptive phase (default IEEJ).
        rewrite_target: One of title / subtitle / billboard / abstract /
                        takeaway / qr_label.

    Returns:
        One concatenated report — easy for the user to scroll through.
    """
    blocks = [
        ("# Phase 1 — composite health score",
         poster_health_report(tex_path)),
        ("# Phase 2 — adaptive (conference preset)",
         poster_adaptive_health_report(tex_path, conference=conference)),
        ("# Phase 3 — root cause diagnosis",
         poster_root_cause_diagnosis(tex_path)),
        ("# Phase 4 — top-5 actions",
         poster_next_5_actions(tex_path)),
        ("# Phase 5 — rewrite candidates",
         poster_rewrite_suggest(rewrite_target)),
    ]
    out = [f"poster_run_full_workflow: {tex_path} (conference={conference}, rewrite_target={rewrite_target})", ""]
    for title, body in blocks:
        out.append(title)
        out.append(body)
        out.append("")
    return "\n".join(out)
