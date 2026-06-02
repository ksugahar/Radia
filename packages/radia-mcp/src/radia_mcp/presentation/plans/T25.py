"""Presentation T25 (Intelligence Layer S4): Workflow agent (presentation 用)
(presentation_run_full_workflow).

Phase 1 intake → Phase 2 adaptive health → Phase 3 root cause + next 5 →
Phase 4 time check → Phase 5 rewrite suggest (top action)。
"""

from __future__ import annotations


def presentation_run_full_workflow(
    pptx_path: str,
    venue: str = "",
    target_duration_min: int = 0,
    skip_phases: str = "",
    skip_tools: str = "",
) -> dict:
    """pptx を 1 コール chain 実行。

    Args:
        pptx_path: .pptx
        venue: '学会発表_短' / '学会発表_長' / '招待講演' / '研究室セミナー' / '授業'
        target_duration_min: 予定発表時間 (分)、0 なら time check skip
        skip_phases / skip_tools: phase / T skip
    """
    skip_set = {s.strip().lower() for s in skip_phases.split(",") if s.strip()}
    output: dict = {"phases": {}, "errors": []}

    # Phase 1 intake
    if "phase1" not in skip_set:
        import pathlib
        p = pathlib.Path(pptx_path)
        output["phases"]["phase1_intake"] = {
            "pptx_exists": p.exists(),
            "pptx_size_mb": round(p.stat().st_size / 1024 / 1024, 2) if p.exists() else 0,
            "venue": venue,
            "target_duration_min": target_duration_min,
        }

    # Phase 2 adaptive health
    health = None
    if "phase2" not in skip_set:
        try:
            from .T27 import presentation_adaptive_health_report
            health = presentation_adaptive_health_report(
                pptx_path, venue=venue, skip=skip_tools,
            )
            output["phases"]["phase2_diagnose"] = {
                "overall_score": health.get("overall_score"),
                "n_critical": health.get("n_critical_adjusted", 0),
                "n_high": health.get("n_high_adjusted", 0),
                "adjustments": len(health.get("adjustments_applied", [])),
                "top_3_issues": health.get("adjusted_priority_issues", [])[:3],
            }
        except Exception as e:
            output["errors"].append(f"phase2: {e}")

    # Phase 3 root cause + next 5
    next_actions = None
    if "phase3" not in skip_set:
        try:
            from .T23 import presentation_root_cause_diagnosis
            rc = presentation_root_cause_diagnosis(pptx_path, skip=skip_tools)
        except Exception as e:
            rc = None
            output["errors"].append(f"phase3 root_cause: {e}")
        try:
            from .T24 import presentation_next_5_actions
            next_actions = presentation_next_5_actions(pptx_path, skip=skip_tools)
        except Exception as e:
            output["errors"].append(f"phase3 next_5: {e}")
        output["phases"]["phase3_synthesize"] = {
            "root_causes_detected": rc.get("n_root_causes", 0) if rc else 0,
            "top_3_root_causes": rc.get("detected_root_causes", [])[:3] if rc else [],
            "top_5_actions": next_actions.get("top_5_actions", []) if next_actions else [],
        }

    # Phase 4 time check
    if ("phase4" not in skip_set and target_duration_min > 0):
        try:
            from ..tools import presentation_estimate_speaking_time
            te = presentation_estimate_speaking_time(pptx_path)
            if isinstance(te, dict):
                est_min = te.get("estimated_minutes") or te.get("estimated_total_minutes") or 0
                output["phases"]["phase4_time_check"] = {
                    "estimated_minutes": est_min,
                    "target_minutes": target_duration_min,
                    "delta_min": round(est_min - target_duration_min, 1),
                    "within_tolerance": abs(est_min - target_duration_min) <= 2,
                    "raw": te,
                }
        except Exception as e:
            output["errors"].append(f"phase4: {e}")

    # Phase 5 rewrite suggest
    if "phase5" not in skip_set and next_actions and next_actions.get("top_5_actions"):
        try:
            from .T26 import presentation_rewrite_suggest
            tool_to_target = {
                "T1": "opening_hook_3sec",
                "T2": "take_home_one_sentence",
                "T13": "slide_title_claim_form",
                "T14": "slide_title_claim_form",
                "T17": "mini_imrad_phase_skeleton",
                "T20": "equation_slide_compliance_template",
                "T21": "figure_slide_compliance_template",
                "T22": "results_slide_stat_template",
            }
            top = next_actions["top_5_actions"][0]
            tid = top.get("tool_id", "")
            target = tool_to_target.get(tid)
            if target:
                rw = presentation_rewrite_suggest(target=target, n_candidates=3)
                output["phases"]["phase5_suggest"] = {
                    "for_top_action": tid,
                    "target": target,
                    "candidates": rw.get("candidates", []),
                    "instruction": rw.get("instruction", ""),
                }
            else:
                output["phases"]["phase5_suggest"] = {
                    "for_top_action": tid,
                    "note": f"{tid} 用 target 未定義、T26 直接参照。",
                }
        except Exception as e:
            output["errors"].append(f"phase5: {e}")

    overall = "Workflow 完了。"
    if health and health.get("n_critical_adjusted", 0) >= 1:
        overall = (
            f"⚠ CRITICAL {health['n_critical_adjusted']} 件。"
            "phase3 top_5 から対応。"
        )
    elif health and health.get("n_high_adjusted", 0) >= 2:
        overall = f"HIGH {health['n_high_adjusted']} 件。修正で前進。"
    elif health:
        overall = (f"overall {health.get('overall_score', '-')}/10 — "
                    "致命傷無し。リハーサルで最終確認。")

    output["overall_summary"] = overall
    output["execution_summary"] = {
        "phases_run": list(output["phases"].keys()),
        "phases_skipped": sorted(skip_set),
        "n_errors": len(output["errors"]),
    }
    output["hint"] = (
        "presentation T25 は grant T49 / paper T18 の pptx 版。"
        "venue + target_duration で adaptive 診断 + time check。"
    )
    output["source"] = "Intelligence Layer S4 (presentation 版)。(v0.25.0)"
    return output
