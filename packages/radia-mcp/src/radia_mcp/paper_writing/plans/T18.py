"""Paper T18 (Intelligence Layer S4): Workflow agent — paper 用 1 コール
(paper_writing_run_full_workflow).

paper 用の Phase chain:
  Phase 1 intake → Phase 2 adaptive health → Phase 3 root cause + next 5 →
  Phase 4 (reviewer #2 simulate) → Phase 5 rewrite suggest (top action)

revision phase では response letter draft も自動生成。
"""

from __future__ import annotations


def paper_writing_run_full_workflow(
    tex_or_text: str,
    bib: str = "",
    abstract: str = "",
    author_last_names: str = "",
    journal_tier: str = "",
    phase: str = "",
    reviewer_comments: str = "",
    skip_phases: str = "",
    skip_tools: str = "",
) -> dict:
    """paper 用 Phase 1-5 を 1 コール chain 実行。

    Args:
        tex_or_text / bib / abstract / author_last_names: paper T8 へ passthrough
        journal_tier / phase: paper T20 adaptive 用
        reviewer_comments: revision phase 時に response letter 生成用
        skip_phases: カンマ区切り (例 "phase4,phase5")
        skip_tools: T8 の skip 引数

    Returns:
        dict: phases / overall_summary / errors
    """
    if tex_or_text is None:
        tex_or_text = ""

    skip_phases_set = {
        s.strip().lower() for s in skip_phases.split(",") if s.strip()
    }

    output: dict = {"phases": {}, "errors": []}

    # ---- Phase 1 intake ----
    if "phase1" not in skip_phases_set:
        output["phases"]["phase1_intake"] = {
            "char_count": len(tex_or_text),
            "has_bib": bool(bib.strip()),
            "has_abstract": bool(abstract.strip()),
            "context": {"journal_tier": journal_tier, "phase": phase},
        }

    # ---- Phase 2 adaptive health ----
    health = None
    if "phase2" not in skip_phases_set and tex_or_text:
        try:
            from .T20 import paper_writing_adaptive_health_report
            health = paper_writing_adaptive_health_report(
                tex_or_text, bib=bib, abstract=abstract,
                author_last_names=author_last_names,
                journal_tier=journal_tier, phase=phase, skip=skip_tools,
            )
            output["phases"]["phase2_diagnose"] = {
                "overall_score": health.get("overall_score"),
                "n_critical": health.get("n_critical_adjusted", 0),
                "n_high": health.get("n_high_adjusted", 0),
                "adjustments": len(health.get("adjustments_applied", [])),
                "top_3_issues": (
                    health.get("adjusted_priority_issues", [])[:3]
                ),
            }
        except Exception as e:
            output["errors"].append(f"phase2 failed: {e}")

    # ---- Phase 3 root cause + next 5 ----
    next_actions = None
    if "phase3" not in skip_phases_set and tex_or_text:
        try:
            from .T16 import paper_writing_root_cause_diagnosis
            rc = paper_writing_root_cause_diagnosis(
                tex_or_text, bib=bib, abstract=abstract,
                author_last_names=author_last_names, skip=skip_tools,
            )
        except Exception as e:
            rc = None
            output["errors"].append(f"phase3 root_cause failed: {e}")

        try:
            from .T17 import paper_writing_next_5_actions
            next_actions = paper_writing_next_5_actions(
                tex_or_text, bib=bib, abstract=abstract,
                author_last_names=author_last_names, skip=skip_tools,
            )
        except Exception as e:
            output["errors"].append(f"phase3 next_actions failed: {e}")

        output["phases"]["phase3_synthesize"] = {
            "root_causes_detected": rc.get("n_root_causes", 0) if rc else 0,
            "top_3_root_causes": (
                rc.get("detected_root_causes", [])[:3] if rc else []
            ),
            "top_5_actions": (
                next_actions.get("top_5_actions", [])
                if next_actions else []
            ),
        }

    # ---- Phase 4 reviewer #2 simulate (revision 時特に有用) ----
    if "phase4" not in skip_phases_set and tex_or_text:
        try:
            from .T9 import paper_writing_reviewer_2_trigger_summary
            rev2 = paper_writing_reviewer_2_trigger_summary(tex_or_text)
            observations = rev2.get("observations", {})
            active = [
                trigger for trigger in observations.get("triggers", [])
                if (trigger.get("contribution") or 0) > 0
            ]
            output["phases"]["phase4_reviewer_2_simulate"] = {
                "total_triggers": len(active),
                "top_tier_count": len(observations.get("top_triggers", [])),
                "risk_percent": observations.get("risk_percent"),
                "unknown_triggers": observations.get("unknown_triggers", []),
                "summary": (rev2.get("comments") or [""])[0],
            }
        except Exception as e:
            output["errors"].append(f"phase4 failed: {e}")

    # ---- Phase 5 rewrite suggest for top action ----
    if ("phase5" not in skip_phases_set and next_actions
            and next_actions.get("top_5_actions")):
        try:
            from .T19 import paper_writing_rewrite_suggest
            tool_to_target = {
                "T1": "abstract_structure",
                "T2": "contribution_bullets",
                "T3": "statistical_reporting_fix",
                "T4": "limitation_paragraph",
                "T9": "reviewer_2_response_prep",
                "T10": "title_claim_form",
                "T12": "statistical_reporting_fix",
                "T14": "discussion_4_elements",
            }
            top = next_actions["top_5_actions"][0]
            tid = top.get("tool_id", "")
            target = tool_to_target.get(tid)
            if target:
                rw = paper_writing_rewrite_suggest(
                    target=target, n_candidates=3
                )
                output["phases"]["phase5_suggest"] = {
                    "for_top_action": tid,
                    "target": target,
                    "candidates": rw.get("candidates", []),
                    "instruction": rw.get("instruction", ""),
                }
            else:
                output["phases"]["phase5_suggest"] = {
                    "for_top_action": tid,
                    "note": f"{tid} 用 target 未定義、T19 直接参照。",
                }
        except Exception as e:
            output["errors"].append(f"phase5 failed: {e}")

    # ---- Phase 6 (revision のみ): response letter draft ----
    if (phase == "revision" and reviewer_comments
            and "phase6" not in skip_phases_set):
        try:
            from ..tools import paper_writing_generate_response_letter
            # 粗い split (paragraph 単位)
            first_chunk = reviewer_comments.split("\n\n")[0][:500]
            resp = paper_writing_generate_response_letter(
                reviewer_comment=first_chunk,
                agreement_stance="agree",
                response_draft="",
                page_line_ref="",
            )
            output["phases"]["phase6_response_letter"] = {
                "note": "revision phase で reviewer_comments 指定時に生成。"
                        "多くの comments があれば caller が loop。",
                "sample_response": str(resp)[:600],
            }
        except Exception as e:
            output["errors"].append(f"phase6 failed: {e}")

    # Summary
    overall = "Workflow 完了。"
    if health and health.get("n_critical_adjusted", 0) >= 1:
        overall = (
            f"⚠ CRITICAL {health['n_critical_adjusted']} 件。"
            "phase3 の top_5_actions から対応推奨。"
        )
    elif health and health.get("n_high_adjusted", 0) >= 2:
        overall = (
            f"HIGH 課題 {health['n_high_adjusted']} 件。修正で大きく前進。"
        )
    elif health:
        overall = (
            f"overall {health.get('overall_score', '-')}/10 — 致命傷無し。"
            "最終パスは音読で reviewer 視点で判定。"
        )

    output["overall_summary"] = overall
    output["execution_summary"] = {
        "phases_run": list(output["phases"].keys()),
        "phases_skipped": sorted(skip_phases_set),
        "n_errors": len(output["errors"]),
    }
    output["hint"] = (
        "paper T18 は grant T49 の paper 版。Phase 1-6 chain。"
        "revision phase で reviewer_comments 指定時は Phase 6 "
        "response letter draft も自動生成 (first chunk のみ)。"
        "reviewer_2_simulate (Phase 4) は既存 T9 を統合。"
    )
    output["source"] = (
        "Intelligence Layer S4 (paper 版)。paper Phase 1-6 orchestrator。"
        "(paper-writing 2026-04 採用、v0.24.0)。"
    )

    return output
