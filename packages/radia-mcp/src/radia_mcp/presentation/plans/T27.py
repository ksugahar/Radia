"""Presentation T27 (Intelligence Layer S6): Adaptive health report (presentation 用)
(presentation_adaptive_health_report).

presentation_context (venue: 学会発表_短 / 学会発表_長 / 招待講演 / 研究室セミナー /
授業) で severity を adjust。
"""

from __future__ import annotations


ADAPTIVE_RULES: list[dict] = [
    # 学会発表 短 (10-15 min): 密度 + take-home 最重視
    {
        "context_match": {"venue": ("学会発表_短", "conference_short", "poster")},
        "rules": [
            {"tool": "T2", "score_below": 8, "upgrade_to": "CRITICAL"},
            {"tool": "T11", "score_below": 6, "upgrade_to": "HIGH"},
            {"tool": "T15", "score_below": 6, "upgrade_to": "HIGH"},
            # backup / speaker note は緩め
            {"tool": "T7", "skip": True},
        ],
    },
    # 学会発表 長 (20-30 min): discussion + Q&A が重要
    {
        "context_match": {"venue": ("学会発表_長", "conference_long", "regular")},
        "rules": [
            {"tool": "T17", "score_below": 7, "upgrade_to": "HIGH"},
            {"tool": "T22", "score_below": 7, "upgrade_to": "HIGH"},
        ],
    },
    # 招待講演 (30-60 min): vision + context 広め
    {
        "context_match": {"venue": ("招待講演", "keynote", "invited")},
        "rules": [
            {"tool": "T1", "score_below": 8, "upgrade_to": "CRITICAL"},
            {"tool": "T2", "score_below": 8, "upgrade_to": "CRITICAL"},
            # 密度は緩め (時間あるので)
            {"tool": "T11", "score_below": 5, "downgrade_to": "MEDIUM"},
            {"tool": "T16", "score_below": 5, "downgrade_to": "MEDIUM"},
        ],
    },
    # 研究室セミナー: 深い議論前提、formality 緩め
    {
        "context_match": {"venue": ("研究室セミナー", "lab_seminar", "informal")},
        "rules": [
            # logo / progress indicator 不要
            {"tool": "T4", "skip": True},
            {"tool": "T5", "skip": True},
            # 理系 compliance は必要、その他緩和
            {"tool": "T16", "score_below": 4, "downgrade_to": "MEDIUM"},
        ],
    },
    # 授業: 多彩な slide + 学生向け簡素化
    {
        "context_match": {"venue": ("授業", "lecture", "class")},
        "rules": [
            {"tool": "T6", "score_below": 6, "upgrade_to": "HIGH"},  # visual ratio
            # 数式 compliance 厳しめ (学生理解のため)
            {"tool": "T20", "score_below": 8, "upgrade_to": "HIGH"},
            # 統計 / DMP は緩め
            {"tool": "T22", "skip": True},
        ],
    },
]


def _matches(rule_ctx, ctx):
    for k, expected in rule_ctx.items():
        actual = ctx.get(k)
        if not actual:
            return False
        actual_lc = actual.lower() if isinstance(actual, str) else actual
        if isinstance(expected, tuple):
            return any(
                (actual_lc == e or (isinstance(actual, str) and e in actual_lc))
                for e in expected
            )
        return actual_lc == (expected.lower() if isinstance(expected, str) else expected)


def presentation_adaptive_health_report(
    pptx_path: str, venue: str = "", skip: str = ""
) -> dict:
    """pptx health_report の severity を venue で adjust。

    venue: 学会発表_短 / 学会発表_長 / 招待講演 / 研究室セミナー / 授業 / ''
    """
    try:
        from .T12 import presentation_health_report
    except Exception as e:
        return {"error": f"health_report import failed: {e}",
                "source": "presentation T27 adaptive"}

    health = presentation_health_report(pptx_path, skip=skip)
    ctx = {"venue": venue}
    detailed_scores = health.get("detailed_scores", {})
    priority = list(health.get("priority_issues", []))
    adjustments: list[dict] = []
    skipped_tools: list[str] = []

    for adapt in ADAPTIVE_RULES:
        if not _matches(adapt["context_match"], ctx):
            continue
        for rule in adapt["rules"]:
            tid = rule["tool"]
            if rule.get("skip"):
                skipped_tools.append(tid)
                priority = [p for p in priority if p.get("tool") != tid]
                adjustments.append({
                    "tool": tid, "context": adapt["context_match"],
                    "action": "skip (venue-irrelevant)",
                })
                continue
            score = detailed_scores.get(tid)
            if score is None:
                continue
            if "score_below" in rule and score < rule["score_below"]:
                for p in priority:
                    if p.get("tool") == tid:
                        old_sev = p.get("severity")
                        new_sev = (
                            rule.get("downgrade_to")
                            or rule.get("upgrade_to")
                            or old_sev
                        )
                        if new_sev != old_sev:
                            p["severity"] = new_sev
                            p["adjusted_from"] = old_sev
                            adjustments.append({
                                "tool": tid,
                                "context": adapt["context_match"],
                                "action": f"{old_sev} → {new_sev}",
                                "score": score,
                            })
                        break

    sev_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 2, "UNKNOWN": 3, "LOW": 4}
    priority.sort(key=lambda x: (sev_rank.get(x.get("severity"), 99),
                                   -1 * (x.get("score") or 0)))

    n_c = sum(1 for p in priority if p.get("severity") == "CRITICAL")
    n_h = sum(1 for p in priority if p.get("severity") == "HIGH")

    comments = []
    if not venue:
        comments.append(
            "venue 未指定。固定閾値 health_report と同じ。"
            "venue 指定 (学会発表_短/長 / 招待講演 / 研究室セミナー / 授業) で adjust。"
        )
    else:
        comments.append(f"Venue: {venue}")
    if adjustments:
        comments.append(f"Adjustments: {len(adjustments)} 件")
    if skipped_tools:
        comments.append(f"Skipped: {skipped_tools}")
    if n_c >= 1:
        comments.append(f"Adjusted CRITICAL: {n_c}")

    return {
        "context": ctx,
        "adjusted_priority_issues": priority,
        "adjustments_applied": adjustments,
        "skipped_tools_due_to_context": skipped_tools,
        "n_critical_adjusted": n_c,
        "n_high_adjusted": n_h,
        "overall_score": health.get("overall_score"),
        "tools_run": health.get("tools_run"),
        "tools_skipped": health.get("tools_skipped"),
        "comments": comments,
        "hint": "presentation T27 は grant T51 / paper T20 の pptx 版。5 venue group。",
        "source": "Intelligence Layer S6 (presentation 版)。(v0.25.0)",
    }
