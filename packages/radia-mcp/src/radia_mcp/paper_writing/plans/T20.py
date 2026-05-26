"""Paper T20 (Intelligence Layer S6): Adaptive health report (paper 用)
(paper_writing_adaptive_health_report).

target journal (IEEE Trans / Nature / Science / IEEJ 等) や publication phase
(first submission / revision / camera-ready) で severity 判定を adjust。
"""

from __future__ import annotations


ADAPTIVE_RULES: list[dict] = [
    # Nature / Science 系: novelty + broad impact 最重視
    {
        "context_match": {"journal_tier": ("top_tier", "nature", "science")},
        "rules": [
            {"tool": "T1", "score_below": 8, "upgrade_to": "CRITICAL"},
            {"tool": "T2", "score_below": 8, "upgrade_to": "CRITICAL"},
            {"tool": "T10", "score_below": 7, "upgrade_to": "HIGH"},
            {"tool": "T11", "score_below": 6, "upgrade_to": "HIGH"},
        ],
    },
    # IEEE Trans 系: technical depth + 統計報告
    {
        "context_match": {"journal_tier": ("ieee", "ieee_trans")},
        "rules": [
            {"tool": "T3", "score_below": 7, "upgrade_to": "HIGH"},
            {"tool": "T12", "score_below": 7, "upgrade_to": "HIGH"},
            {"tool": "T11", "score_below": 6, "upgrade_to": "HIGH"},
        ],
    },
    # IEEJ / 国内論文誌: IMRAD + 日本語作文技術
    {
        "context_match": {"journal_tier": ("ieej", "domestic_jp")},
        "rules": [
            {"tool": "T10", "score_below": 6, "upgrade_to": "HIGH"},
            {"tool": "T7", "score_below": 6, "upgrade_to": "HIGH"},
            # journal fit は国内誌では緩め
            {"tool": "T15", "score_below": 5, "downgrade_to": "MEDIUM"},
        ],
    },
    # Revision phase: reviewer comment 対応が最優先
    {
        "context_match": {"phase": "revision"},
        "rules": [
            {"tool": "T9", "score_below": 8, "upgrade_to": "CRITICAL"},
            {"tool": "T4", "score_below": 7, "upgrade_to": "HIGH"},
        ],
    },
    # Camera-ready: proofread 系が重要
    {
        "context_match": {"phase": "camera_ready"},
        "rules": [
            # この段階で大幅 rewrite はしない、表記 / 形式中心
            {"tool": "T1", "skip": True},
            {"tool": "T2", "skip": True},
        ],
    },
    # Preprint (arXiv 等): 緩め
    {
        "context_match": {"phase": "preprint"},
        "rules": [
            {"tool": "T15", "skip": True},  # journal fit 不要
            {"tool": "T10", "score_below": 5, "downgrade_to": "MEDIUM"},
        ],
    },
]


def _matches_context(rule_context: dict, ctx: dict) -> bool:
    for k, expected in rule_context.items():
        actual = ctx.get(k)
        if not actual:
            return False
        actual_lc = actual.lower() if isinstance(actual, str) else actual
        if isinstance(expected, tuple):
            return any(
                (actual_lc == e or (isinstance(actual, str) and e in actual_lc))
                for e in expected
            )
        return actual_lc == expected.lower()


def paper_writing_adaptive_health_report(
    tex_or_text: str,
    bib: str = "",
    abstract: str = "",
    author_last_names: str = "",
    journal_tier: str = "",
    phase: str = "",
    skip: str = "",
) -> dict:
    """paper T8 health_report の severity 判定を context で adjust。

    Args:
        tex_or_text / bib / abstract / author_last_names / skip: passthrough
        journal_tier: 'top_tier' / 'nature' / 'science' / 'ieee' /
                      'ieee_trans' / 'ieej' / 'domestic_jp' / ''
        phase: 'first_submission' / 'revision' / 'camera_ready' /
               'preprint' / ''
    """
    if tex_or_text is None:
        tex_or_text = ""

    try:
        from .T8 import paper_writing_health_report
    except Exception as e:
        return {"error": f"health_report import failed: {e}",
                "source": "paper T20 adaptive health report"}

    health = paper_writing_health_report(
        tex_or_text, bib=bib, abstract=abstract,
        author_last_names=author_last_names, skip=skip,
    )

    ctx = {"journal_tier": journal_tier, "phase": phase}
    detailed_scores = health.get("detailed_scores", {})
    priority_issues = list(health.get("priority_issues", []))
    adjustments_applied: list[dict] = []
    skipped_tools: list[str] = []

    for adapt in ADAPTIVE_RULES:
        if not _matches_context(adapt["context_match"], ctx):
            continue
        for rule in adapt["rules"]:
            tid = rule["tool"]
            if rule.get("skip"):
                skipped_tools.append(tid)
                priority_issues = [
                    p for p in priority_issues if p.get("tool") != tid
                ]
                adjustments_applied.append({
                    "tool": tid,
                    "context": adapt["context_match"],
                    "action": "skip (phase-irrelevant)",
                })
                continue
            score = detailed_scores.get(tid)
            if score is None:
                continue
            if "score_below" in rule and score < rule["score_below"]:
                for p in priority_issues:
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
                            adjustments_applied.append({
                                "tool": tid,
                                "context": adapt["context_match"],
                                "action": f"{old_sev} → {new_sev}",
                                "score": score,
                            })
                        break

    sev_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 2, "UNKNOWN": 3, "LOW": 4}
    priority_issues.sort(
        key=lambda x: (sev_rank.get(x.get("severity"), 99),
                       -1 * (x.get("score") or 0))
    )

    n_critical = sum(1 for p in priority_issues if p.get("severity") == "CRITICAL")
    n_high = sum(1 for p in priority_issues if p.get("severity") == "HIGH")

    comments: list[str] = []
    if not journal_tier and not phase:
        comments.append(
            "journal_tier / phase が未指定。固定閾値 health_report と同じ。"
            "context 指定で adjust が走る。"
        )
    else:
        ctx_str = ", ".join(f"{k}={v}" for k, v in ctx.items() if v)
        comments.append(f"Context: {ctx_str}")

    if adjustments_applied:
        comments.append(f"Adjustments: {len(adjustments_applied)} 件")
        for adj in adjustments_applied[:5]:
            comments.append(f"  • {adj['tool']}: {adj['action']}")

    if skipped_tools:
        comments.append(f"Skipped: {skipped_tools}")

    if n_critical >= 1:
        comments.append(f"Adjusted CRITICAL: {n_critical} 件")
    elif n_high >= 1:
        comments.append(f"Adjusted HIGH: {n_high} 件")

    return {
        "context": ctx,
        "adjusted_priority_issues": priority_issues,
        "adjustments_applied": adjustments_applied,
        "skipped_tools_due_to_context": skipped_tools,
        "n_critical_adjusted": n_critical,
        "n_high_adjusted": n_high,
        "overall_score": health.get("overall_score"),
        "tools_run": health.get("tools_run"),
        "tools_skipped": health.get("tools_skipped"),
        "comments": comments,
        "hint": (
            "paper T20 は grant T51 の paper 版。"
            "journal tier (top/IEEE/IEEJ) や phase (first/revision/CR/preprint) "
            "で severity を adjust。reviewer 視点の重み付けを反映。"
        ),
        "source": (
            "Intelligence Layer S6 (paper 版)。paper 用 6 group × ~15 rule "
            "を embed、context-aware に severity 判定。"
            "(paper-writing 2026-04 採用、v0.24.0)。"
        ),
    }
