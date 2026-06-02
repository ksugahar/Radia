"""T12: presentation_health_report — Plan B 診断の統合レポート。

Plan B Tier 1 (T1-T5) + Tier 2 (T6-T11) を束ね、修正優先度付きで summary
を返す。発表準備レビュー時の「入口」として使う。

grant_writing T16 / paper_writing T8 と同じ return shape。
"""
from __future__ import annotations

import pathlib

from .T1 import presentation_opening_hook_strength
from .T2 import presentation_takehome_strength
from .T3 import presentation_check_pie_3d_charts
from .T4 import presentation_check_logo_on_every_slide
from .T5 import presentation_check_progress_indicator
from .T6 import presentation_visual_text_ratio_score
from .T7 import presentation_speaker_note_ratio
from .T8 import presentation_font_consistency
from .T9 import presentation_arrow_usage
from .T10 import presentation_check_underline_in_pptx
from .T11 import presentation_slide_density_balance
# Tier 3 (v0.21.0) — pptx_path 単独で動く tool (T13-T22 全て) を統合
from .T13 import presentation_slide_titles_outline_coherence
from .T14 import presentation_title_body_alignment_check
from .T15 import presentation_single_message_per_slide_semantic
from .T16 import presentation_text_density_per_slide_western_style
from .T17 import presentation_mini_imrad_structure_check
from .T18 import presentation_rikei_minimalism_score
from .T19 import presentation_chart_simplification_check
# Tier 4 (v0.21.0) — 理系特化
from .T20 import presentation_equation_slide_compliance
from .T21 import presentation_figure_slide_compliance
from .T22 import presentation_results_slide_statistical_evidence


_HINT = (
    "pptx path を 1 つ与えるだけで発表準備の多面性を診断。"
    "score 最適化ゲームに陥らず、priority_issues の CRITICAL/HIGH から"
    "手を付ける。宮野 S1-S15 と skill.md と照合。"
)

_SOURCE = (
    "Plan B 統合レポート v0.21.0 — "
    "T1-T11 (Hook/Takehome/Chart/Logo/Progress/Visual/Notes/Font/Arrow/"
    "Underline/Density) + T13-T19 (titles outline / title-body / "
    "single message / 欧米 text / mini-IMRAD / 理系 minimalism / chart) + "
    "T20-T22 (理系特化: equation / figure / results stats) の 21 観点"
)


def _severity_from_score(score) -> str:
    if score is None:
        return "UNKNOWN"
    if score < 4:
        return "CRITICAL"
    if score < 6:
        return "HIGH"
    if score < 8:
        return "MEDIUM"
    return "LOW"


def presentation_health_report(
    pptx_path: str,
    skip: str = "",
) -> dict:
    """presentation Plan B の全 T1-T11 を束ねた統合レポート。

    grant_writing T16 / paper_writing T8 と同じ設計。発表スライドの全体像を
    1 コールで把握するための入口 tool。各 subscore に severity
    (CRITICAL/HIGH/MEDIUM/LOW) を付け、修正優先度付きの priority_issues
    リストを返す。

    Args:
        pptx_path: 対象 .pptx のパス。
        skip: カンマ区切りで除外する T 番号 (例: "T4,T8" でロゴ/フォント診断を
              スキップ)。該当観点を既に手動確認済の場合などに有用。

    Returns:
        dict:
          - overall_score: 0-10 の粗い総合スコア
          - overall_severity: CRITICAL/HIGH/MEDIUM/LOW
          - summary_comment: 1-2 文の overall 判定
          - detailed_scores: 各 tool の score
          - priority_issues: severity 順にソートされた issue list
          - tools_run / tools_skipped
          - hint / source
    """
    skip_set = {s.strip().upper() for s in skip.split(",") if s.strip()}

    # ---- 前段: pptx 自体が存在しない場合の早期 return ----
    p = pathlib.Path(pptx_path)
    if not p.exists():
        return {
            "overall_score": 0.0,
            "score_max": 10,
            "overall_severity": "CRITICAL",
            "summary_comment": (
                f"pptx が見つからない: {pptx_path}。"
                "path を確認して再実行すること。"
            ),
            "detailed_scores": {},
            "priority_issues": [{
                "tool": "T12",
                "name": "health_report",
                "severity": "CRITICAL",
                "reason": f"file not found: {pptx_path}",
                "comments": [],
            }],
            "tools_run": [],
            "tools_skipped": sorted(skip_set),
            "hint": _HINT,
            "source": _SOURCE,
        }

    tools_to_run: list[tuple[str, str, callable]] = [
        # Tier 1
        ("T1", "opening_hook_strength", presentation_opening_hook_strength),
        ("T2", "takehome_strength", presentation_takehome_strength),
        ("T3", "check_pie_3d_charts", presentation_check_pie_3d_charts),
        ("T4", "check_logo_on_every_slide", presentation_check_logo_on_every_slide),
        ("T5", "check_progress_indicator", presentation_check_progress_indicator),
        # Tier 2
        ("T6", "visual_text_ratio_score", presentation_visual_text_ratio_score),
        ("T7", "speaker_note_ratio", presentation_speaker_note_ratio),
        ("T8", "font_consistency", presentation_font_consistency),
        ("T9", "arrow_usage", presentation_arrow_usage),
        ("T10", "check_underline_in_pptx", presentation_check_underline_in_pptx),
        ("T11", "slide_density_balance", presentation_slide_density_balance),
        # Tier 3 (v0.21.0)
        ("T13", "slide_titles_outline_coherence",
         presentation_slide_titles_outline_coherence),
        ("T14", "title_body_alignment_check",
         presentation_title_body_alignment_check),
        ("T15", "single_message_per_slide_semantic",
         presentation_single_message_per_slide_semantic),
        ("T16", "text_density_per_slide_western_style",
         presentation_text_density_per_slide_western_style),
        ("T17", "mini_imrad_structure_check",
         presentation_mini_imrad_structure_check),
        ("T18", "rikei_minimalism_score",
         presentation_rikei_minimalism_score),
        ("T19", "chart_simplification_check",
         presentation_chart_simplification_check),
        # Tier 4 (理系特化)
        ("T20", "equation_slide_compliance",
         presentation_equation_slide_compliance),
        ("T21", "figure_slide_compliance",
         presentation_figure_slide_compliance),
        ("T22", "results_slide_statistical_evidence",
         presentation_results_slide_statistical_evidence),
    ]

    detailed_scores: dict[str, float] = {}
    detailed_results: dict[str, dict] = {}
    priority_issues: list[dict] = []

    for tid, name, func in tools_to_run:
        if tid in skip_set:
            continue
        try:
            r = func(pptx_path)
        except Exception as e:
            priority_issues.append({
                "tool": tid,
                "name": name,
                "severity": "CRITICAL",
                "reason": f"tool crashed: {e}",
                "comments": [],
            })
            continue

        score = r.get("score")
        detailed_scores[tid] = score
        detailed_results[tid] = r
        severity = _severity_from_score(score)

        # LOW は priority_issues に入れない (修正不要な acknowledgment)
        if severity != "LOW":
            priority_issues.append({
                "tool": tid,
                "name": name,
                "severity": severity,
                "score": score,
                "comments": r.get("comments", [])[:3],
            })

    # severity 順にソート (CRITICAL > HIGH > MEDIUM > UNKNOWN > LOW)
    sev_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 2, "UNKNOWN": 3, "LOW": 4}
    priority_issues.sort(key=lambda x: (sev_rank.get(x["severity"], 99),
                                         -1 * (x.get("score") or 0)))

    # overall = None を除いた単純平均
    valid_scores = [s for s in detailed_scores.values() if s is not None]
    overall_score = (round(sum(valid_scores) / len(valid_scores), 1)
                     if valid_scores else 0.0)
    overall_severity = _severity_from_score(overall_score)

    # summary comment
    n_critical = sum(1 for i in priority_issues if i["severity"] == "CRITICAL")
    n_high = sum(1 for i in priority_issues if i["severity"] == "HIGH")
    if n_critical >= 2:
        summary = (f"CRITICAL 課題が {n_critical} 件。"
                   "Hook/Takehome/視覚比のいずれかが根本的に欠けている可能性。"
                   "priority_issues の CRITICAL から順に対応。")
    elif n_critical == 1:
        summary = (f"CRITICAL 課題が 1 件 ({priority_issues[0]['name']})。"
                   "まずここを直すと発表の印象が大きく変わる。")
    elif n_high >= 2:
        summary = (f"HIGH 課題が {n_high} 件。致命傷は無いが、"
                   "修正で理解度が上がる余地がある。")
    elif n_high == 1:
        summary = (f"HIGH 課題が 1 件 ({priority_issues[0]['name']})。"
                   "それ以外は概ね整合。")
    else:
        summary = ("CRITICAL / HIGH 課題なし。発表スライドの骨格は整っている。"
                   "最終パスは tool を閉じて通しリハで最終判定。")

    return {
        "overall_score": overall_score,
        "score_max": 10,
        "overall_severity": overall_severity,
        "summary_comment": summary,
        "detailed_scores": detailed_scores,
        "priority_issues": priority_issues,
        "tools_run": list(detailed_scores.keys()),
        "tools_skipped": sorted(skip_set),
        "hint": _HINT,
        "source": _SOURCE,
    }
