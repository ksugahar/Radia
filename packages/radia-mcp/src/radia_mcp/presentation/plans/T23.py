"""Presentation T23 (Intelligence Layer S1): Root cause synthesis
(presentation_root_cause_diagnosis).
"""

from __future__ import annotations


ROOT_CAUSE_PATTERNS: list[dict] = [
    {
        "name": "story_collapse",
        "label": "ストーリーが通っていない (title 並べて目次として破綻)",
        "tool_signals": ["T13", "T17", "T15"],  # outline / IMRAD / single msg
        "threshold_avg_below": 6.0,
        "remediation": (
            "T13 slide_titles_outline_coherence でタイトル並びを確認、"
            "T17 mini_imrad_structure で 7 phase の順序整合、"
            "T15 single_message で 1 slide 1 主張を徹底。"
        ),
    },
    {
        "name": "crowded_slides",
        "label": "スライド詰め込みすぎ (reviewer 認知負荷 overflow)",
        "tool_signals": ["T11", "T15", "T16"],
        "threshold_avg_below": 6.0,
        "remediation": (
            "T11 slide_density_balance + T15 single_message + "
            "T16 欧米 text-heavy で density 制御。1 slide ≤50 字目標。"
        ),
    },
    {
        "name": "weak_opening_closing",
        "label": "Hook / Take-home が弱い (印象に残らない)",
        "tool_signals": ["T1", "T2"],
        "threshold_avg_below": 6.0,
        "remediation": (
            "T1 opening_hook_strength で冒頭 2 slide を強化、"
            "T2 takehome_strength で末尾の 1 文 message を明確化。"
        ),
    },
    {
        "name": "evidence_thin_slides",
        "label": "Evidence 薄い (理系プレゼンの core 欠落)",
        "tool_signals": ["T18", "T20", "T21", "T22"],
        "threshold_avg_below": 6.0,
        "remediation": (
            "T18 理系 minimalism で figure + 数値 + 出典、"
            "T20 equation + T21 figure + T22 results stats で "
            "各 slide 種別の compliance を充足。"
        ),
    },
    {
        "name": "title_body_drift",
        "label": "Title の対象・観点が body の内容と対応していない",
        "tool_signals": ["T14"],
        "threshold_avg_below": 5.0,
        "remediation": (
            "T14 title_body_alignment で title-body overlap 確認、"
            "対象＋観点の短い具体語句に書き換え。"
        ),
    },
    {
        "name": "chart_clutter",
        "label": "Chart 雑多 (3D/凡例多/clutter)",
        "tool_signals": ["T19", "T3"],
        "threshold_avg_below": 6.0,
        "remediation": (
            "T19 chart_simplification で 3D / 凡例 ≥7 / "
            "1 slide 3+ charts を修正、T3 pie_3d で円グラフ形状。"
        ),
    },
    {
        "name": "pace_mismatch",
        "label": "時間配分が崩れている",
        "tool_signals": ["T5"],
        "threshold_avg_below": 5.0,
        "remediation": (
            "T5 check_progress_indicator で進行表示、"
            "speaker_note と slide 数から time estimate 再計算。"
        ),
    },
]


def presentation_root_cause_diagnosis(pptx_path: str, skip: str = "") -> dict:
    """pptx の health_report を横断 pattern matching、根本原因診断。"""
    try:
        from .T12 import presentation_health_report
    except Exception as e:
        return {"error": f"health_report import failed: {e}",
                "source": "presentation T23 root cause"}

    health = presentation_health_report(pptx_path, skip=skip)
    detailed_scores = health.get("detailed_scores", {})

    detected: list[dict] = []
    for p in ROOT_CAUSE_PATTERNS:
        signals = p["tool_signals"]
        rel = [detailed_scores[t] for t in signals
               if t in detailed_scores and detailed_scores[t] is not None]
        if not rel:
            continue
        avg = sum(rel) / len(rel)
        if avg <= p["threshold_avg_below"]:
            sev = "CRITICAL" if avg <= 3 else "HIGH" if avg <= 5 else "MEDIUM"
            detected.append({
                "pattern_name": p["name"], "label": p["label"],
                "severity": sev, "tool_signals": signals,
                "tool_scores": {t: detailed_scores.get(t) for t in signals},
                "avg_score": round(avg, 2),
                "remediation": p["remediation"],
            })

    sev_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 2}
    detected.sort(key=lambda x: (sev_rank.get(x["severity"], 9), x["avg_score"]))

    comments = []
    if not detected:
        comments.append("Root cause 検出無し。個別 tool 確認。")
    else:
        comments.append(f"Root cause {len(detected)} 件検出")
        for r in detected[:3]:
            comments.append(
                f"[{r['severity']}] {r['label']}: avg {r['avg_score']} "
                f"→ {r['remediation'][:70]}…"
            )

    return {
        "detected_root_causes": detected,
        "n_root_causes": len(detected),
        "health_report_overall_score": health.get("overall_score"),
        "all_pattern_definitions": [
            {"name": p["name"], "label": p["label"],
             "tool_signals": p["tool_signals"]}
            for p in ROOT_CAUSE_PATTERNS
        ],
        "comments": comments,
        "hint": "presentation T23 は grant T46 / paper T16 の pptx 版。7 典型 pattern。",
        "source": "Intelligence Layer S1 (presentation 版)。(v0.25.0)",
    }
