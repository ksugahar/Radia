"""Presentation T24 (Intelligence Layer S2): next-5-actions
(presentation_next_5_actions).
"""

from __future__ import annotations


TOOL_ESTIMATES: dict[str, dict] = {
    "T1": {"effort": 2, "impact": 5, "name": "opening hook",
           "rewrite_example": "冒頭 2 slide で Hook (問題提起 + 数字) を置く。"},
    "T2": {"effort": 2, "impact": 5, "name": "take-home strength",
           "rewrite_example": "末尾 slide に 1 文の take-home message を明示。"},
    "T3": {"effort": 1, "impact": 3, "name": "pie/3D chart",
           "rewrite_example": "3D pie → 2D bar に変更 (Knaflic Ch.3)。"},
    "T4": {"effort": 1, "impact": 2, "name": "logo consistency",
           "rewrite_example": "logo を全 slide に配置 or 全 slide から削除、統一。"},
    "T5": {"effort": 1, "impact": 2, "name": "progress indicator",
           "rewrite_example": "ページ番号 or progress bar で進行可視化。"},
    "T6": {"effort": 2, "impact": 3, "name": "visual/text ratio",
           "rewrite_example": "text-heavy slide に図を追加、1 slide visual 40% 以上。"},
    "T7": {"effort": 3, "impact": 3, "name": "speaker notes",
           "rewrite_example": "speaker note を全 slide に、発表時の台本として機能。"},
    "T8": {"effort": 2, "impact": 2, "name": "font consistency",
           "rewrite_example": "font 種別を 1-2 種に統一、size 階層 (title ≥32pt, body ≥24pt)。"},
    "T9": {"effort": 2, "impact": 3, "name": "arrow usage",
           "rewrite_example": "矢印は因果 / 時間 / 階層のどれかに統一、装飾矢印削減。"},
    "T10": {"effort": 1, "impact": 2, "name": "underline",
            "rewrite_example": "pptx 内 underline を 3-6/slide に削減、太字で代替。"},
    "T11": {"effort": 3, "impact": 4, "name": "slide density",
            "rewrite_example": "density 均等化、過密 slide を 2 slide に分割。"},
    "T13": {"effort": 3, "impact": 5, "name": "title outline coherence",
            "rewrite_example": "全 title を並べて目次化し、対象＋観点を短く具体化。"},
    "T14": {"effort": 2, "impact": 4, "name": "title-body alignment",
            "rewrite_example": "title を body に対応する対象＋観点の短い具体語句へ変更。"},
    "T15": {"effort": 3, "impact": 5, "name": "single message semantic",
            "rewrite_example": "1 slide に figure + 段落 + bullet が混在する場合は 2 slide に分割。"},
    "T16": {"effort": 3, "impact": 4, "name": "western text-heavy flag",
            "rewrite_example": "1 slide ≥100 字の欧米式を日本理系基準 ≤50 字に圧縮。"},
    "T17": {"effort": 3, "impact": 4, "name": "mini-IMRAD structure",
            "rewrite_example": "title → bg → problem → methods → results → discussion → conclusion 7 phase に再編。"},
    "T18": {"effort": 2, "impact": 4, "name": "理系 minimalism",
            "rewrite_example": "各 slide に figure + 数値 + 出典 + axis 単位の 5 要素を確保。"},
    "T19": {"effort": 2, "impact": 3, "name": "chart simplification",
            "rewrite_example": "3D / 凡例 ≥7 / 1 slide 3+ chart を修正。"},
    "T20": {"effort": 2, "impact": 4, "name": "equation compliance",
            "rewrite_example": "数式 slide に 記号定義 + 物理意味 + 次元単位 + 出典 の 4 要素。"},
    "T21": {"effort": 2, "impact": 4, "name": "figure compliance",
            "rewrite_example": "figure slide に 軸 + 単位 + 凡例 + caption + 出典 の 5 要素。"},
    "T22": {"effort": 3, "impact": 5, "name": "results statistical evidence",
            "rewrite_example": "Results slide に error bar + n + 統計検定 + effect size の 4 要素。"},
}


def _sev_val(s): return {"CRITICAL": 5, "HIGH": 4, "MEDIUM": 3, "LOW": 1}.get(s, 1)


def presentation_next_5_actions(pptx_path: str, skip: str = "") -> dict:
    try:
        from .T12 import presentation_health_report
    except Exception as e:
        return {"error": f"health_report import failed: {e}",
                "source": "presentation T24 next 5"}

    health = presentation_health_report(pptx_path, skip=skip)
    priority = health.get("priority_issues", [])

    actions = []
    for issue in priority:
        tid = issue.get("tool", "")
        sev = issue.get("severity", "LOW")
        score = issue.get("score") or 0
        est = TOOL_ESTIMATES.get(tid, {})
        effort, impact = est.get("effort", 3), est.get("impact", 3)
        name = est.get("name", issue.get("name", tid))
        rewrite = est.get("rewrite_example", f"{tid} の comments を確認。")
        ps = _sev_val(sev) * impact / max(effort, 1) * max(10 - score, 1)
        actions.append({
            "tool_id": tid, "name": name, "severity": sev,
            "current_score": round(score, 1) if score else None,
            "estimated_effort": effort, "estimated_impact": impact,
            "priority_score": round(ps, 1),
            "rewrite_example": rewrite,
            "comments_from_health_report": issue.get("comments", [])[:2],
        })

    actions.sort(key=lambda a: -a["priority_score"])
    top_5 = actions[:5]

    comments = []
    if not top_5:
        comments.append("Priority issue 無し。音読で最終判定。")
    else:
        comments.append(f"Top {len(top_5)} actions:")
        for i, a in enumerate(top_5, 1):
            comments.append(
                f"  {i}. [{a['severity']}] {a['name']} "
                f"(score {a['current_score']}, eff {a['estimated_effort']}, "
                f"imp {a['estimated_impact']}) → {a['rewrite_example'][:70]}…"
            )

    return {
        "top_5_actions": top_5,
        "all_action_count": len(actions),
        "overall_score": health.get("overall_score"),
        "overall_severity": health.get("overall_severity"),
        "comments": comments,
        "hint": "presentation T24 は grant T47 / paper T17 の pptx 版。",
        "source": "Intelligence Layer S2 (presentation 版)。(v0.25.0)",
    }
