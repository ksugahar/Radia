"""Paper T17 (Intelligence Layer S2): Priority next-5-actions generator
(paper_writing_next_5_actions).

paper 版 next-5-actions。paper T8 health_report の priority_issues を
severity * impact / effort * gap で再 sort、top-5 に rewrite example を添える。
"""

from __future__ import annotations


TOOL_ESTIMATES: dict[str, dict] = {
    "T1": {
        "effort": 2, "impact": 5, "name": "abstract strength",
        "rewrite_example": (
            "Abstract を Bg→Methods→Results→Contribution の 4 要素構成に再編。"
            "contribution は必ず末尾 1-2 文で明示。background ≤25%。"
        ),
    },
    "T2": {
        "effort": 2, "impact": 5, "name": "contribution clarity",
        "rewrite_example": (
            "Introduction 末尾に「本論文の contribution は (1)...(2)...(3)...」"
            "の 3 点列挙を明示。各 contribution に数値裏付けを添える。"
        ),
    },
    "T3": {
        "effort": 3, "impact": 4, "name": "claim quantification",
        "rewrite_example": (
            "各主張に数値 (%, 倍率, p, CI, effect size) を添える。"
            "「improves」「better」だけの主張は補強必要。"
        ),
    },
    "T4": {
        "effort": 1, "impact": 4, "name": "limitation statement",
        "rewrite_example": (
            "Discussion 中段に「本研究の limitation は...」1 段落必須。"
            "書かないと reviewer が勝手に掘る (Wallwork §9.12)。"
        ),
    },
    "T5": {
        "effort": 3, "impact": 3, "name": "related work density",
        "rewrite_example": (
            "Related Work の citation density を見直し、直近 5 年の "
            "論文を 3-5 本追加。自己引用比率 ≤30%。"
        ),
    },
    "T6": {
        "effort": 2, "impact": 3, "name": "figure referencing coverage",
        "rewrite_example": (
            "全 \\label が本文で \\ref されるよう確認。"
            "unreferenced figure は削除 or text で言及追加。"
        ),
    },
    "T7": {
        "effort": 2, "impact": 3, "name": "given-new ordering",
        "rewrite_example": (
            "段落間で既知→新情報の流れを確認 (Wallwork §3.4-3.6)。"
            "新 topic 段落の第 1 文は前段落との繋ぎ必須。"
        ),
    },
    "T9": {
        "effort": 3, "impact": 5, "name": "reviewer #2 triggers",
        "rewrite_example": (
            "reviewer_2_trigger_summary の 4-tier から Top/High を対応。"
            "strong adjective 減、自己引用比率確認、統計報告 compliance 強化。"
        ),
    },
    "T10": {
        "effort": 2, "impact": 5, "name": "Title/Abstract/Conclusion triangle",
        "rewrite_example": (
            "3 セクションを順に音読、keyword overlap を確認。"
            "Title の数値が Abstract/Conclusion に裏付けあるか必ず check。"
        ),
    },
    "T11": {
        "effort": 2, "impact": 4, "name": "reproducibility (6 軸)",
        "rewrite_example": (
            "code/data availability + preregistration + methods replicability "
            "+ seed + COI の 6 軸を methods / acknowledgments に articulate。"
        ),
    },
    "T12": {
        "effort": 3, "impact": 5, "name": "statistical reporting",
        "rewrite_example": (
            "各 p 値に effect size + 95% CI + sample size n の 4 要素 seq 添付。"
            "p < 0.05 単独は desk reject リスク (APA 7th / CONSORT 違反)。"
        ),
    },
    "T13": {
        "effort": 3, "impact": 3, "name": "citation health",
        "rewrite_example": (
            "recency (median age ≤5年) / 地理分散 (3+ region) / "
            "著者集中 (top1 <20%) / self-citation (<30%) の 4 軸 green 化。"
        ),
    },
    "T14": {
        "effort": 3, "impact": 4, "name": "discussion 4 要素",
        "rewrite_example": (
            "Discussion に interpretation + limitations + future_work + "
            "generalization の 4 要素段落を配置。Results 再記述で終わらせない。"
        ),
    },
    "T15": {
        "effort": 2, "impact": 3, "name": "journal fit",
        "rewrite_example": (
            "Target journal の Aims & Scope と論文 keywords の overlap 確認。"
            "score < 4 なら別 journal 検討 or cover letter で橋渡し。"
        ),
    },
}


def _severity_value(sev: str) -> int:
    return {"CRITICAL": 5, "HIGH": 4, "MEDIUM": 3, "LOW": 1}.get(sev, 1)


def paper_writing_next_5_actions(
    tex_or_text: str,
    bib: str = "",
    abstract: str = "",
    author_last_names: str = "",
    skip: str = "",
) -> dict:
    """paper health_report の priority_issues を impact / effort で再 sort、
    top-5 に rewrite example 付きで返す。"""
    if tex_or_text is None:
        tex_or_text = ""

    try:
        from .T8 import paper_writing_health_report
    except Exception as e:
        return {"error": f"health_report import failed: {e}",
                "source": "paper T17 next 5 actions"}

    health = paper_writing_health_report(
        tex_or_text, bib=bib, abstract=abstract,
        author_last_names=author_last_names, skip=skip,
    )
    priority_issues = health.get("priority_issues", [])

    actions: list[dict] = []
    for issue in priority_issues:
        tid = issue.get("tool", "")
        sev = issue.get("severity", "LOW")
        score = issue.get("score") or 0
        est = TOOL_ESTIMATES.get(tid, {})
        effort = est.get("effort", 3)
        impact = est.get("impact", 3)
        name = est.get("name", issue.get("name", tid))
        rewrite = est.get("rewrite_example",
                           f"{tid} の comments を確認し対応。")
        priority_score = (
            _severity_value(sev) * impact / max(effort, 1)
            * max(10 - score, 1)
        )
        actions.append({
            "tool_id": tid,
            "name": name,
            "severity": sev,
            "current_score": round(score, 1) if score else None,
            "estimated_effort": effort,
            "estimated_impact": impact,
            "priority_score": round(priority_score, 1),
            "rewrite_example": rewrite,
            "comments_from_health_report": issue.get("comments", [])[:2],
        })

    actions.sort(key=lambda a: -a["priority_score"])
    top_5 = actions[:5]

    comments: list[str] = []
    if not top_5:
        comments.append("Priority issue 無し (全 tool LOW)。"
                         "最終パスは音読で著者の声で判定。")
    else:
        comments.append(
            f"Top {len(top_5)} actions (priority = severity * impact / "
            "effort * gap):"
        )
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
        "hint": (
            "paper T17 は grant T47 の paper 版。"
            "priority_score で ROI 最大の action から順に対応可能。"
            "rewrite_example は template、context に応じて caller が調整。"
        ),
        "source": (
            "Intelligence Layer S2 (paper 版)。T8 health_report の "
            "priority_issues を rewrite example 付きで再 sort。"
            "(paper-writing 2026-04 採用、v0.24.0)。"
        ),
    }
