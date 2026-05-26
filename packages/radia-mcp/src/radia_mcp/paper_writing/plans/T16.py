"""Paper T16 (Intelligence Layer S1): Cross-tool root cause synthesis
(paper_writing_root_cause_diagnosis).

paper_writing 版の grant T46 に相当。複数 tool の low score を correlate して
典型的な「論文の根本原因」を診断する。

Embedded root cause patterns (paper 固有 8 種):
  1. unsupported_claim      — claim quantification + stats + limitation 弱
  2. out_of_date_literature — related_work_density + citation health 弱
  3. triangle_misalignment  — T10 triangle + given-new 弱
  4. imrad_collapse         — IMRAD balance + discussion structure 弱
  5. reviewer_2_trigger     — reviewer_2_trigger + strong_adjective 多
  6. reproducibility_drift  — T11 reproducibility + stats 弱
  7. contribution_fog       — T2 contribution + abstract_strength 弱
  8. journal_mismatch       — journal_fit + abstract_background_ratio 弱
"""

from __future__ import annotations


ROOT_CAUSE_PATTERNS: list[dict] = [
    {
        "name": "unsupported_claim",
        "label": "主張が裏付け不足 (reviewer #2 筆頭トリガ)",
        "tool_signals": ["T3", "T12", "T4"],
        "threshold_avg_below": 6.0,
        "remediation": (
            "claim に対し統計/数値裏付け (T12 p+effect+CI+n) を必ず添え、"
            "limitation (T4) を 1 段落で明示、claim quantification (T3) "
            "で各主張に defensible な数値を付ける。"
        ),
    },
    {
        "name": "out_of_date_literature",
        "label": "文献が古い/狭い (literature review 不足印象)",
        "tool_signals": ["T5", "T13"],
        "threshold_avg_below": 6.0,
        "remediation": (
            "related_work_density (T5) を再確認、citation_health_4_axes (T13) で "
            "recency (median age ≤5年) / 地理分散 / 著者集中 / 自己引用 を "
            "全軸 green に。直近 5 年の論文を 3-5 本追加推奨。"
        ),
    },
    {
        "name": "triangle_misalignment",
        "label": "Title/Abstract/Conclusion が別の物語を語っている",
        "tool_signals": ["T10", "T7"],
        "threshold_avg_below": 6.0,
        "remediation": (
            "title_abstract_conclusion_triangle (T10) で 3 辺の overlap 確認、"
            "given_new_ordering (T7) で段落間の繋ぎを整える。"
            "title 過大広告 / conclusion 段差は信頼の一発勝負。"
        ),
    },
    {
        "name": "imrad_collapse",
        "label": "IMRAD 構造が崩れている (特に Discussion 薄)",
        "tool_signals": ["T14"],
        "threshold_avg_below": 5.0,
        "remediation": (
            "discussion_structure_4_elements (T14) で "
            "interpretation/limitations/future_work/generalization の 4 要素を "
            "articulate。Discussion を Results 再記述で終わらせない。"
        ),
    },
    {
        "name": "reviewer_2_trigger",
        "label": "Reviewer #2 が引き寄せる weakness 多発",
        "tool_signals": ["T9"],
        "threshold_avg_below": 5.0,
        "remediation": (
            "reviewer_2_trigger_summary (T9) の 4-tier classification を "
            "読み、Top-tier (致命的) と High-tier から順に対応。"
            "strong_adjective / unsupported claim / 自己引用過剰 等。"
        ),
    },
    {
        "name": "reproducibility_drift",
        "label": "再現性 compliance 不足 (実験系 reviewer 警戒)",
        "tool_signals": ["T11", "T12"],
        "threshold_avg_below": 6.0,
        "remediation": (
            "reproducibility_open_science_check (T11) 6 軸と "
            "statistical_reporting_compliance (T12) 4 要素で journal 基準を "
            "充足。code/data availability + preregistration + seed + COI。"
        ),
    },
    {
        "name": "contribution_fog",
        "label": "Contribution が曖昧 (「何が新しいか」分からない)",
        "tool_signals": ["T1", "T2"],
        "threshold_avg_below": 6.0,
        "remediation": (
            "contribution_clarity_score (T2) で contribution 箇所を明示、"
            "abstract_strength (T1) で abstract 後半に contribution 1-2 文 "
            "必須。reviewer は 10 秒で contribution を取る。"
        ),
    },
    {
        "name": "journal_mismatch",
        "label": "Target journal との scope 不一致 (desk reject リスク)",
        "tool_signals": ["T15"],
        "threshold_avg_below": 5.0,
        "remediation": (
            "journal_fit_assessment (T15) で aims & scope overlap を確認、"
            "abstract_background_ratio で導入部が journal aims に合うか check。"
            "別 journal 検討、または cover letter で scope 橋渡しを 1 段落。"
        ),
    },
]


def paper_writing_root_cause_diagnosis(
    tex_or_text: str,
    bib: str = "",
    abstract: str = "",
    author_last_names: str = "",
    skip: str = "",
) -> dict:
    """health_report 結果を横断 pattern matching し、論文の根本原因を診断。

    Args:
        tex_or_text / bib / abstract / author_last_names / skip:
            paper T8 health_report にそのまま passthrough

    Returns:
        dict: detected_root_causes / health_report_summary / comments
    """
    if tex_or_text is None:
        tex_or_text = ""

    try:
        from .T8 import paper_writing_health_report
    except Exception as e:
        return {
            "error": f"health_report import failed: {e}",
            "source": "paper T16 root cause diagnosis",
        }

    health = paper_writing_health_report(
        tex_or_text,
        bib=bib,
        abstract=abstract,
        author_last_names=author_last_names,
        skip=skip,
    )
    detailed_scores: dict = health.get("detailed_scores", {})

    detected_root_causes: list[dict] = []
    for pattern in ROOT_CAUSE_PATTERNS:
        signals = pattern["tool_signals"]
        relevant_scores = [
            detailed_scores[t] for t in signals
            if t in detailed_scores and detailed_scores[t] is not None
        ]
        if not relevant_scores:
            continue
        avg = sum(relevant_scores) / len(relevant_scores)
        if avg <= pattern["threshold_avg_below"]:
            severity = (
                "CRITICAL" if avg <= 3 else
                "HIGH" if avg <= 5 else
                "MEDIUM"
            )
            detected_root_causes.append({
                "pattern_name": pattern["name"],
                "label": pattern["label"],
                "severity": severity,
                "tool_signals": signals,
                "tool_scores": {t: detailed_scores.get(t) for t in signals},
                "avg_score": round(avg, 2),
                "remediation": pattern["remediation"],
            })

    sev_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 2}
    detected_root_causes.sort(
        key=lambda x: (sev_rank.get(x["severity"], 9), x["avg_score"])
    )

    comments: list[str] = []
    if not detected_root_causes:
        comments.append(
            "Root cause pattern 検出無し。個別 tool の low score があれば "
            "それぞれ対応。health_report priority_issues を併せて確認。"
        )
    else:
        n_critical = sum(1 for r in detected_root_causes if r["severity"] == "CRITICAL")
        n_high = sum(1 for r in detected_root_causes if r["severity"] == "HIGH")
        comments.append(
            f"Root cause {len(detected_root_causes)} 件検出 "
            f"(CRITICAL {n_critical}, HIGH {n_high})"
        )
        for r in detected_root_causes[:3]:
            comments.append(
                f"[{r['severity']}] {r['label']}: avg {r['avg_score']} "
                f"(tools: {r['tool_signals']}) → {r['remediation'][:70]}…"
            )

    return {
        "detected_root_causes": detected_root_causes,
        "n_root_causes": len(detected_root_causes),
        "health_report_overall_score": health.get("overall_score"),
        "health_report_overall_severity": health.get("overall_severity"),
        "all_pattern_definitions": [
            {"name": p["name"], "label": p["label"],
             "tool_signals": p["tool_signals"]}
            for p in ROOT_CAUSE_PATTERNS
        ],
        "comments": comments,
        "hint": (
            "paper T16 は grant T46 の paper 版。8 典型 pattern で "
            "論文の構造的問題を診断。T8 health_report の補完として使う。"
        ),
        "source": (
            "Intelligence Layer S1 (paper 版)。paper T1-T15 の score を "
            "横断 synthesis し、paper-specific な典型 root cause を診断。"
            "(paper-writing 2026-04 採用、v0.24.0)。"
        ),
    }
