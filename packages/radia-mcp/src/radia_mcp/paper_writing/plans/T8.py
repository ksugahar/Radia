"""T8: paper_writing_health_report — Plan B 統合レポート。

paper_writing Plan B T1-T7 を束ね、修正優先度付きの summary を返す。
論文ドラフトを 1 コールでレビューするための入口。

設計方針 (grant_writing T16 と整合):
- severity: CRITICAL (<4) / HIGH (<6) / MEDIUM (<8) / LOW (else)
- priority_issues は LOW を除外し、severity 順でソート
- bib / abstract が空なら自動で T5 / T1 を skip
- skip="T3,T5" 形式で手動スキップ可能
"""
from __future__ import annotations

from .T1 import paper_writing_abstract_strength
from .T2 import paper_writing_contribution_clarity_score
from .T3 import paper_writing_claim_quantification
from .T4 import paper_writing_limitation_statement_presence
from .T5 import paper_writing_related_work_density
from .T6 import paper_writing_figure_referencing_coverage
from .T7 import paper_writing_given_new_ordering
# Tier 3 (v0.19.0) — single-text 動作の tool のみ統合
# T12 stats / T13 citation health は他引数必須、T15 journal fit は aims 必須のため除外
from .T10 import paper_writing_title_abstract_conclusion_triangle
from .T11 import paper_writing_reproducibility_open_science_check
from .T14 import paper_writing_discussion_structure_4_elements


def _severity_from_score(score: float | None) -> str:
    if score is None:
        return "UNKNOWN"
    if score < 4:
        return "CRITICAL"
    if score < 6:
        return "HIGH"
    if score < 8:
        return "MEDIUM"
    return "LOW"


def paper_writing_health_report(
    tex_or_text: str,
    bib: str = "",
    abstract: str = "",
    author_last_names: str = "",
    skip: str = "",
) -> dict:
    """paper_writing Plan B の全 T1-T7 を束ねた統合レポート。

    grant_writing T16 と同じ shape。priority_issues (CRITICAL/HIGH/MEDIUM/LOW)、
    detailed_scores、summary_comment を返し、論文ドラフトを 1 コールで
    レビューするための入口として使う。

    Args:
        tex_or_text: 論文原稿 (LaTeX or プレーンテキスト)。
        bib: .bib ファイル内容。空なら T5 (related_work_density) を自動 skip。
        abstract: Abstract 本文。空なら T1 (abstract_strength) を自動 skip。
        author_last_names: 自分+共著者の姓 (T5 用、カンマ区切り)。
        skip: カンマ区切りで除外する T 番号 (例: "T3,T5")。

    Returns:
        dict:
          - overall_score: 0-10 の平均 (skip/None 除外)
          - score_max: 10
          - overall_severity: CRITICAL/HIGH/MEDIUM/LOW
          - summary_comment: 1-2 文の統合判定
          - detailed_scores: 各 tool の score ({"T1": 8.0, ...})
          - priority_issues: severity 順にソートされた issue list
          - tools_run / tools_skipped: 実行/除外 tool id
          - hint / source
    """
    if tex_or_text is None:
        tex_or_text = ""
    if bib is None:
        bib = ""
    if abstract is None:
        abstract = ""
    if author_last_names is None:
        author_last_names = ""

    skip_set = {s.strip().upper() for s in skip.split(",") if s.strip()}

    # 自動 skip: abstract 空 → T1, bib 空 → T5
    if not abstract.strip():
        skip_set.add("T1")
    if not bib.strip():
        skip_set.add("T5")

    tools_to_run: list[tuple[str, str, callable, tuple, dict]] = [
        ("T1", "abstract_strength",
         paper_writing_abstract_strength, (abstract,), {}),
        ("T2", "contribution_clarity_score",
         paper_writing_contribution_clarity_score, (tex_or_text,), {}),
        ("T3", "claim_quantification",
         paper_writing_claim_quantification, (tex_or_text,), {}),
        ("T4", "limitation_statement_presence",
         paper_writing_limitation_statement_presence, (tex_or_text,), {}),
        ("T5", "related_work_density",
         paper_writing_related_work_density, (tex_or_text,),
         {"bib": bib, "author_last_names": author_last_names}),
        ("T6", "figure_referencing_coverage",
         paper_writing_figure_referencing_coverage, (tex_or_text,), {}),
        ("T7", "given_new_ordering",
         paper_writing_given_new_ordering, (tex_or_text,), {}),
        # Tier 3 (v0.19.0)
        ("T10", "title_abstract_conclusion_triangle",
         paper_writing_title_abstract_conclusion_triangle,
         (tex_or_text,), {}),
        ("T11", "reproducibility_open_science_check",
         paper_writing_reproducibility_open_science_check,
         (tex_or_text,), {}),
        ("T14", "discussion_structure_4_elements",
         paper_writing_discussion_structure_4_elements,
         (tex_or_text,), {}),
    ]

    detailed_scores: dict[str, float] = {}
    detailed_results: dict[str, dict] = {}
    priority_issues: list[dict] = []

    for tid, name, func, args, kwargs in tools_to_run:
        if tid in skip_set:
            continue
        try:
            r = func(*args, **kwargs)
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

        # LOW = 肯定のみ、priority_issues から除外
        if severity != "LOW":
            priority_issues.append({
                "tool": tid,
                "name": name,
                "severity": severity,
                "score": score,
                "comments": r.get("comments", [])[:3],  # top-3 advisor comments
            })

    # severity でソート (CRITICAL > HIGH > MEDIUM > UNKNOWN > LOW)
    sev_rank = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 2, "UNKNOWN": 3, "LOW": 4}
    priority_issues.sort(
        key=lambda x: (sev_rank.get(x["severity"], 99),
                       -1 * (x.get("score") or 0))
    )

    # overall_score = 単純平均 (None 除外)
    valid_scores = [s for s in detailed_scores.values() if s is not None]
    overall_score = (round(sum(valid_scores) / len(valid_scores), 1)
                     if valid_scores else 0.0)
    overall_severity = _severity_from_score(overall_score)

    # summary_comment
    n_critical = sum(1 for i in priority_issues if i["severity"] == "CRITICAL")
    n_high = sum(1 for i in priority_issues if i["severity"] == "HIGH")
    if n_critical >= 2:
        summary = (
            f"CRITICAL 課題が {n_critical} 件。"
            "Abstract / Contribution / Claim 裏付け / Limitation / "
            "Related Work / Figure / Given-New のいずれかが根本的に欠けている。"
            "priority_issues の CRITICAL から順に対応。"
        )
    elif n_critical == 1:
        summary = (
            f"CRITICAL 課題 {priority_issues[0]['name']} が 1 件。"
            "まずここを直すと論文の印象が大きく変わる。"
        )
    elif n_high >= 2:
        summary = (
            f"HIGH 課題 {n_high} 件。致命傷は無いが、"
            "修正で査読通過確度が上がる余地がある。"
        )
    elif n_high == 1:
        summary = (
            f"HIGH 課題 {priority_issues[0]['name']} が 1 件。"
            "それ以外は概ね整合。"
        )
    else:
        summary = (
            "構造は整っている。最終パスは音読で "
            "著者の声で Given-New / Abstract / Contribution を最終判定。"
        )

    return {
        "overall_score": overall_score,
        "score_max": 10,
        "overall_severity": overall_severity,
        "summary_comment": summary,
        "detailed_scores": detailed_scores,
        "priority_issues": priority_issues,
        "tools_run": list(detailed_scores.keys()),
        "tools_skipped": sorted(skip_set),
        "hint": (
            "overall_score は参考値。priority_issues の CRITICAL/HIGH を "
            "1 つずつ skill.md §章と照合。"
            "abstract/bib が別ファイルなら別々に渡す。"
        ),
        "source": (
            "Plan B 統合レポート v0.19.0 — "
            "T1 Abstract / T2 Contribution / T3 Claim / T4 Limitation / "
            "T5 Related Work / T6 Figure Coverage / T7 Given-New / "
            "T10 Triangle / T11 Reproducibility / T14 Discussion 構造 "
            "を束ねる。"
        ),
    }
