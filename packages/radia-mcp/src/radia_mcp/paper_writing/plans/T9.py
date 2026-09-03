"""Plan B T9: reviewer-2 trigger summary (paper_writing_reviewer_2_trigger_summary).

paper_writing 配下の 5-7 属性ツールを weighted union で統合し、
「悪意ある reviewer-2 が最も突いてきそうな point」を attack-probability
順に列挙する meta 診断。

設計方針:
- score は 0-10 の粗い温度計。**逆スケール**: 高 score = safe (reviewer
  攻撃 surface 小)、低 score = reviewer が攻撃する可能性が高い。
- PRIMARY OUTPUT は comments (reviewer-2 の声で書かれた attack
  vector 警告)。
- observations は raw measurement (各 trigger の weight / count /
  contribution)。自動化 numeric 判定には使わない。
- 他の T3-T7 + cross_lint / tools.py の既存ツールを**呼び出す
  だけ**で、新しい lint ロジックは実装しない (薄い orchestrator)。

呼び出す下位ツール:
  - paper_writing_claim_quantification        (T3) -> unquantified_count
  - paper_writing_related_work_density         (T5) -> self_cite_ratio
  - paper_writing_limitation_statement_presence (T4) -> limitation_paragraphs
  - paper_writing_figure_referencing_coverage  (T6) -> unreferenced_labels
  - paper_writing_find_undefined_acronyms      (cross_lint) -> undefined_count
  - paper_writing_count_weak_expressions       (tools.py) -> total_weak_expressions
  - paper_writing_check_english_redflags       (tools.py) -> total_issues

skill.md NG #1-#7 / Wallwork §8.12 (quantification) / §9.12 (limitation) /
§14 (related work) / reviewer-2 archetype (MCP コミュニティ事例) に基づく。
"""

from __future__ import annotations


# ------------------------------------------------------------------
# Trigger weights
# ------------------------------------------------------------------
# (trigger_name, weight, human-readable label)
_TRIGGERS_SPEC = [
    ("unquantified_hype",  3.0, "Unquantified hype claims"),
    ("self_cite_excess",   2.5, "Self-cite 過剰 (> 20%)"),
    ("no_limitation",      2.5, "Limitation 不在"),
    ("unreferenced_figures", 2.0, "Unreferenced figures"),
    ("undefined_acronyms", 1.5, "Undefined acronyms"),
    ("weak_expressions",   1.5, "Weak hedges"),
    ("english_redflags",   1.0, "English redflags"),
]

# A "typical max" per trigger when normalising raw counts to 0-1.
# Counts >= _TYPICAL_MAX が「攻撃 probability が概ね頭打ち」
# になる目安値。粗いノブなので調整余地あり。
_TYPICAL_MAX = {
    "unquantified_hype":    5.0,   # 5 件以上の hype は "bad enough"
    "self_cite_excess":     1.0,   # bool
    "no_limitation":        1.0,   # bool
    "unreferenced_figures": 3.0,   # 3 件以上で頭打ち
    "undefined_acronyms":   5.0,
    "weak_expressions":     5.0,
    "english_redflags":     5.0,
}

# Self-cite が「過剰」とみなす下限比率 (skill.md NG #5 / Wallwork §14)
_SELF_CITE_THRESHOLD = 0.20


def _safe_get(d: dict, *keys, default=None):
    """Nested dict から安全に値を取り出す。途中で None / key 欠損なら default。"""
    cur = d
    for k in keys:
        if not isinstance(cur, dict):
            return default
        if k not in cur:
            return default
        cur = cur[k]
    return cur


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def paper_writing_reviewer_2_trigger_summary(
    tex_or_text: str,
    bib: str = "",
    author_last_names: str = "",
) -> dict:
    """悪意ある reviewer-2 が最も突いてくるポイントを weighted union で列挙。

    paper_writing 内の 5-7 属性ツール (T3/T4/T5/T6/cross_lint/tools.py) を
    統合し、reviewer attack probability 順に trigger をソート。
    skill.md NG #1-#7 + reviewer-2 archetype に基づく。

    引数:
        tex_or_text: 論文原稿 (LaTeX 全文 / プレーンテキスト)。
        bib: .bib ファイル内容。省略時は self-cite check を skip。
        author_last_names: 自分+共著者の姓をカンマ区切り (self-cite 判定用)。

    返り値:
        dict: score / score_max / comments / observations / hint / source
            score は逆スケール (高 = safe、低 = reviewer が攻撃する可能性高)。
    """
    if tex_or_text is None:
        tex_or_text = ""
    if bib is None:
        bib = ""
    if author_last_names is None:
        author_last_names = ""

    # 遅延 import (循環依存防止、テスト時の mock 容易化)
    from .T3 import paper_writing_claim_quantification
    from .T4 import paper_writing_limitation_statement_presence
    from .T5 import paper_writing_related_work_density
    from .T6 import paper_writing_figure_referencing_coverage
    from ..cross_lint import paper_writing_find_undefined_acronyms
    from ..tools import (
        paper_writing_count_weak_expressions,
        paper_writing_check_english_redflags,
    )

    # --------------------------------------------------------------
    # 1. 下位ツールを呼ぶ (例外は "その trigger は 0" として扱う)
    # --------------------------------------------------------------
    counts: dict[str, float] = {}
    raw_extra: dict[str, dict] = {}
    unknown_triggers: set[str] = set()

    # unquantified_hype (T3)
    try:
        r_t3 = paper_writing_claim_quantification(tex_or_text)
        counts["unquantified_hype"] = float(
            _safe_get(r_t3, "observations", "unquantified_count", default=0)
            or 0
        )
        raw_extra["unquantified_hype"] = {
            "source_tool": "paper_writing_claim_quantification",
            "example": (
                _safe_get(r_t3, "observations", "unquantified_examples", default=[])
                or []
            )[:1],
        }
    except Exception as e:  # noqa: BLE001
        counts["unquantified_hype"] = 0.0
        unknown_triggers.add("unquantified_hype")
        raw_extra["unquantified_hype"] = {"error": str(e)}

    # self_cite_excess (T5). bib 空なら skip (count=0)。
    if bib.strip():
        try:
            r_t5 = paper_writing_related_work_density(
                tex_or_text,
                bib=bib,
                author_last_names=author_last_names,
            )
            ratio = float(
                _safe_get(r_t5, "observations", "self_cite_ratio", default=0.0)
                or 0.0
            )
            # bool trigger: excess = (ratio > 0.20)
            counts["self_cite_excess"] = 1.0 if ratio > _SELF_CITE_THRESHOLD else 0.0
            raw_extra["self_cite_excess"] = {
                "source_tool": "paper_writing_related_work_density",
                "self_cite_ratio": ratio,
                "threshold": _SELF_CITE_THRESHOLD,
                "skipped": False,
            }
        except Exception as e:  # noqa: BLE001
            counts["self_cite_excess"] = 0.0
            unknown_triggers.add("self_cite_excess")
            raw_extra["self_cite_excess"] = {"error": str(e), "skipped": False}
    else:
        counts["self_cite_excess"] = 0.0
        unknown_triggers.add("self_cite_excess")
        raw_extra["self_cite_excess"] = {
            "source_tool": "paper_writing_related_work_density",
            "skipped": True,
            "skipped_reason": "bib 未指定のため self-cite 判定不可",
        }

    # no_limitation (T4): limitation_paragraphs が空なら trigger on
    try:
        r_t4 = paper_writing_limitation_statement_presence(tex_or_text)
        lim_paras = (
            _safe_get(r_t4, "observations", "limitation_paragraphs", default=[])
            or []
        )
        disc_found = bool(
            _safe_get(r_t4, "observations", "discussion_found", default=False)
        )
        # Discussion 未検出なら trigger を ON にはしない (検査対象外)
        if disc_found and not lim_paras:
            counts["no_limitation"] = 1.0
        else:
            counts["no_limitation"] = 0.0
        raw_extra["no_limitation"] = {
            "source_tool": "paper_writing_limitation_statement_presence",
            "discussion_found": disc_found,
            "limitation_paragraph_count": len(lim_paras),
        }
    except Exception as e:  # noqa: BLE001
        counts["no_limitation"] = 0.0
        unknown_triggers.add("no_limitation")
        raw_extra["no_limitation"] = {"error": str(e)}

    # unreferenced_figures (T6)
    try:
        r_t6 = paper_writing_figure_referencing_coverage(tex_or_text)
        unref = (
            _safe_get(r_t6, "observations", "unreferenced_labels", default=[])
            or []
        )
        counts["unreferenced_figures"] = float(len(unref))
        raw_extra["unreferenced_figures"] = {
            "source_tool": "paper_writing_figure_referencing_coverage",
            "unreferenced_labels": unref[:5],
        }
    except Exception as e:  # noqa: BLE001
        counts["unreferenced_figures"] = 0.0
        unknown_triggers.add("unreferenced_figures")
        raw_extra["unreferenced_figures"] = {"error": str(e)}

    # undefined_acronyms (cross_lint)
    try:
        r_acro = paper_writing_find_undefined_acronyms(tex_or_text)
        counts["undefined_acronyms"] = float(
            r_acro.get("undefined_count", 0) or 0
        )
        raw_extra["undefined_acronyms"] = {
            "source_tool": "paper_writing_find_undefined_acronyms",
            "sample": [
                u.get("acronym")
                for u in (r_acro.get("undefined", []) or [])[:3]
            ],
        }
    except Exception as e:  # noqa: BLE001
        counts["undefined_acronyms"] = 0.0
        unknown_triggers.add("undefined_acronyms")
        raw_extra["undefined_acronyms"] = {"error": str(e)}

    # weak_expressions (tools.py)
    try:
        r_weak = paper_writing_count_weak_expressions(tex_or_text)
        counts["weak_expressions"] = float(
            r_weak.get("total_weak_expressions", 0) or 0
        )
        raw_extra["weak_expressions"] = {
            "source_tool": "paper_writing_count_weak_expressions",
        }
    except Exception as e:  # noqa: BLE001
        counts["weak_expressions"] = 0.0
        unknown_triggers.add("weak_expressions")
        raw_extra["weak_expressions"] = {"error": str(e)}

    # english_redflags (tools.py)
    try:
        r_rf = paper_writing_check_english_redflags(tex_or_text)
        # 仕様書では "total_red_flags"、実装では "total_issues"。
        # 互換のため両方見る。
        rf_n = (
            r_rf.get("total_red_flags")
            or r_rf.get("total_issues")
            or 0
        )
        counts["english_redflags"] = float(rf_n)
        raw_extra["english_redflags"] = {
            "source_tool": "paper_writing_check_english_redflags",
        }
    except Exception as e:  # noqa: BLE001
        counts["english_redflags"] = 0.0
        unknown_triggers.add("english_redflags")
        raw_extra["english_redflags"] = {"error": str(e)}

    # --------------------------------------------------------------
    # 2. Weighted union (contribution = weight * normalised count)
    # --------------------------------------------------------------
    triggers: list[dict] = []
    risk_score_raw = 0.0
    total_possible = 0.0
    for name, weight, label in _TRIGGERS_SPEC:
        c = counts.get(name, 0.0)
        if name in unknown_triggers:
            triggers.append({
                "name": name,
                "label": label,
                "weight": weight,
                "count": None,
                "normalised": None,
                "contribution": None,
                "status": "unknown",
                "extra": raw_extra.get(name, {}),
            })
            continue
        typ = _TYPICAL_MAX.get(name, 5.0)
        if typ <= 0:
            typ = 1.0
        norm = min(1.0, c / typ)
        contrib = weight * norm
        risk_score_raw += contrib
        total_possible += weight  # weight * 1.0
        triggers.append({
            "name": name,
            "label": label,
            "weight": weight,
            "count": c,
            "normalised": round(norm, 3),
            "contribution": round(contrib, 3),
            "extra": raw_extra.get(name, {}),
        })

    # risk_percent: 0% = 完全 safe、100% = reviewer 確実に攻撃
    if total_possible > 0:
        risk_percent = min(100.0, (risk_score_raw / total_possible) * 100.0)
    else:
        risk_percent = 0.0

    # score は **逆スケール** (高 = safe)
    # risk_percent 0 → score 10, risk_percent 100 → score 0
    score = round(max(0.0, 10.0 - risk_percent / 10.0), 1)

    # contribution 降順 (attack 確率高い順)
    triggers_sorted = sorted(
        triggers,
        key=lambda t: (
            t["contribution"] if t["contribution"] is not None else -1.0
        ),
        reverse=True,
    )
    top_triggers = [
        t["name"] for t in triggers_sorted if (t["contribution"] or 0) > 0
    ][:3]

    observations = {
        "risk_score_raw": round(risk_score_raw, 3),
        "risk_percent": round(risk_percent, 1),
        "triggers": triggers_sorted,
        "top_triggers": top_triggers,
        "unknown_triggers": sorted(unknown_triggers),
    }

    # --------------------------------------------------------------
    # 3. Comments (PRIMARY — reviewer-2 の声で)
    # --------------------------------------------------------------
    comments: list[str] = _build_comments(
        triggers_sorted,
        top_triggers,
        counts,
        raw_extra,
        score,
        risk_percent,
    )

    hint = (
        "5-7 属性の weighted union。実際の reviewer 攻撃は分野・journal 文化で "
        "変動 (pure EE より CS top conference で hype に辛い等)。"
        "risk_percent は絶対値ではなく『相対的に怪しい度合い』の指標。"
    )

    source = (
        "Plan B reviewer-2 archetype (v0.14.0) — skill.md NG #1-#7 / "
        "Wallwork §8.12 (quantification) / §9.12 (limitation) / §14 "
        "(related work) + MCP コミュニティ事例"
    )

    return {
        "score": float(score),
        "score_max": 10,
        "comments": comments,
        "observations": observations,
        "hint": hint,
        "source": source,
    }


# ------------------------------------------------------------------
# Comment builder
# ------------------------------------------------------------------
def _build_comments(
    triggers_sorted: list[dict],
    top_triggers: list[str],
    counts: dict[str, float],
    raw_extra: dict[str, dict],
    score: float,
    risk_percent: float,
) -> list[str]:
    """上位 trigger を reviewer-2 の声で comment 化。2-5 件。"""
    comments: list[str] = []

    if not top_triggers:
        # すべての trigger が 0 寄与 (safe)
        comments.append(
            "reviewer-2 attack surface は最小化済み。主要 trigger "
            "(hype/self-cite/limitation/figure/acronym) すべて clear。"
        )
        comments.append(
            "ただし本 tool は heuristic。分野特有の『novelty 弱い』"
            "『実験が少ない』系の trigger は検出できない — "
            "手動で skill.md NG #1-#7 を読み返すこと。"
        )
        return comments[:5]

    # 各 top trigger に対して reviewer-2 風の警告を生成
    for name in top_triggers:
        n = int(counts.get(name, 0.0))
        extra = raw_extra.get(name, {})
        line = _comment_for(name, n, extra)
        if line:
            comments.append(line)

    # 軽い summary (risk % を添える)
    if risk_percent >= 50:
        comments.append(
            f"総合 risk {int(round(risk_percent))}% — reviewer-2 は "
            f"上記 {len(top_triggers)} 箇所を優先 attack してくる。"
            "rebuttal を書く前に修正せよ。"
        )
    elif risk_percent >= 20:
        comments.append(
            f"総合 risk {int(round(risk_percent))}% — "
            "主要 trigger が散在。上から順に潰すと attack surface が激減する。"
        )
    else:
        comments.append(
            f"総合 risk {int(round(risk_percent))}% — "
            "攻撃面は小さい。残り trigger は cosmetic。"
        )

    return comments[:5]


def _comment_for(name: str, n: int, extra: dict) -> str:
    """trigger 名に応じた reviewer-2 voice comment を返す。"""
    if name == "unquantified_hype":
        return (
            f"Reviewer 2 の声: 「{n} 箇所の hype 表現に数値/引用根拠がない。"
            "'Our method is significantly better' — better than what? "
            "by how much?」 → paper_writing_claim_quantification 参照、"
            "各 hype 文に % か \\cite を添える。"
        )
    if name == "self_cite_excess":
        ratio = extra.get("self_cite_ratio", 0.0)
        pct = int(round(ratio * 100))
        return (
            f"Reviewer 2 の声: 「self-cite 率 {pct}% は suspicious。"
            "Are you hiding relevant prior work?」 → 競合研究 5-10 本を "
            "Intro に足す (skill.md NG #5 / Wallwork §14)。"
        )
    if name == "no_limitation":
        return (
            "Reviewer 2 の声: 「Where is the limitation discussion? "
            "Authors don't seem to understand their own work's "
            "weaknesses.」 → Discussion 中盤に 150 字以上で制約を明示 "
            "(skill.md §O / Wallwork §9.12)。"
        )
    if name == "unreferenced_figures":
        labels = extra.get("unreferenced_labels", [])
        head = labels[0] if labels else "?"
        return (
            f"Reviewer 2 の声: 「Figure/Table {head} is never referenced "
            f"in the text ({n} labels total). Why is it here?」 → "
            "本文に \\ref を追加 or 図を削る (Tufte / Wallwork §17)。"
        )
    if name == "undefined_acronyms":
        sample = extra.get("sample", [])
        acro = sample[0] if sample else "ACRO"
        return (
            f"Reviewer 2 の声: 「What does '{acro}' mean? "
            f"Never defined ({n} 件)。」 → 初出に (full form) を追加。"
        )
    if name == "weak_expressions":
        return (
            f"Reviewer 2 の声: 「'may', 'might', 'tends to' 系 hedge が "
            f"{n} 件。Authors seem unsure of their own claims.」 → "
            "Conclusion/Contribution 文の hedge を削る (Reviewer 2 trigger)。"
        )
    if name == "english_redflags":
        return (
            f"Reviewer 2 の声: 「English redflags ({n} 件) — "
            "'In the paper' / 'we will' / 'Please refer to' 等。"
            "Needs native proofreading before submission.」 → "
            "paper_writing_check_english_redflags で個別修正。"
        )
    return ""
