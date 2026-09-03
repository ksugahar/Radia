"""Plan B T10: Title-Abstract-Conclusion 三位一体整合性
(paper_writing_title_abstract_conclusion_triangle).

論文の Title / Abstract / Conclusion は **同じ物語** を語るべきだが、
書き手は別々のセクションとして扱いがちで:

  - Title が過大広告 (abstract に書かれていない範囲を主張)
  - Conclusion が abstract と異なる主張
  - Abstract が title の core keyword を欠く

を起こす。reviewer は最初に Title → Abstract → (本文飛ばし) → Conclusion
を読むので、3 つの整合性は「信頼の一発勝負」。

設計方針:
- 3 セクションの主要 noun phrase / 数値 claim / methodology terms を抽出
- 三角形の各辺 (Title↔Abstract, Abstract↔Conclusion, Title↔Conclusion)
  で keyword overlap を計算
- score は 3 辺の harmonic mean (一辺でも弱ければ全体が下がる)
- 数値の食い違い (title:10x improvement vs conclusion:5x) を別途 flag
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# Section extraction (LaTeX or markdown)
# ------------------------------------------------------------------
_TITLE_PATTERNS = [
    re.compile(r"\\title\{([^{}]+)\}"),
    re.compile(r"^title\s*[:：]\s*(.+?)$", re.MULTILINE | re.IGNORECASE),
    re.compile(r"^#\s+(.+?)\s*$", re.MULTILINE),
]

_ABSTRACT_PATTERNS = [
    re.compile(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", re.DOTALL),
    re.compile(r"\\abstract\{([^{}]+)\}", re.DOTALL),
    re.compile(
        r"^#{1,3}\s*Abstract\s*\n+(.+?)(?=\n#{1,3}\s|\Z)",
        re.MULTILINE | re.DOTALL | re.IGNORECASE,
    ),
    re.compile(
        r"^#{1,3}\s*要[\s ]*旨\s*\n+(.+?)(?=\n#{1,3}\s|\Z)",
        re.MULTILINE | re.DOTALL,
    ),
]

_CONCLUSION_PATTERNS = [
    re.compile(
        r"\\section\*?\{(?:Conclusions?|結論|まとめ|結語)\}(.*?)"
        r"(?=\\section\b|\\bibliography|\Z)",
        re.DOTALL | re.IGNORECASE,
    ),
    re.compile(
        r"^#{1,3}\s*(?:Conclusions?|結論|まとめ|結語)\s*\n+(.+?)"
        r"(?=\n#{1,3}\s|\Z)",
        re.MULTILINE | re.DOTALL | re.IGNORECASE,
    ),
]


def _extract_first(patterns: list[re.Pattern], text: str) -> str:
    for pat in patterns:
        m = pat.search(text)
        if m:
            return m.group(1).strip()
    return ""


def _extract_keywords(text: str, top_k: int = 15) -> list[str]:
    """Extract normalized content words across Japanese and English."""
    keywords: list[str] = []
    seen: set[str] = set()
    # カタカナ 3+
    for m in re.finditer(r"[ァ-ヴー]{3,}", text):
        k = m.group(0)
        if k not in seen:
            seen.add(k); keywords.append(k)
    # 漢字 2+
    for m in re.finditer(r"[一-鿿]{2,}", text):
        k = m.group(0)
        if k not in seen and len(k) >= 2:
            seen.add(k); keywords.append(k)
    stop_words = {
        "the", "this", "that", "these", "those", "and", "or", "but",
        "for", "from", "with", "without", "into", "over", "under", "of",
        "in", "on", "at", "to", "by", "via", "using", "based", "paper",
        "study", "result", "results", "however", "therefore", "proposed",
    }
    for m in re.finditer(r"\b[A-Za-z][A-Za-z0-9\-]{2,}\b", text):
        k = m.group(0).casefold()
        if k in stop_words:
            continue
        if k not in seen:
            seen.add(k); keywords.append(k)
    return keywords[:top_k]


def _extract_numeric_claims(text: str) -> list[dict]:
    """数値+単位 / 倍率 / % の claim を抽出。"""
    claims: list[dict] = []
    # X% improvement / X 倍 / X-fold / X% reduction
    pat = re.compile(
        r"(\d+(?:\.\d+)?)\s*"
        r"(%|％|倍|fold|x|×|times)"
        r"\s*([A-Za-z]+|改善|向上|短縮|削減|低下|高速化)?",
        re.IGNORECASE,
    )
    for m in pat.finditer(text):
        claims.append({
            "value": m.group(1),
            "unit": m.group(2),
            "object": m.group(3) or "",
            "raw": m.group(0),
        })
    # 単独数値 + 単位
    pat2 = re.compile(
        r"(\d+(?:\.\d+)?)\s*"
        r"(kHz|MHz|GHz|Hz|kW|W|V|A|%|％|nm|μm|µm|mm|cm|m|km)"
    )
    for m in pat2.finditer(text):
        claims.append({
            "value": m.group(1),
            "unit": m.group(2),
            "object": "",
            "raw": m.group(0),
        })
    return claims


def _overlap_score(set_a: set[str], set_b: set[str]) -> tuple[float, list[str]]:
    """Jaccard-like overlap (0-1) と共通 keywords。"""
    if not set_a or not set_b:
        return (0.0, [])
    common = set_a & set_b
    score = len(common) / min(len(set_a), len(set_b))
    return (round(score, 3), sorted(common))


def paper_writing_title_abstract_conclusion_triangle(
    text: str
) -> dict:
    """Title / Abstract / Conclusion の三角形整合性を診断。

    3 辺 (Title↔Abstract, Abstract↔Conclusion, Title↔Conclusion) で
    keyword overlap を計算し、harmonic mean を score 化。
    数値 claim の食い違いも別途 flag。

    Args:
        text: 論文全文 (LaTeX or markdown)

    Returns:
        dict: score / triangle / numeric_conflicts / sections / comments
    """
    if text is None:
        text = ""

    title = _extract_first(_TITLE_PATTERNS, text)
    abstract = _extract_first(_ABSTRACT_PATTERNS, text)
    conclusion = _extract_first(_CONCLUSION_PATTERNS, text)

    sections = {
        "title_chars": len(title),
        "abstract_chars": len(abstract),
        "conclusion_chars": len(conclusion),
        "title_head": title[:120],
        "abstract_head": abstract[:200],
        "conclusion_head": conclusion[:200],
    }

    # Missing section detection
    missing: list[str] = []
    if not title:
        missing.append("title")
    if not abstract:
        missing.append("abstract")
    if not conclusion:
        missing.append("conclusion")

    if len(missing) >= 2:
        return {
            "score": None,
            "missing_sections": missing,
            "sections": sections,
            "comments": [
                f"3 セクションのうち {len(missing)} 個 ({missing}) が "
                "検出できないため診断不能。\\title{} / \\begin{abstract} / "
                "\\section{Conclusions} のような明示的 markup を推奨。"
            ],
            "hint": "section 抽出失敗。markup 確認推奨。",
            "source": "T10 triangle (paper-writing 2026-04)",
        }

    # Keyword extraction
    title_kw = set(_extract_keywords(title, top_k=10))
    abstract_kw = set(_extract_keywords(abstract, top_k=20))
    conclusion_kw = set(_extract_keywords(conclusion, top_k=20))

    # Triangle overlaps
    ta_score, ta_common = _overlap_score(title_kw, abstract_kw)
    ac_score, ac_common = _overlap_score(abstract_kw, conclusion_kw)
    tc_score, tc_common = _overlap_score(title_kw, conclusion_kw)

    triangle = {
        "title_abstract": {
            "overlap": ta_score,
            "common_keywords": ta_common[:10],
            "title_only": sorted(title_kw - abstract_kw)[:5],
        },
        "abstract_conclusion": {
            "overlap": ac_score,
            "common_keywords": ac_common[:10],
        },
        "title_conclusion": {
            "overlap": tc_score,
            "common_keywords": tc_common[:10],
        },
    }

    # Harmonic mean (一辺でも弱ければ大幅減点)
    eps = 1e-6
    edges = [ta_score + eps, ac_score + eps, tc_score + eps]
    harmonic = 3 / sum(1 / e for e in edges) - eps
    score = round(max(0.0, min(1.0, harmonic)) * 10, 1)

    # Numeric conflicts
    title_nums = _extract_numeric_claims(title)
    abstract_nums = _extract_numeric_claims(abstract)
    conclusion_nums = _extract_numeric_claims(conclusion)

    numeric_conflicts: list[dict] = []
    # Title の数値が abstract / conclusion に出てこない場合 = 過大広告
    for tn in title_nums:
        ab_has = any(tn["raw"] == a["raw"] for a in abstract_nums)
        co_has = any(tn["raw"] == c["raw"] for c in conclusion_nums)
        if not ab_has and not co_has:
            numeric_conflicts.append({
                "type": "title_number_unsupported",
                "title_claim": tn["raw"],
                "issue": (
                    f"Title の数値 '{tn['raw']}' が abstract / conclusion に "
                    "出てこない。Title 過大広告の疑い。"
                ),
            })

    # Abstract と Conclusion で数値が食い違う場合
    abs_nums_only = {a["raw"] for a in abstract_nums}
    con_nums_only = {c["raw"] for c in conclusion_nums}
    only_in_abs = abs_nums_only - con_nums_only
    only_in_con = con_nums_only - abs_nums_only
    if only_in_abs and only_in_con:
        # かつ似た object (改善 / 削減 等) の数値が両方にあるなら矛盾候補
        # ここでは単純に「両方にユニーク数値がある」場合だけ警告
        numeric_conflicts.append({
            "type": "numeric_value_divergence",
            "only_in_abstract": sorted(only_in_abs)[:5],
            "only_in_conclusion": sorted(only_in_con)[:5],
            "issue": (
                "Abstract と Conclusion で数値 claim が食い違っている可能性。"
                "同じ指標で異なる値が報告されていないか確認。"
            ),
        })

    # ---- comments ----
    comments: list[str] = []

    if missing:
        comments.append(f"検出できなかったセクション: {missing}。")

    if ta_score < 0.3:
        comments.append(
            f"Title-Abstract overlap が {ta_score:.2f} と低い。"
            "Title の core keyword が abstract に書かれていない可能性。"
            f"Title only keywords: {sorted(title_kw - abstract_kw)[:5]}"
        )
    if ac_score < 0.3:
        comments.append(
            f"Abstract-Conclusion overlap が {ac_score:.2f} と低い。"
            "Conclusion が abstract と別の物語を語っている可能性。"
            "reviewer は abstract で期待した内容と conclusion を照合する。"
        )
    if tc_score < 0.3:
        comments.append(
            f"Title-Conclusion overlap が {tc_score:.2f} と低い。"
            "論文の入口 (title) と出口 (conclusion) が繋がっていない。"
        )

    if numeric_conflicts:
        comments.append(
            f"数値 claim 不整合 {len(numeric_conflicts)} 件。"
            "Title 過大広告 / abstract-conclusion 値段差は reviewer の "
            "信頼を一発で失う。"
        )
        for nc in numeric_conflicts[:3]:
            comments.append(f"⚠ {nc['issue']}")

    if score >= 8 and not numeric_conflicts:
        comments.append(
            f"三角形整合性 {score}/10 — Title/Abstract/Conclusion が "
            "同じ物語を語れている。"
        )

    if not comments or comments == [f"検出できなかったセクション: {missing}。"]:
        comments.append(
            "三角形 keyword overlap は OK。最終確認は 3 セクションを "
            "順に音読し、reviewer 視点で物語整合性を判定。"
        )

    comments = comments[:8]

    hint = (
        "Harmonic mean を score にしているので、3 辺のうち最も弱い辺が "
        "全体を引っ張る。Title 過大広告は reviewer #2 トリガー筆頭の 1 つ。"
        "本ツールは keyword overlap ベースの hint — 意味論的整合性は "
        "別途音読で確認。"
    )

    return {
        "score": score,
        "score_max": 10,
        "missing_sections": missing,
        "sections": sections,
        "triangle": triangle,
        "numeric_conflicts": numeric_conflicts,
        "comments": comments,
        "hint": hint,
        "source": (
            "Day & Gastel 'How to Write and Publish a Scientific Paper' "
            "(7th ed.) Chapter 5 + Wallwork 'English for Writing Research "
            "Papers' Chapter 12 — Title/Abstract/Conclusion triangle "
            "consistency principle. (paper-writing 2026-04 採用)。"
        ),
    }
