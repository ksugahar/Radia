"""Plan B T15: Journal fit assessment
(paper_writing_journal_fit_assessment).

Submit 前に **target journal の aims & scope** と **論文の abstract /
keywords / title** を比較し、desk reject 確率を下げるための fit 診断。

editor が 30 秒で desk reject を判断する典型理由:
  - "Out of scope" — 論文 topic が journal の aims に合わない
  - "Not novel enough for this journal" — IF/prestige mismatch
  - "Wrong methodology emphasis" — theory journal に実験論文 等

本ツールは **caller が journal aims & scope を渡す** 形 (LLM/web 検索結果
を貼り付け) で fit を診断。

設計方針:
- abstract / title / keywords を抽出
- aims_text と main_text の token overlap
- per_aim_keyword の hit 状況
- score は overlap 率 (0-10) + missing_critical_terms の検出
- LLM-based ではなく keyword overlap で粗く出す。最終判断は人間。
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# Section extraction
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
]

_KEYWORDS_PATTERNS = [
    re.compile(r"\\keywords\{([^{}]+)\}"),
    re.compile(
        r"(?:keywords|key\s+words|キーワード)\s*[:：]?\s*"
        r"(.+?)(?:\n\n|\n#|\Z)",
        re.IGNORECASE | re.DOTALL,
    ),
]


def _extract_first(patterns: list[re.Pattern], text: str) -> str:
    for pat in patterns:
        m = pat.search(text)
        if m:
            return m.group(1).strip()
    return ""


def _extract_meaningful_terms(
    text: str, top_k: int = 30, min_kanji: int = 2, min_kata: int = 3
) -> list[str]:
    """Extract noun-like terms (kanji compounds, katakana, English uppercase)."""
    terms: list[str] = []
    seen: set[str] = set()
    # 漢字 2+
    for m in re.finditer(rf"[一-鿿]{{{min_kanji},}}", text):
        k = m.group(0)
        if k not in seen:
            seen.add(k); terms.append(k)
    # カタカナ 3+
    for m in re.finditer(rf"[ァ-ヴー]{{{min_kata},}}", text):
        k = m.group(0)
        if k not in seen:
            seen.add(k); terms.append(k)
    # 英大文字始まり (Pascal/UPPER)
    for m in re.finditer(r"\b[A-Z][A-Za-z0-9\-]{2,}\b", text):
        k = m.group(0)
        if k not in seen:
            seen.add(k); terms.append(k)
    # 英小文字 noun-like (5+ chars, 一般語除外)
    common_en = {
        "paper", "study", "research", "result", "method", "model",
        "approach", "system", "analysis", "based", "novel", "this",
        "that", "with", "from", "have", "been", "their", "where",
        "which", "while", "between",
    }
    for m in re.finditer(r"\b[a-z][a-z\-]{4,}\b", text):
        k = m.group(0)
        if k not in seen and k.lower() not in common_en:
            seen.add(k); terms.append(k)
    return terms[:top_k]


def paper_writing_journal_fit_assessment(
    text: str,
    journal_name: str,
    aims_and_scope: str,
    additional_notes: str = "",
) -> dict:
    """target journal の aims & scope と論文の fit を診断。

    Args:
        text: 論文 (LaTeX or markdown、abstract / title が含まれていること)
        journal_name: target journal 名 (例: "IEEE TPS")
        aims_and_scope: journal の aims & scope 本文 (web から貼り付け)
        additional_notes: 追加の制約 (例: "Special issue on quantum systems")

    Returns:
        dict: score / paper_terms / aims_terms / overlap / missing_critical /
              comments
    """
    if text is None:
        text = ""
    if aims_and_scope is None:
        aims_and_scope = ""

    if not aims_and_scope.strip():
        return {
            "score": None,
            "comments": [
                "aims_and_scope が空。target journal の Web ページから "
                "「Aims and Scope」セクションを貼り付けて再投入してください。"
            ],
            "hint": "input が不十分。",
            "source": "T15 journal fit (paper-writing 2026-04)",
        }

    title = _extract_first(_TITLE_PATTERNS, text)
    abstract = _extract_first(_ABSTRACT_PATTERNS, text)
    keywords_block = _extract_first(_KEYWORDS_PATTERNS, text)

    if not title and not abstract:
        return {
            "score": None,
            "comments": [
                "論文 text から title / abstract が抽出できなかった。"
                "\\title{} / \\begin{abstract} markup を確認、または "
                "abstract 本文を直接 text 引数に含めてください。"
            ],
            "hint": "section 抽出失敗。",
            "source": "T15 journal fit",
        }

    # Combine paper-side identification text
    paper_signal_text = " ".join([
        title, abstract, keywords_block, additional_notes
    ])

    paper_terms = _extract_meaningful_terms(paper_signal_text, top_k=40)
    aims_terms = _extract_meaningful_terms(aims_and_scope, top_k=40)

    paper_set = set(paper_terms)
    aims_set = set(aims_terms)

    common = paper_set & aims_set
    paper_only = paper_set - aims_set
    aims_only = aims_set - paper_set

    # Coverage from journal perspective: how much of journal's vocabulary
    # is reflected in paper
    journal_coverage = len(common) / len(aims_set) if aims_set else 0
    # Coverage from paper perspective: how much of paper's vocabulary
    # is in journal scope
    paper_coverage = len(common) / len(paper_set) if paper_set else 0
    # Symmetric F1-like
    if journal_coverage + paper_coverage > 0:
        f1 = 2 * journal_coverage * paper_coverage / (
            journal_coverage + paper_coverage
        )
    else:
        f1 = 0.0

    score = round(f1 * 10, 1)

    # ---- comments ----
    comments: list[str] = []

    comments.append(
        f"Journal '{journal_name}' fit: F1={score}/10 "
        f"(journal_coverage={journal_coverage:.2f}, "
        f"paper_coverage={paper_coverage:.2f})。"
    )

    if score < 4:
        comments.append(
            f"⚠ score {score} は LOW — desk reject リスク高。"
            "Out-of-scope 判定の可能性。"
            "別 journal 検討、または abstract / introduction で "
            "本研究と journal aims の橋渡し文を 1-2 文追加。"
        )
    elif score < 6:
        comments.append(
            f"score {score} は MODERATE。submit 可能だが、cover letter で "
            "明示的に「本研究は journal の X 領域に該当」を述べると安全。"
        )
    else:
        comments.append(
            f"score {score} は GOOD — fit は良好。"
        )

    # Show top common keywords
    common_sorted = sorted(common)
    if common_sorted:
        comments.append(
            f"共通 keywords (top 10): {common_sorted[:10]}"
        )

    # Identify potentially critical missing terms (keywords in aims that
    # are noun-like and frequently appear → not in paper)
    if aims_only:
        # Sort aims_only by their frequency in aims_and_scope
        freq_in_aims = {
            t: aims_and_scope.lower().count(t.lower())
            for t in aims_only
        }
        top_missing = sorted(
            aims_only, key=lambda t: -freq_in_aims[t]
        )[:8]
        comments.append(
            f"Journal aims にあるが論文に不在の term (頻度順 top 8): "
            f"{top_missing}。これらが目立つ term なら abstract 微修正で "
            "言及するだけで fit score が上がる。"
        )

    if paper_only and len(paper_only) >= 10:
        # Show terms paper uses that journal doesn't — could indicate scope drift
        paper_only_sorted = sorted(paper_only)[:8]
        comments.append(
            f"論文側のみで journal aims に無い term (例): {paper_only_sorted}。"
            "本研究の独自 contribution であれば良いが、journal scope と "
            "齟齬があれば cover letter で justify。"
        )

    if not title:
        comments.append(
            "ℹ Title が抽出できなかった。\\title{...} markup 確認。"
        )

    comments = comments[:7]

    hint = (
        "F1 score = harmonic mean of paper_coverage と journal_coverage。"
        "片方だけ高くてもダメ (e.g. journal terms 全部使ってるが論文の "
        "contribution が無関係 = paper_coverage 低)。"
        "本ツールは keyword overlap ベース — 意味論的 fit は別途 editor 視点で確認。"
        "投稿前に journal の published article と本論文の topic 距離も別途 check。"
    )

    return {
        "score": score,
        "score_max": 10,
        "journal_name": journal_name,
        "f1": round(f1, 3),
        "journal_coverage": round(journal_coverage, 3),
        "paper_coverage": round(paper_coverage, 3),
        "n_paper_terms": len(paper_set),
        "n_aims_terms": len(aims_set),
        "n_common": len(common),
        "common_keywords": sorted(common)[:20],
        "paper_only_keywords": sorted(paper_only)[:20],
        "aims_only_keywords": sorted(aims_only)[:20],
        "comments": comments,
        "hint": hint,
        "source": (
            "Editor desk reject 統計 (e.g. Bohannon 2013 Science 342:60-65 "
            "に類するもの) + journal manuscript central 編集者 guidance。"
            "本ツールは aims-paper keyword overlap で粗く fit 推定する "
            "pre-submission check。"
            "(paper-writing 2026-04 採用)。"
        ),
    }
