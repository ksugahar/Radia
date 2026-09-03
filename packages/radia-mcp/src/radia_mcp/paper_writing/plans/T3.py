"""Plan B T3: Unquantified Claim 検出 (paper_writing_claim_quantification).

skill.md NG #7 / Wallwork §8.12 / reviewer-2 archetype に基づく粗い温度計。
"significantly better", "substantially faster" のような hype adjective を
含むが、同一文内に数値 / 引用 / 図参照が無い文を flag する。

設計方針:
- score は 0-10 の粗い温度計。hype 文に裏付けが無いほど score が下がる。
- PRIMARY OUTPUT は comments (人間レビュアー風の助言文)。
- observations は raw measurement。自動化 numeric 判定には使わない。
- hype word list は粗い。文脈 (Conclusion vs Related Work) で許容度が
  変わるため、comments と手判断で最終確認する。
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# Patterns
# ------------------------------------------------------------------
# hype adjective / adverb / verb (skill.md NG #7 / reviewer-2 archetype)
_HYPE_WORDS = [
    "significantly",
    "substantially",
    "faster",
    "better",
    "efficient",
    "superior",
    "excellent",
    "outperform",
    "remarkable",
    "notable",
    "effective",
]

_HYPE_RE = re.compile(
    r"\b(?:" + "|".join(_HYPE_WORDS) + r")\w*\b",
    re.IGNORECASE,
)

# 数値: 12%, 3.2×, 4x, 10fold, 10 times, 10倍, 100 kHz, 3.2 MHz, 5 ms, 2 s,
# 10 mm, 20 μm など.  \d+[\.×x%] を含み、かつ単位付き数値も捕捉する.
_NUMBER_RE = re.compile(
    r"\b\d+\.\d+\b"
    r"|\b\d+(?:\.\d+)?\s*(?:%|倍|fold|times|×|kHz|MHz|GHz|Hz|ms|"
    r"mm|μm|um|nm|dB|s)(?![A-Za-z])",
    re.IGNORECASE,
)

# 引用: LaTeX \cite{...} または [1] / [12,13] 形式の reference bracket
_CITE_RE = re.compile(
    r"\\cite[a-zA-Z]*\s*\{[^}]*\}"
    r"|\[\s*\d+(?:\s*[,\-]\s*\d+)*\s*\]"
)

# 図表参照: \ref / \autoref / Fig. 3 / 図 3
_FIGREF_RE = re.compile(
    r"\\(?:ref|autoref|cref|Cref|eqref)\s*\{[^}]*\}"
    r"|\bFig(?:ure)?\.?\s*\d+"
    r"|\bTable\s*\d+"
    r"|図\s*\d+"
    r"|表\s*\d+",
    re.IGNORECASE,
)

# 文分割: . / 。 / ． (英日混在)
_SENTENCE_SPLIT_RE = re.compile(r"(?<=[.。．])\s+|(?<=[。．])")


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def _split_sentences(text: str) -> list[str]:
    """英日混在の単純な文分割。空文字は除去する。

    次のケースでは分割しない:
    - 小数点  (e.g. "3.2×")       → 前後が数字
    - 略語    (Fig., Eq., et al.) → placeholder に退避
    """
    if not text:
        return []
    # 改行→空白に正規化
    normalized = re.sub(r"\s+", " ", text).strip()
    if not normalized:
        return []

    # Protect abbreviations that end with '.' followed by a digit or word.
    protected = normalized
    abbrevs = [
        (r"\bFig\.", "FIGABBR"),
        (r"\bFigs\.", "FIGSABBR"),
        (r"\bEq\.", "EQABBR"),
        (r"\bEqs\.", "EQSABBR"),
        (r"\bet al\.", "ETALABBR"),
        (r"\be\.g\.", "EGABBR"),
        (r"\bi\.e\.", "IEABBR"),
        (r"\bvs\.", "VSABBR"),
        (r"\bRef\.", "REFABBR"),
        (r"\bNo\.", "NOABBR"),
    ]
    for pat, tok in abbrevs:
        protected = re.sub(pat, tok, protected, flags=re.IGNORECASE)

    # Protect decimals: digit.digit  →  digit<DOT>digit
    protected = re.sub(r"(\d)\.(\d)", r"\1<DOT>\2", protected)

    # Split on . / 。 / ．
    parts = re.split(r"(?<=[.。．])\s+", protected)

    # Restore placeholders
    restored: list[str] = []
    for p in parts:
        p = p.replace("<DOT>", ".")
        for pat, tok in abbrevs:
            # pat like r"\bFig\." → literal replacement
            literal = pat.replace(r"\b", "").replace(r"\.", ".")
            p = p.replace(tok, literal)
        p = p.strip()
        if p:
            restored.append(p)
    return restored


def _first_hype(sentence: str) -> str | None:
    """文中最初の hype word を返す。見つからなければ None。"""
    m = _HYPE_RE.search(sentence)
    return m.group(0) if m else None


def _truncate(s: str, n: int = 120) -> str:
    return s if len(s) <= n else s[: n - 1] + "…"


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def paper_writing_claim_quantification(text: str) -> dict:
    """Unquantified hype claims を検出。

    'significantly better' のような hype adjective を含むが、同一文内に
    数値 / 引用 / 図参照が無い文を flag する。
    skill.md NG #7 / Wallwork §8.12 / reviewer-2 archetype。

    引数:
        text: 論文原稿 (LaTeX またはプレーンテキスト)

    返り値:
        dict: score / score_max / comments / observations / hint / source
    """
    if text is None:
        text = ""

    sentences = _split_sentences(text)
    total_sentences = len(sentences)

    hype_sentences: list[dict] = []
    unquantified_examples: list[dict] = []
    unquantified_count = 0

    for idx, sent in enumerate(sentences):
        hype_word = _first_hype(sent)
        if not hype_word:
            continue
        has_number = bool(_NUMBER_RE.search(sent))
        has_cite = bool(_CITE_RE.search(sent))
        has_figref = bool(_FIGREF_RE.search(sent))
        entry = {
            "pos": idx,
            "sentence": sent,
            "hype_word": hype_word,
            "has_number": has_number,
            "has_cite": has_cite,
            "has_figref": has_figref,
        }
        hype_sentences.append(entry)
        if not (has_number or has_cite or has_figref):
            unquantified_count += 1
            if len(unquantified_examples) < 5:
                unquantified_examples.append(
                    {
                        "pos": idx,
                        "sentence": _truncate(sent, 120),
                        "hype_word": hype_word,
                    }
                )

    observations = {
        "total_sentences": total_sentences,
        "hype_sentences": hype_sentences,
        "unquantified_count": unquantified_count,
        "unquantified_examples": unquantified_examples,
    }

    # ---- score ----
    score = max(0.0, 10.0 - float(unquantified_count))

    # ---- comments ----
    comments: list[str] = []

    if total_sentences == 0:
        comments.append("文検出ゼロ。入力を確認。")
        comments.append(
            "skill.md NG #7 / Wallwork §8.12 は定量裏付けを要求する。"
            "本ツールは hype word の周囲に数値 / \\cite / 図参照があるかを見る。"
        )
    elif unquantified_count > 0:
        first_ex = unquantified_examples[0]["sentence"] if unquantified_examples else ""
        comments.append(
            f"定量化されていない hype 文が {unquantified_count} 件。"
            f"例: 「{first_ex}」。"
            "具体数値 (%, ×, kHz) or \\cite / Fig. 番号で裏付けないと "
            "reviewer-2 trigger (skill.md NG #7, Wallwork §8.12)。"
        )
        if unquantified_count >= 3:
            comments.append(
                f"hype 多用 ({unquantified_count} 件) は paper 全体の tone を"
                "「誇大宣伝」に傾ける。Conclusion/Abstract で特に削減し、"
                "Related Work ではむしろ中立語 ('differs', 'contrasts') を使う。"
            )
        if len(unquantified_examples) > 1:
            extras = "; ".join(
                f"「{e['sentence']}」" for e in unquantified_examples[1:3]
            )
            comments.append(f"他の例: {extras}")
        comments.append(
            "書き換え例: 'significantly faster' → '3.2× faster (Fig. 5)' / "
            "'better accuracy' → 'accuracy improved by 4.1% (Table 2)' / "
            "'outperforms prior work' → 'outperforms \\cite{prev} by 12%'"
        )
    else:
        # unquantified_count == 0 and total_sentences > 0
        if hype_sentences:
            comments.append(
                f"hype 表現 ({len(hype_sentences)} 件) はすべて定量要素付き "
                "(数値 / \\cite / Fig. 参照)。reviewer-2 防衛済み。"
            )
        else:
            comments.append(
                "hype 表現は検出されなかった。reviewer-2 の誇大宣伝攻撃は"
                "防衛済みだが、主張の強度自体が不足していないかは別途確認。"
            )
        comments.append(
            "ただし hype word list は粗い。'state-of-the-art', 'novel' 等は"
            "本チェックに含まれないため、skill.md NG #7 を再読して"
            "人間判断で補完すること。"
        )

    # 2-5 個担保
    if len(comments) < 2:
        comments.append(
            "score は粗い温度計。skill.md NG #7 / Wallwork §8.12 を"
            "直接参照して人間判断で補完すること。"
        )
    comments = comments[:5]

    hint = (
        "hype word list は粗い。文脈 (Conclusion vs Related Work) で"
        "許容度が変わる — comments と手判断で最終確認。"
    )
    source = "skill.md NG #7 / Wallwork §8.12 / reviewer-2 archetype"

    return {
        "score": float(score),
        "score_max": 10,
        "comments": comments,
        "observations": observations,
        "hint": hint,
        "source": source,
    }
