"""Plan B T7: Given-New 情報配置診断 (paper_writing_given_new_ordering).

Wallwork §3.4-3.6 / skill.md §D / Halliday theme-rheme に基づく粗い温度計。
各文は末尾に新情報 (rheme) を置き、次文冒頭でそれを既知情報 (theme)
として受け直す — この連結が英語論文トップ誌で最も効く linguistic
technique とされる。

本ツールは purely heuristic。文尾の substantive noun / 数値 / 大文字
固有名詞と、次文冒頭 (40 char) のそれらを比較し、一致 pair 比率を
score 化する。proper-noun 抽出は粗く、最終判断は comments を読み、
実際に読んで追えるかで人間が判定する。

設計方針:
- score は 0-10 の粗い温度計。9 以上は誤差範囲、6 以下は改善サイン。
- PRIMARY OUTPUT は comments (人間レビュアー風の助言文)。
- observations は raw measurement。自動化 numeric 判定には使わない。
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# Language detection (T1 と同等 — _shared.language 未整備のため local)
# ------------------------------------------------------------------
def _detect_language(text: str) -> str:
    """CJK 比率 > 0.3 を "ja"、それ以外を "en" とする。"""
    if not text:
        return "en"
    non_ws = [c for c in text if not c.isspace()]
    if not non_ws:
        return "en"
    cjk = 0
    for c in non_ws:
        o = ord(c)
        if (0x3040 <= o <= 0x30FF) or (0x4E00 <= o <= 0x9FFF):
            cjk += 1
    return "ja" if (cjk / len(non_ws)) > 0.3 else "en"


def _is_japanese(text: str) -> bool:
    return _detect_language(text) == "ja"


# ------------------------------------------------------------------
# Paragraph / sentence splitters
# ------------------------------------------------------------------
_PARAGRAPH_SPLIT = re.compile(r"\n\s*\n+")

_SENT_SPLIT_EN = re.compile(r"(?<=[.!?])\s+")
_SENT_SPLIT_JA = re.compile(r"(?<=[。．！？])")


def _split_paragraphs(text: str) -> list[str]:
    if not text:
        return []
    return [p.strip() for p in _PARAGRAPH_SPLIT.split(text) if p.strip()]


def _split_sentences(paragraph: str, lang: str) -> list[str]:
    paragraph = paragraph.strip()
    if not paragraph:
        return []
    if lang == "ja":
        parts = _SENT_SPLIT_JA.split(paragraph)
    else:
        parts = _SENT_SPLIT_EN.split(paragraph)
    return [s.strip() for s in parts if s.strip()]


# ------------------------------------------------------------------
# Entity extractors
# ------------------------------------------------------------------
# EN: 数字 / 大文字から始まる 2 語以上の proper noun 連 / 単一 Capitalized 語
_DIGIT_RE = re.compile(r"\d+(?:\.\d+)?")
_PROPER_NOUN_RE_EN = re.compile(
    r"\b(?:[A-Z][A-Za-z0-9\-]{1,}(?:\s+[A-Z][A-Za-z0-9\-]{1,}){0,3})\b"
)
# 簡易 stopword (sentence 先頭の Capitalized word が意味的 new-info で
# ないケースを除外する)
_STOPWORDS_EN = {
    "the", "a", "an", "this", "that", "these", "those", "we", "our",
    "their", "it", "its", "he", "she", "they", "i", "you",
    "in", "on", "at", "by", "for", "of", "to", "with", "from", "as",
    "and", "or", "but", "so", "then", "thus", "however", "moreover",
    "furthermore", "finally", "first", "second", "third", "therefore",
    "because", "since", "if", "although", "though", "while", "when",
    "is", "are", "was", "were", "be", "been", "being", "has", "have",
    "had", "do", "does", "did", "can", "could", "may", "might", "must",
    "should", "would", "will", "shall",
}

# JA: Katakana 2 字以上の連、Kanji 2 字以上の連
_KATAKANA_RE = re.compile(r"[゠-ヿー]{2,}")
_KANJI_RE = re.compile(r"[一-鿿]{2,}")


def _strip_latex(text: str) -> str:
    """LaTeX コマンドを削除して素の proper noun を見やすくする。"""
    text = re.sub(r"\\[a-zA-Z]+\*?(?:\{[^}]*\})?", " ", text)
    text = re.sub(r"[\{\}]", " ", text)
    return text


def _extract_entities_en(text: str) -> set[str]:
    """EN: proper-noun 相当 + 数字を抽出。stopword-only な hit は除外。"""
    cleaned = _strip_latex(text)
    hits: set[str] = set()
    for m in _PROPER_NOUN_RE_EN.finditer(cleaned):
        phrase = m.group(0).strip()
        # 全 token が stopword のみ → skip
        tokens = [t.lower() for t in re.split(r"\s+", phrase) if t]
        if not tokens:
            continue
        # 冒頭が Capitalized だが実は stopword (文頭の "The", "We" 等) は除外
        if all(t in _STOPWORDS_EN for t in tokens):
            continue
        # 単一の単語で stopword なら除外
        if len(tokens) == 1 and tokens[0] in _STOPWORDS_EN:
            continue
        hits.add(phrase)
    for m in _DIGIT_RE.finditer(cleaned):
        hits.add(m.group(0))
    return hits


def _extract_entities_ja(text: str) -> set[str]:
    """JA: Katakana compound / Kanji compound / 数字を抽出。"""
    cleaned = _strip_latex(text)
    hits: set[str] = set()
    for m in _KATAKANA_RE.finditer(cleaned):
        hits.add(m.group(0))
    for m in _KANJI_RE.finditer(cleaned):
        hits.add(m.group(0))
    for m in _DIGIT_RE.finditer(cleaned):
        hits.add(m.group(0))
    return hits


def _extract_entities(text: str, lang: str) -> set[str]:
    if lang == "ja":
        return _extract_entities_ja(text)
    return _extract_entities_en(text)


# ------------------------------------------------------------------
# Pair analysis
# ------------------------------------------------------------------
def _pair_ok(s1: str, s2: str, lang: str) -> tuple[bool, str, str, list[str]]:
    """(ok, tail, head, shared_entities)。

    - tail = s1 の末尾 80 char
    - head = s2 の先頭 40 char
    - tail から抽出した entity のうち、head に 1 つでも出てくれば ok=True
    """
    tail = s1[-80:] if len(s1) > 80 else s1
    head = s2[:40] if len(s2) > 40 else s2
    tail_ents = _extract_entities(tail, lang)
    head_ents = _extract_entities(head, lang)
    # head を substring としても許容 (分かち書き違いを吸収)
    shared: list[str] = []
    for e in tail_ents:
        if e in head_ents or e in head:
            shared.append(e)
    return (bool(shared), tail, head, shared)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def paper_writing_given_new_ordering(
    text: str,
    language: str = "auto",
) -> dict:
    """Wallwork §3.4-3.6: Given-New 情報配置ルールの heuristic 評価。

    各文末尾の new info (numbers / 固有名詞) が次文冒頭で再登場
    (given) しているかを確認する。英語論文トップ誌で最も効く
    linguistic technique (skill.md §D)。

    Args:
        text: 診断対象本文。paragraph は空行区切り。
        language: "auto" / "en" / "ja"。auto は CJK 比率で判定。

    Returns:
        dict: score / score_max / comments / observations / hint / source
    """
    if text is None:
        text = ""

    # --- language ---
    if language == "auto":
        lang = "ja" if _is_japanese(text) else "en"
    elif language in ("en", "ja"):
        lang = language
    else:
        lang = "ja" if _is_japanese(text) else "en"

    # --- split ---
    paragraphs = _split_paragraphs(text)

    total_pairs = 0
    ok_pairs = 0
    poor_pairs: list[dict] = []
    pair_idx = 0

    for para in paragraphs:
        sents = _split_sentences(para, lang)
        if len(sents) < 2:
            continue
        for i in range(len(sents) - 1):
            s1 = sents[i]
            s2 = sents[i + 1]
            total_pairs += 1
            ok, tail, head, _shared = _pair_ok(s1, s2, lang)
            if ok:
                ok_pairs += 1
            else:
                if len(poor_pairs) < 5:
                    poor_pairs.append({
                        "pair_idx": pair_idx,
                        "s1_tail": tail,
                        "s2_head": head,
                        "reason": (
                            "tail entity が次文冒頭 40 char 内で再登場せず"
                        ),
                    })
            pair_idx += 1

    ok_ratio = (ok_pairs / total_pairs) if total_pairs else 0.0

    observations = {
        "language": lang,
        "total_sentence_pairs": total_pairs,
        "ok_pairs": ok_pairs,
        "poor_pairs": poor_pairs,
        "ok_ratio": round(ok_ratio, 3),
    }

    # --- score ---
    if total_pairs == 0:
        score = 5.0
    else:
        score = round(ok_ratio * 10, 1)
        # safety clamp
        score = max(0.0, min(10.0, score))

    # --- comments (PRIMARY) ---
    comments: list[str] = []

    if total_pairs == 0:
        comments.append(
            "単文のみ - given-new は無関係。"
            "複数文からなる段落で再診断を。"
        )
        comments.append(
            "Wallwork §3.4-3.6 は「前文末尾の new info を次文冒頭で "
            "given として受け直せ」と指示する。この診断は "
            "paragraph 内 pair に対してのみ有意。"
        )
    else:
        pct = int(round(ok_ratio * 100))
        if ok_ratio < 0.3:
            first_poor = poor_pairs[0] if poor_pairs else None
            if first_poor is not None:
                comments.append(
                    f"Given-New 連結が弱い ({pct}%)。"
                    "読者は各文で出てきた新情報を次文で追いかけられない "
                    "(Wallwork §3.4)。"
                    f"例: {first_poor['s1_tail']} → {first_poor['s2_head']}"
                )
            else:
                comments.append(
                    f"Given-New 連結が弱い ({pct}%)。"
                    "読者は各文で出てきた新情報を次文で追いかけられない "
                    "(Wallwork §3.4)。"
                )
            comments.append(
                "各文末尾に key noun / 数値を置き、次文冒頭で同じ語"
                "(または代名詞・指示詞) で受け直す "
                "(skill.md §D 'Given-New 情報配置')。"
            )
            comments.append(
                "典型修正: \"X yields Y. Z depends on W.\" → "
                "\"X yields Y. This Y, in turn, determines Z, which...\""
            )
        elif ok_ratio < 0.6:
            n = total_pairs - ok_pairs
            comments.append(
                f"Given-New 連結が中程度 ({pct}%)。"
                f"段落内の論理が追いにくい箇所 {n} 件を手直し推奨。"
            )
            if poor_pairs:
                fp = poor_pairs[0]
                comments.append(
                    f"最初の断絶箇所: \"...{fp['s1_tail']}\" → "
                    f"\"{fp['s2_head']}...\"。"
                    "前文末尾の term を次文頭で明示的に受けると連結が回復する。"
                )
            else:
                comments.append(
                    "前文末尾の term を次文頭で明示的に受けると連結が回復する "
                    "(skill.md §D)。"
                )
        else:
            comments.append(
                f"Given-New 連結は良好 ({pct}%)。"
                "新情報が次文の起点として機能している "
                "(Wallwork §3.4-3.6)。"
            )
            if poor_pairs:
                fp = poor_pairs[0]
                comments.append(
                    f"残り {len(poor_pairs)} 件の弱連結はレビューで確認を。"
                    f"例: \"...{fp['s1_tail']}\" → \"{fp['s2_head']}...\"."
                )
            else:
                comments.append(
                    "段落内の theme-rheme 連鎖は保たれている "
                    "(Halliday theme-rheme)。"
                )

    # 2-5 個担保
    if len(comments) < 2:
        comments.append(
            "score は heuristic。proper noun 抽出は粗い。"
            "最終判断は実際に読んで追えるかで。"
        )
    comments = comments[:5]

    hint = (
        "purely heuristic tool。keyword match base で proper noun 抽出は粗い。"
        "最終判断は comments を読み、実際に読んで追えるかで判定。"
    )

    source = (
        "Wallwork §3.4-3.6 / skill.md §D 'Given-New 情報配置' / "
        "Halliday theme-rheme"
    )

    return {
        "score": float(score),
        "score_max": 10,
        "comments": comments,
        "observations": observations,
        "hint": hint,
        "source": source,
    }
