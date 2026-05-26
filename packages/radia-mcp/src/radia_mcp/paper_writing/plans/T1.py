"""Plan B T1: Abstract 強度診断 (paper_writing_abstract_strength).

Wallwork §13.16 / skill.md §A §J に基づく粗い温度計。
Abstract の 4 要素 (problem / method / result-with-number / impact) に
加え、strong-verb・定量数値・banned-phrase を観測する。

設計方針:
- score は 0-10 の粗い温度計。9 以上は誤差範囲、6 以下は改善サイン。
- PRIMARY OUTPUT は comments (人間レビュアー風の助言文)。
- observations は raw measurement。自動化した numeric な正誤判定には使わない。
- paper_writing_check_abstract_background_ratio と補完関係。
  (あちらは ratio の精密測定、本ツールは 4 要素 + 禁句の総合診断)
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# Keyword patterns
# ------------------------------------------------------------------
# 問題提起 (problem / challenge / gap)
_PROBLEM_PATTERNS_EN = [
    r"\bproblem(s)?\b",
    r"\bchallenge(s|d|ing)?\b",
    r"\black(s|ing|ed)?\b",
    r"\bgap(s)?\b",
    r"\blimitation(s)?\b",
    r"\bdifficult(y|ies)?\b",
    r"\bsuffer(s|ed|ing)?\s+from\b",
    r"\bfail(s|ed|ure)?\b",
    r"\bbottleneck(s)?\b",
]
_PROBLEM_PATTERNS_JA = [
    r"課題",
    r"問題",
    r"限界",
    r"困難",
    r"不十分",
    r"ボトルネック",
]

# 手法提示
_METHOD_PATTERNS_EN = [
    r"\bpropose(s|d)?\b",
    r"\bpresent(s|ed)?\b",
    r"\bintroduce(s|d)?\b",
    r"\bdevelop(s|ed)?\b",
    r"\bdesign(s|ed)?\b",
    r"\bconstruct(s|ed)?\b",
    r"\bformulate(s|d)?\b",
]
_METHOD_PATTERNS_JA = [
    r"提案",
    r"開発",
    r"構築",
    r"設計",
    r"定式化",
]

# 数値 + 単位 (% 含む) — 定量結果
# 例: 94.2%, 3.2x, 50 kHz, 120 GB, 1.5 h, 0.1 mm, 12 dB, 1e-3
_UNIT_PATTERN = (
    r"(?:%|×|x\b|‐times|倍|"
    r"[kKMGT]?(?:Hz|B|bps|eV|V|A|W|J|N|Pa|Wh|bar)|"
    r"(?:n|u|µ|m|c|k|M|G|T)?m\b|"
    r"(?:n|u|µ|m)?s\b|"
    r"(?:m|u|µ|n)?g\b|"
    r"dB|deg|°|rad|"
    r"Ω|ohm|"
    r"sec(?:ond)?s?\b|min(?:ute)?s?\b|hour(?:s)?\b|h\b)"
)
_NUM_UNIT_PATTERN = (
    r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?\s*" + _UNIT_PATTERN
)
_NUM_PERCENT_PATTERN = r"\d+(?:\.\d+)?\s*%"
_NUM_X_PATTERN = r"\d+(?:\.\d+)?\s*(?:×|x\b|倍)"

# インパクト / 応用
_IMPACT_PATTERNS_EN = [
    r"\benable(s|d|ing)?\b",
    r"\bimprove(s|d|ment)?\b",
    r"\breduce(s|d|ction)?\b",
    r"\bincrease(s|d)?\b",
    r"\boutperform(s|ed|ing)?\b",
    r"\bachieve(s|d)?\b",
    r"\boffer(s|ed)?\b",
    r"\bmake(s)?\s+\w+\s+feasible\b",
    r"\bspeedup(s)?\b",
]
_IMPACT_PATTERNS_JA = [
    r"実現",
    r"向上",
    r"削減",
    r"達成",
    r"可能にす",
    r"有用",
    r"寄与",
]

# Strong opening verbs (能動・主張動詞)
_STRONG_OPENING_VERBS_EN = [
    r"\bwe\s+propose\b",
    r"\bwe\s+present\b",
    r"\bwe\s+show\b",
    r"\bwe\s+demonstrate\b",
    r"\bwe\s+introduce\b",
    r"\bwe\s+develop\b",
    r"\bwe\s+report\b",
    r"\bwe\s+derive\b",
    r"\bthis\s+paper\s+proposes\b",
    r"\bthis\s+paper\s+presents\b",
    r"\bthis\s+paper\s+introduces\b",
    r"\bthis\s+work\s+proposes\b",
    r"\bthis\s+work\s+presents\b",
]
_STRONG_OPENING_VERBS_JA = [
    r"提案する",
    r"提案し",
    r"示す",
    r"示し",
    r"構築する",
    r"構築し",
    r"開発する",
    r"開発し",
]

# Banned phrases (Wallwork §13.2 / skill.md §A)
_BANNED_PHRASES_EN = [
    "It is well-known",
    "It is well known",
    "In this paper",
    "Recent advances",
    "There are many",
]
_BANNED_PHRASES_JA = [
    "本稿では",  # 冒頭 1 文目に限定して判定
    "近年",
]

# Background / review の signal 語 (sentence 冒頭判定用)
_BACKGROUND_OPENERS_EN = [
    r"^\s*(?:in\s+)?recent(?:ly)?\b",
    r"^\s*(?:over\s+)?the\s+(?:past|last)\s+\w+\s+years\b",
    r"^\s*it\s+(?:is|has\s+been)\s+(?:well[-\s]?known|widely)\b",
    r"^\s*(?:many|several|various|numerous)\s+",
    r"^\s*there\s+(?:are|have\s+been)\s+",
    r"^\s*traditionally\b",
    r"^\s*historically\b",
    r"^\s*conventional(?:ly)?\b",
    r"^\s*in\s+the\s+(?:field|area|literature)\b",
    r"^\s*previous(?:ly)?\s+(?:studies|works?|research)\b",
]
_BACKGROUND_OPENERS_JA = [
    r"^\s*近年",
    r"^\s*従来",
    r"^\s*一般に",
    r"^\s*これまで",
    r"^\s*広く",
    r"^\s*多くの研究",
    r"^\s*先行研究",
]


# ------------------------------------------------------------------
# Helpers
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
        # Hiragana 0x3040-0x309F, Katakana 0x30A0-0x30FF,
        # CJK Unified Ideographs 0x4E00-0x9FFF
        o = ord(c)
        if (0x3040 <= o <= 0x30FF) or (0x4E00 <= o <= 0x9FFF):
            cjk += 1
    return "ja" if (cjk / len(non_ws)) > 0.3 else "en"


def _split_sentences(text: str, lang: str) -> list[str]:
    """. / 。/ ． で文を分割。"""
    if lang == "ja":
        parts = re.split(r"[。．]", text)
    else:
        # EN: 文末句点 + 次の文字が空白 or 終端
        parts = re.split(r"(?<=[.!?])\s+", text)
    return [s.strip() for s in parts if s.strip()]


def _any_match(patterns: list[str], text: str, flags: int = 0) -> bool:
    return any(re.search(p, text, flags) for p in patterns)


def _has_quantitative(text: str) -> bool:
    """数値 + 単位 / % / × の組み合わせを検出。"""
    if re.search(_NUM_PERCENT_PATTERN, text):
        return True
    if re.search(_NUM_X_PATTERN, text):
        return True
    if re.search(_NUM_UNIT_PATTERN, text):
        return True
    return False


def _find_banned_phrases(text: str, lang: str) -> list[str]:
    hits: list[str] = []
    # EN: 全文中どこに出ても拾う (case-insensitive)
    for p in _BANNED_PHRASES_EN:
        if re.search(re.escape(p), text, re.IGNORECASE):
            hits.append(p)
    # JA: 「本稿では」は冒頭 1 文目限定、「近年」も冒頭 1 文目限定
    sentences = _split_sentences(text, lang)
    first = sentences[0] if sentences else ""
    for p in _BANNED_PHRASES_JA:
        if first.startswith(p) or re.match(r"^\s*" + re.escape(p), first):
            hits.append(p)
    return hits


def _has_strong_opening_verb(first_sentence: str, lang: str) -> bool:
    if not first_sentence:
        return False
    if lang == "ja":
        return _any_match(_STRONG_OPENING_VERBS_JA, first_sentence)
    return _any_match(_STRONG_OPENING_VERBS_EN, first_sentence, re.IGNORECASE)


def _background_ratio(sentences: list[str], lang: str) -> float:
    if not sentences:
        return 0.0
    patterns = _BACKGROUND_OPENERS_JA if lang == "ja" else _BACKGROUND_OPENERS_EN
    flags = 0 if lang == "ja" else re.IGNORECASE
    bg = sum(
        1 for s in sentences
        if any(re.search(p, s, flags) for p in patterns)
    )
    return bg / len(sentences)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def paper_writing_abstract_strength(
    abstract: str,
    language: str = "auto",
) -> dict:
    """Abstract の強度を 4 要素 (problem / method / result-with-number / impact)
    + strong-verb + quantitative number + banned-phrase で診断。

    Wallwork §13.16 / skill.md §A §J に基づく粗い温度計。

    Args:
        abstract: Abstract 本文。
        language: "auto" (既定) / "en" / "ja"。auto は CJK 比率 > 0.3 で ja。

    Returns:
        dict: score / score_max / comments / observations / hint / source
    """
    if abstract is None:
        abstract = ""
    text = abstract.strip()

    # --- language detection ---
    if language == "auto":
        lang = _detect_language(text)
    elif language in ("en", "ja"):
        lang = language
    else:
        lang = _detect_language(text)

    sentences = _split_sentences(text, lang)
    first_sentence = sentences[0] if sentences else ""

    # --- observations ---
    if lang == "ja":
        has_problem = _any_match(_PROBLEM_PATTERNS_JA, text)
        has_method = _any_match(_METHOD_PATTERNS_JA, text)
        has_impact = _any_match(_IMPACT_PATTERNS_JA, text)
    else:
        has_problem = _any_match(
            _PROBLEM_PATTERNS_EN, text, re.IGNORECASE
        )
        has_method = _any_match(
            _METHOD_PATTERNS_EN, text, re.IGNORECASE
        )
        has_impact = _any_match(
            _IMPACT_PATTERNS_EN, text, re.IGNORECASE
        )

    has_quantitative_result = _has_quantitative(text)
    strong_opening_verb = _has_strong_opening_verb(first_sentence, lang)
    banned_phrases = _find_banned_phrases(text, lang)
    background_ratio = _background_ratio(sentences, lang)

    observations: dict = {
        "language": lang,
        "has_problem": has_problem,
        "has_method": has_method,
        "has_quantitative_result": has_quantitative_result,
        "has_impact": has_impact,
        "strong_opening_verb": strong_opening_verb,
        "banned_phrases": banned_phrases,
        "background_ratio": round(background_ratio, 3),
    }
    if lang == "en":
        # 空白区切りの単純カウント
        observations["word_count"] = len(text.split()) if text else 0
    else:
        # 非空白文字数 (日本文は句読点含む)
        observations["char_count"] = sum(
            1 for c in text if not c.isspace()
        )

    # --- score ---
    score = 0.0
    if has_problem:
        score += 2
    if has_method:
        score += 2
    if has_quantitative_result:
        score += 2
    if has_impact:
        score += 2
    if strong_opening_verb:
        score += 1
    if len(banned_phrases) == 0:
        score += 1
    if background_ratio > 0.4:
        score -= 1
    score = max(0.0, min(10.0, score))

    # --- comments (PRIMARY) ---
    comments: list[str] = []

    if not has_problem:
        comments.append(
            "Abstract に問題提起が不在。"
            "「Recent X approaches suffer from Y」「現状の Z は W の点で不十分である」"
            "のような課題文を冒頭に。"
        )
    if not has_method:
        comments.append(
            "手法提示が不在。"
            "「We propose X...」「本稿では Y を開発する」の形で"
            "提案内容を明示。"
        )
    if not has_quantitative_result:
        comments.append(
            "定量結果が無い。"
            "「accuracy 94.2%」「reduces latency by 3.2×」等の数値で"
            "効果を裏付ける (Wallwork §13.16 最重要)。"
        )
    if not has_impact:
        comments.append(
            "インパクト文が不在。"
            "「This enables...」「実世界の Z で有用」の形で"
            "applicable domain or downstream benefit を示す。"
        )
    if banned_phrases:
        joined = ", ".join(f"「{p}」" for p in banned_phrases)
        comments.append(
            f"禁句ヒット: {joined} は desk-reject 率を上げる (Wallwork §13.2)。"
            "削る or 書き換え。"
        )
    if not strong_opening_verb and first_sentence:
        comments.append(
            "冒頭動詞が弱い。"
            "受身「...has been studied」や it-is から始めず、"
            "能動動詞 (propose / present / 提案する) で開く。"
        )
    if background_ratio > 0.4:
        comments.append(
            f"Abstract の半分以上 ({background_ratio:.0%}) が既存研究レビュー。"
            "自分の寄与に舵を切る "
            "(2 文以内で背景を済ませ、残りを提案・結果・インパクトへ)。"
        )

    # 全要素揃い + banned 無 → 肯定コメント
    all_elements = (
        has_problem and has_method
        and has_quantitative_result and has_impact
    )
    if all_elements and not banned_phrases:
        comments.append(
            "4 要素 (問題/手法/定量結果/インパクト) 揃い、"
            "banned phrase なし。Abstract として十分。"
        )

    # 短すぎる入力のフォールバック
    if not comments:
        if not text:
            comments.append(
                "Abstract が空。"
                "最低でも 4 文 (問題/手法/結果/意義) を埋めて再診断を。"
            )
            comments.append(
                "skill.md §A §J / Wallwork §13.16 の 4 要素テンプレートを参照。"
            )
        else:
            comments.append(
                "Abstract の判定困難。"
                "skill.md §A §J / Wallwork §13.16 の基準を直接参照し、"
                "人間判断で補完すること。"
            )

    # 最大 6, 最少 2
    comments = comments[:6]
    if len(comments) < 2:
        comments.append(
            "score は粗い指標。"
            "skill.md §J / Wallwork §13 の実例と照合して最終確認を推奨。"
        )

    hint = (
        "score は 4 要素 + 補助指標の粗い合算。"
        "最終的な Abstract 品質は語彙選択・論理の流れが決める。"
        "comments と skill.md §J / Wallwork §13 の実例と照合を。"
    )

    source = (
        "Wallwork §13.16 / skill.md §A §J / "
        "paper_writing_check_abstract_background_ratio と補完関係"
    )

    return {
        "score": float(score),
        "score_max": 10,
        "comments": comments,
        "observations": observations,
        "hint": hint,
        "source": source,
    }
