"""presentation 用の JA-lint re-exports (v0.13.0).

grant_writing に実装されている和文 lint をスライド原稿 (台本) にも
適用できるよう、presentation 名前空間に露出する薄い wrapper。
pptx 埋め込みテキスト or speaker note を抽出後、そのまま feed する
運用を想定。
"""
from __future__ import annotations

from .._shared.translationese import check_translationese as _translationese
from ._text_checks import (
    grant_writing_check_notation_variants as _notation,
    grant_writing_find_undefined_acronyms as _acronyms,
    grant_writing_acronym_usage_audit as _acronym_audit,
    grant_writing_check_kanji_ratio as _kanji_ratio,
    grant_writing_lint_bedrock as _bedrock,
    grant_writing_check_misuse_japanese as _misuse,
    grant_writing_suggest_redundancy_fixes as _redundancy,
)


def presentation_check_notation_variants(text: str) -> dict:
    """スライドテキストの表記ゆれ検出 (re-export)。
    extract_pptx_text の返り値を feed する運用。"""
    return _notation(text)


def presentation_find_undefined_acronyms(
    text: str,
    whitelist: str = "",
    min_len: int = 2,
    max_len: int = 6,
    context_window: int = 200,
) -> dict:
    """スライド内略語の初出定義 check (re-export)。"""
    return _acronyms(
        text,
        whitelist=whitelist,
        min_len=min_len,
        max_len=max_len,
        context_window=context_window,
    )


def presentation_acronym_usage_audit(
    text: str,
    whitelist: str = "",
    min_uses_for_abbrev: int = 3,
    min_len: int = 2,
    max_len: int = 6,
) -> dict:
    """略語使用頻度監査 (re-export)。"""
    return _acronym_audit(
        text,
        whitelist=whitelist,
        min_uses_for_abbrev=min_uses_for_abbrev,
        min_len=min_len,
        max_len=max_len,
    )


def presentation_check_kanji_ratio(
    text: str,
    min_ratio: float = 0.18,
    max_ratio: float = 0.40,
) -> dict:
    """スライド台本の漢字比率 check (re-export)。
    口頭発表台本はやや低め (0.22 前後) の方が朗読しやすい傾向。"""
    return _kanji_ratio(text, min_ratio=min_ratio, max_ratio=max_ratio)


def presentation_lint_bedrock(text: str) -> dict:
    """台本・スライド注釈の bedrock lint (木下 10 原則、re-export)。"""
    return _bedrock(text)


def presentation_check_misuse_japanese(text: str) -> dict:
    """台本の現代誤用検出 (re-export)。"""
    return _misuse(text)


def presentation_suggest_redundancy_fixes(text: str) -> dict:
    """冗長表現 25 パターンの置換候補 (re-export)。"""
    return _redundancy(text)


def presentation_translationese_check(text: str) -> dict:
    """台本・スライド本文の直訳調・AI調の候補 (自他動詞の誤り、英語注記、定型句)。

    2026-09-02 の科研費原稿通読で、bedrock・誤用・表記ゆれをすべて通過した
    文に「〜へ発展する」「接着層（glue）」「劇的な差」が残っていた実例が種。
    extract_pptx_text や speaker note の返り値をそのまま feed する。
    """
    return _translationese(text)
