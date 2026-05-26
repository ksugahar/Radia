"""日本語判定 helper (CJK 比率ベース)。"""
from __future__ import annotations


def cjk_char_count(text: str) -> int:
    """text 内の CJK 文字 (ひらがな・カタカナ・漢字) の総数を返す。"""
    if not text:
        return 0
    count = 0
    for c in text:
        if (
            "぀" <= c <= "ゟ"  # Hiragana
            or "゠" <= c <= "ヿ"  # Katakana
            or "一" <= c <= "鿿"  # CJK Unified Ideographs
        ):
            count += 1
    return count


def is_japanese(text: str, threshold: float = 0.3) -> bool:
    """非空白 char 中の CJK 比率 >= threshold なら日本語と判定。

    0.3 は経験値 (Wallwork §13 で英語 abstract 判定の閾値として採用、
    skill.md と既存 paper_writing の validate_abstract_length と整合)。
    """
    if not text:
        return False
    total = sum(1 for c in text if not c.isspace())
    if total == 0:
        return False
    return cjk_char_count(text) / total >= threshold
