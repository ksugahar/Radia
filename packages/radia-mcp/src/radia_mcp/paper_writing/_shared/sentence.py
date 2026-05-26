"""文分割 helper (和文 / mixed ja-en)。

既存コードでは `re.split(r"[。．]", text)` を 8 箇所で別々に呼んでいた。
本 helper に集約し、将来の境界ルール変更 (例: 句読点の除外、カッコ
内句点の扱い) が 1 箇所で完結するように。
"""
from __future__ import annotations

import re

_RE_JA_SPLITTER = re.compile(r"[。．]")
_RE_MIXED_SPLITTER = re.compile(r"[。．\.!?！？]")


def split_ja_sentences(text: str) -> list[str]:
    """和文の文分割 (句点 。／．)。空白のみの断片は除外。"""
    if not text:
        return []
    return [s.strip() for s in _RE_JA_SPLITTER.split(text) if s.strip()]


def split_mixed_sentences(text: str) -> list[str]:
    """和文+英文 mixed の文分割 (。．.!?！？ すべて)。"""
    if not text:
        return []
    return [s.strip() for s in _RE_MIXED_SPLITTER.split(text) if s.strip()]
