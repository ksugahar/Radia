"""弱気修飾 (hedge) の共通パターン集。

木下『理科系の作文技術』原則 8 (緩衝表現削除) に基づき、和文で
曖昧にぼかす表現を集中管理。grant_writing / paper_writing /
presentation で共通利用する。
"""
from __future__ import annotations

import re

# pattern name -> regex string (no capture group for simpler len(findall))
HEDGE_PATTERNS: dict[str, str] = {
    "と考えられる": r"と考えられ(?:る|ます)",
    "と思われる": r"と思われ(?:る|ます)",
    "と見られる": r"と見られ(?:る|ます)",
    "ではないかと": r"ではないかと(?:思われる|思う)?",
    "かもしれない": r"かもしれ(?:ない|ません)",
    "と見てよい": r"と見てよい",
    "ように思う": r"ように思(?:う|われる|います)",
    "可能性がある": r"可能性が(?:ある|あります)",
    "示唆される": r"示唆され(?:る|ます)",
    "期待される": r"期待され(?:る|ます)",
}


def scan_hedges(text: str) -> tuple[int, dict[str, int]]:
    """Return (total_count, per_pattern_counts)."""
    by_pat: dict[str, int] = {}
    total = 0
    for name, pat in HEDGE_PATTERNS.items():
        n = len(re.findall(pat, text))
        if n:
            by_pat[name] = n
            total += n
    return total, by_pat
