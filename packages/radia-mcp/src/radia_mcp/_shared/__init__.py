"""共有 helpers (v0.13.0 から導入)。

全 writing-family モジュール (grant_writing / paper_writing / presentation)
で共通利用する低レベル helper を集約。モジュール間で regex パターンや
判定ロジックが diverge するのを防ぐ。

- hedges:    弱気修飾 (と考えられる / 可能性がある 等) の集中管理
- language:  日本語テキスト判定 (CJK 比率)
- sentence:  文分割 (句点 ．。, mixed ja/en)
"""
from .hedges import HEDGE_PATTERNS, scan_hedges  # noqa: F401
from .language import is_japanese, cjk_char_count  # noqa: F401
from .sentence import split_ja_sentences, split_mixed_sentences  # noqa: F401
