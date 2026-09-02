"""paper_writing JA-lint re-exports.

Originally these tools lived in mcp-server-document.grant_writing
(LAB-private).  When paper_writing was promoted to radia-mcp (public
PyPI) on 2026-05-26 the 8 JA-lint helpers were INLINE-COPIED into
``_ja_lint.py`` so the public package is self-contained.

The canonical source remains mcp-server-document.grant_writing for
ongoing maintenance of grant-specific JA-lint rules.  Any divergence
between the two copies is a knowledge debt and should be reconciled
by re-extracting from grant_writing whenever the grant copy changes.

This file is the thin wrapper that exposes the 8 helpers under the
``paper_writing_*`` namespace so authors of either grant proposals or
journal papers can call the SAME tool names.
"""
from __future__ import annotations

from .._shared.translationese import check_translationese as _translationese
from ._ja_lint import (
    grant_writing_check_notation_variants as _notation,
    grant_writing_find_undefined_acronyms as _acronyms,
    grant_writing_acronym_usage_audit as _acronym_audit,
    grant_writing_check_kanji_ratio as _kanji_ratio,
    grant_writing_lint_bedrock as _bedrock,
    grant_writing_check_misuse_japanese as _misuse,
    grant_writing_suggest_redundancy_fixes as _redundancy,
    grant_writing_check_subject_predicate_distance as _subj_pred,
)


def paper_writing_check_notation_variants(text: str) -> dict:
    """和文表記ゆれ検出 (grant_writing 実装の re-export)。

    論文和文でも形式名詞ゆれ (事/こと)、送り仮名、カタカナ長音、括弧、
    単位スペース、識別子 case、ハイフン形を統一する必要がある。
    """
    return _notation(text)


def paper_writing_find_undefined_acronyms(
    text: str,
    whitelist: str = "",
    min_len: int = 2,
    max_len: int = 6,
    context_window: int = 200,
) -> dict:
    """Latin 略語の初出定義 check (grant_writing 実装の re-export)。

    論文でも初出 full-form ルールは同じ (特に Abstract / Introduction)。
    """
    return _acronyms(
        text,
        whitelist=whitelist,
        min_len=min_len,
        max_len=max_len,
        context_window=context_window,
    )


def paper_writing_acronym_usage_audit(
    text: str,
    whitelist: str = "",
    min_uses_for_abbrev: int = 3,
    min_len: int = 2,
    max_len: int = 6,
) -> dict:
    """略語の使用頻度監査 (grant_writing 実装の re-export)。

    初出 full-form + 3 回未満なら略語導入せずルールを論文にも適用。
    """
    return _acronym_audit(
        text,
        whitelist=whitelist,
        min_uses_for_abbrev=min_uses_for_abbrev,
        min_len=min_len,
        max_len=max_len,
    )


def paper_writing_check_kanji_ratio(
    text: str,
    min_ratio: float = 0.18,
    max_ratio: float = 0.40,
) -> dict:
    """漢字比率 check (本多『日本語の作文技術』第四章 re-export)。

    和文論文でも 0.25-0.35 が推奨範囲。>0.40 = 官庁漢語感、
    <0.18 = 幼稚・冗長。
    """
    return _kanji_ratio(text, min_ratio=min_ratio, max_ratio=max_ratio)


def paper_writing_lint_bedrock(text: str) -> dict:
    """木下 10 原則 + 本多 + 知的 の bedrock 診断 (re-export)。

    緩衝表現/受身連用/逆茂木型/の連打/一文中「は」2個/文末モノトーン/
    意見事実断定/二重否定 を統合検査。
    """
    return _bedrock(text)


def paper_writing_check_misuse_japanese(text: str) -> dict:
    """『問題な日本語』由来の現代誤用検出 (re-export)。"""
    return _misuse(text)


def paper_writing_suggest_redundancy_fixes(text: str) -> dict:
    """冗長表現 25 パターンの置換候補提示 (re-export)。"""
    return _redundancy(text)


def paper_writing_check_subject_predicate_distance(
    text: str,
    max_chars: int = 40,
) -> dict:
    """主述直結原則 (本多 p.22) — 「は、」「が、」 主題マーカーの
    位置から文末までの距離を計測 (re-export)。"""
    return _subj_pred(text, max_chars=max_chars)


def paper_writing_translationese_check(text: str) -> dict:
    """直訳調・AI調の候補 (自他動詞の誤り、英語注記、定型句、空疎な強調語)。

    和文論文を英語原稿から起こすと「〜を向上する」「劇的に」「〜を可能にする」
    が残りやすい。bedrock・誤用・表記ゆれを通過した原稿にも残る種類の欠陥で、
    2026-09-02 の科研費原稿通読が種になっている (radia_mcp._shared.translationese)。
    """
    return _translationese(text)
