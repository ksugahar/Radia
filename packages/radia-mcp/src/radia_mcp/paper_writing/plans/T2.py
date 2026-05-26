"""Plan B T2: Introduction 末尾 Contribution 明確度診断.

paper_writing skill.md T1 / §C / §P (Swales CARS: Create A Research Space)
に基づく粗い温度計。Introduction 末尾の Contribution リスト / 宣言段落が

    - claim verb で始まっているか
    - 3-5 件の bullet / numbered 列挙か
    - 文法的 parallelism が保たれているか
    - hedge 表現 (might / may / かもしれない 等) が紛れていないか
    - 定量要素 (数値) を含むか

を測定し、score (0-10) と人間レビュアー風 comments を返す。

設計方針:
- score は 0-10 の粗い温度計。9 以上は誤差範囲、6 以下は改善サイン。
- PRIMARY OUTPUT は comments (人間レビュアー風の助言文)。
- observations は raw measurement。自動化 numeric 判定には使わない。
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# パターン
# ------------------------------------------------------------------
# Intro 相当 section 名 (英日)
_INTRO_SECTION_NAMES = (
    r"Introduction"
    r"|Introduction of the paper"
    r"|INTRODUCTION"
    r"|はじめに"
    r"|序論"
    r"|序 論"
    r"|緒言"
)

_SECTION_PATTERN = re.compile(
    r"\\section\*?\{([^}]+)\}",
    re.IGNORECASE,
)

_INTRO_HEAD_PATTERN = re.compile(
    rf"\\section\*?\{{\s*(?:{_INTRO_SECTION_NAMES})\s*\}}",
    re.IGNORECASE,
)

_ITEMIZE_BEGIN_PATTERN = re.compile(
    r"\\begin\{(itemize|enumerate|description)\}"
)
_ITEMIZE_END_PATTERN = re.compile(
    r"\\end\{(itemize|enumerate|description)\}"
)

_ITEM_PATTERN = re.compile(r"\\item\b\s*(.*)")

_CLAIM_VERBS = [
    "propose",
    "present",
    "introduce",
    "derive",
    "prove",
    "demonstrate",
    "show",
    "provide",
    "develop",
    "construct",
    "establish",
]

_HEDGE_WORDS = [
    "may",
    "could",
    "might",
    "potentially",
    "seemingly",
    "arguably",
    "tentatively",
    "かもしれない",
    "可能性がある",
]

# 数値 claim: 12% / 3.2% / 12\% (LaTeX) / 4× / 4x / 3.2x / 10倍 等
_QUANTITATIVE_PATTERN = re.compile(
    r"\d+(?:\.\d+)?\s*\\?%"
    r"|\d+(?:\.\d+)?\s*[×xX](?![a-zA-Z])"
    r"|\d+(?:\.\d+)?\s*(?:倍|fold)"
    r"|\d+(?:\.\d+)?\s*(?:ms|Hz|kHz|MHz|GHz|dB|s|min)"
    r"|\d+(?:\.\d+)?\s*\$?\\times\$?"
)

# 行頭 numbered / sequence marker
_NUMBERED_PREFIX_PATTERN = re.compile(
    r"^\s*(?:"
    r"\d+\s*[.)]"
    r"|第[一二三四五六七]に"
    r"|First[,.]|Second[,.]|Third[,.]|Fourth[,.]"
    r"|Finally[,.]|Moreover[,.]|Furthermore[,.]"
    r")",
    re.MULTILINE,
)

_CLAIM_VERB_RE = re.compile(
    r"\b(?:" + "|".join(_CLAIM_VERBS) + r")\b",
    re.IGNORECASE,
)


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def _find_intro_body(text: str) -> str | None:
    """Introduction セクション本文を返す。見つからなければ None。

    `\\section{Introduction}` から次の `\\section` (any) までを返す。
    """
    m = _INTRO_HEAD_PATTERN.search(text)
    if not m:
        return None
    start = m.end()
    # 次のsectionを探す
    next_m = _SECTION_PATTERN.search(text, pos=start)
    end = next_m.start() if next_m else len(text)
    return text[start:end]


def _extract_itemize_block(body: str) -> tuple[str | None, int]:
    """Intro 本文の中で最初の itemize/enumerate ブロックを返す。

    (block_text, start_pos) を返す。見つからなければ (None, -1)。
    block_text は \\begin{itemize}...\\end{itemize} の本体 (\\item 部分)。
    """
    begin = _ITEMIZE_BEGIN_PATTERN.search(body)
    if not begin:
        return None, -1
    end = _ITEMIZE_END_PATTERN.search(body, pos=begin.end())
    if not end:
        # malformed; take to end of body
        return body[begin.end():], begin.start()
    return body[begin.end():end.start()], begin.start()


def _items_from_itemize(block: str) -> list[str]:
    """`\\item ...` 行から本文のリストを作る。"""
    items: list[str] = []
    # \item で split
    parts = re.split(r"\\item\b", block)
    # 最初の part は \item 前のゴミ
    for p in parts[1:]:
        clean = p.strip()
        if clean:
            items.append(clean)
    return items


def _items_from_numbered_prose(text: str) -> list[str]:
    """番号付き prose (First, Second, 第一に, 1., etc) から項目を抽出。"""
    # find all prefix positions
    positions = [m.start() for m in _NUMBERED_PREFIX_PATTERN.finditer(text)]
    if not positions:
        return []
    items: list[str] = []
    for i, pos in enumerate(positions):
        end = positions[i + 1] if i + 1 < len(positions) else len(text)
        snippet = text[pos:end].strip()
        if snippet:
            items.append(snippet)
    return items


def _first_word(s: str) -> str:
    """先頭の単語 (英語なら alpha, 日本語ならまず 2-3 字) を返す。"""
    s = s.strip()
    # skip leading numbers / punct / backslash commands
    s = re.sub(r"^\\[a-zA-Z]+\s*\{?[^}]*\}?\s*", "", s)
    s = re.sub(r"^\d+\s*[.)]\s*", "", s)
    s = re.sub(r"^第[一二三四五六七]に\s*", "", s)
    s = re.sub(
        r"^(?:First|Second|Third|Fourth|Finally|Moreover|Furthermore)[,.]\s*",
        "",
        s,
    )
    m = re.match(r"[A-Za-z]+", s)
    if m:
        return m.group(0).lower()
    # JP: take 2 chars
    return s[:3]


def _is_parallel(items: list[str]) -> bool:
    """各 item の先頭 POS パターンが揃っているかを粗く判定。

    heuristic:
    - 全て claim verb で始まる → True
    - 全て代名詞 "we" で始まり、次が verb → True
    - 全て冠詞 (a / an / the) で始まる → True
    - 全て同じ名詞句系パターン ("A new", "An efficient", "The proposed") → True
    """
    if len(items) < 2:
        return False

    def _pattern_of(it: str) -> str:
        # strip leading \textbf{} etc
        s = it.strip()
        s = re.sub(r"^\\[a-zA-Z]+\s*\{([^}]*)\}", r"\1", s)
        s = re.sub(r"^\d+\s*[.)]\s*", "", s)
        s = re.sub(r"^第[一二三四五六七]に\s*", "", s)
        s = re.sub(
            r"^(?:First|Second|Third|Fourth|Finally|Moreover|Furthermore)"
            r"[,.]\s*",
            "",
            s,
        )
        s_low = s.lower()
        first = _first_word(s)
        # verb ?
        if first in _CLAIM_VERBS:
            return "VERB"
        # We / Our + ...
        if re.match(r"^\s*(we|our)\b", s_low):
            # We propose / We show → still verb-parallel category "WE_VERB"
            # We followed by a claim verb
            m = re.match(r"^\s*(?:we|our)\s+(\w+)", s_low)
            if m and m.group(1) in _CLAIM_VERBS:
                return "WE_VERB"
            return "WE_OTHER"
        # Article + adj/noun ("A new", "An efficient", "The proposed")
        if re.match(r"^\s*(?:a|an|the)\s+\w+", s_low):
            return "ARTICLE"
        # 日本語始まり
        if re.match(r"[ぁ-んァ-ヶ一-龥]", s):
            return "JP"
        return "OTHER"

    patterns = [_pattern_of(x) for x in items]
    # 全て同じ pattern → parallel
    return len(set(patterns)) == 1


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def paper_writing_contribution_clarity_score(
    tex_or_text: str,
) -> dict:
    """Introduction 末尾の Contribution リスト/段落の明確度診断。

    claim verbs + item count + bullet-parallelism + hedge absence +
    quantitative presence を測定。skill.md T1/§P / Swales CARS に基づく。

    引数:
        tex_or_text: 論文原稿 (LaTeX またはプレーンテキスト)

    返り値:
        dict: score / score_max / comments / observations / hint / source
    """
    if tex_or_text is None:
        tex_or_text = ""

    # ---- Step 1: Intro 抽出 ----
    intro_body = _find_intro_body(tex_or_text)
    fallback_used = False
    if intro_body is None:
        # fallback: 入力全体の最後 300 字
        intro_body = tex_or_text
        fallback_used = True

    # ---- Step 2: Contribution block 抽出 ----
    contribution_block_detected = False
    contribution_text = ""

    # itemize/enumerate 探索
    block, block_start = _extract_itemize_block(intro_body)
    items: list[str] = []

    if block is not None:
        # preceding 100 char paragraph
        preamble_start = max(0, block_start - 100)
        preamble = intro_body[preamble_start:block_start]
        contribution_text = (preamble + "\n" + block).strip()[:500]
        items = _items_from_itemize(block)
        # detection: 存在確認
        contribution_block_detected = True
    else:
        # 末尾 300 字 / 最後の段落
        tail = intro_body[-300:] if len(intro_body) >= 300 else intro_body
        # 最後の段落: 空行で区切られた最後のブロック
        paragraphs = re.split(r"\n\s*\n", intro_body)
        last_para = paragraphs[-1].strip() if paragraphs else ""
        # 末尾 300 字と last_para の短い方 (両方 fallback として結合)
        if len(last_para) > 0 and len(last_para) <= 500:
            contribution_text = last_para[:500]
        else:
            contribution_text = tail.strip()[:500]

        # numbered prose 項目数
        items = _items_from_numbered_prose(contribution_text)

        # contribution-ish signals in tail
        has_signal = False
        signal_patterns = [
            r"contributions?\s+(?:of\s+this\s+paper\s+)?are",
            r"main\s+contributions?",
            r"the\s+contributions?\s+are",
            r"our\s+contributions?",
            r"contributions?\s+are\s+(?:threefold|twofold|fourfold)",
            r"貢献(?:は|として)",
            r"本(?:論文|研究)(?:の)?(?:主な)?貢献",
            r"本(?:論文|研究)(?:は|では)",
        ]
        for p in signal_patterns:
            if re.search(p, contribution_text, re.IGNORECASE):
                has_signal = True
                break

        # claim verb 密度も signal
        cv_matches = _CLAIM_VERB_RE.findall(contribution_text)
        if has_signal or len(cv_matches) >= 2 or len(items) >= 2:
            contribution_block_detected = True

        # fallback 利用時: Intro が全く見つかっていない → detection 弱い
        if fallback_used and not has_signal and len(items) < 2 and len(cv_matches) < 2:
            contribution_block_detected = False

    # ---- Step 3: 観察 ----
    claim_verbs_found: list[str] = []
    if contribution_text:
        seen: set[str] = set()
        for m in _CLAIM_VERB_RE.finditer(contribution_text):
            v = m.group(0).lower()
            if v not in seen:
                seen.add(v)
                claim_verbs_found.append(v)

    item_count = len(items)

    is_parallel_pos = _is_parallel(items) if items else False

    hedge_found: list[str] = []
    lower_text = contribution_text.lower()
    for h in _HEDGE_WORDS:
        # 英単語は word boundary、日本語はそのまま
        if re.match(r"[A-Za-z]", h):
            if re.search(rf"\b{re.escape(h)}\b", lower_text):
                hedge_found.append(h)
        else:
            if h in contribution_text:
                hedge_found.append(h)
    hedge_count = len(hedge_found)

    has_quantitative = bool(_QUANTITATIVE_PATTERN.search(contribution_text))

    observations = {
        "contribution_block_detected": contribution_block_detected,
        "contribution_text": contribution_text,
        "claim_verbs_found": claim_verbs_found,
        "item_count": item_count,
        "is_parallel_pos": is_parallel_pos,
        "hedge_count": hedge_count,
        "has_quantitative": has_quantitative,
    }

    # ---- Step 4: score ----
    if not contribution_block_detected:
        score = 2.0
    else:
        score = 0.0
        # claim verbs: +2 each, max +6
        score += min(len(claim_verbs_found) * 2, 6)
        # item count sweet spot
        if 3 <= item_count <= 5:
            score += 2
        # parallel
        if is_parallel_pos:
            score += 1
        # no hedge
        if hedge_count == 0:
            score += 1
        # quantitative
        if has_quantitative:
            score += 1
        score = max(0.0, min(10.0, score))

    # ---- Step 5: comments ----
    comments: list[str] = []

    if not contribution_block_detected:
        comments.append(
            "Introduction 末尾に Contribution 明示ブロックが見当たらない。"
            "最後の段落で「我々の貢献は以下である: ...」と明示的に宣言する "
            "(skill.md §C)。"
        )
        comments.append(
            "3-5 件の bullet (\\begin{itemize}\\item ...) または "
            "First, Second, Finally, の numbered prose で Contribution を列挙すると、"
            "査読者の最初のページ読みで「何が新しいか」が即座に伝わる。"
        )
    else:
        if hedge_count > 0:
            shown = ", ".join(hedge_found[:3])
            comments.append(
                f"Contribution 内に hedge 表現 ({shown})。"
                "「X を提案する」は断定形で書く。"
                "hedge は Discussion 節で使う (skill.md 原則 8)。"
            )
        if item_count and item_count < 3:
            comments.append(
                f"Contribution 項目が {item_count} 件と少ない。"
                "3-5 件の bullet 列挙が査読者フレンドリ (skill.md §C)。"
            )
        if item_count > 7:
            comments.append(
                f"Contribution が {item_count} 件と多すぎる。"
                "3-5 件に圧縮し、残りは §以降で詳述する。"
            )
        if items and not is_parallel_pos:
            comments.append(
                "Contribution 項目の書き出しが不揃い (文法的 parallelism 崩壊)。"
                "同一 POS で揃える (全項目を verb 開始 or 全項目を noun 開始)。"
            )
        if not has_quantitative:
            comments.append(
                "Contribution に定量要素がない。"
                "「accuracy を 12% 向上」「3.2× faster」のような数値で"
                "引き合いを強化する (skill.md §D reviewer 2 対策)。"
            )
        if len(claim_verbs_found) == 0:
            comments.append(
                "claim verb (propose / demonstrate / show / establish 等) が"
                "検出されない。「We discuss / We consider」は弱い。"
                "各項目で動詞を明示して貢献を断言する (skill.md §C)。"
            )
        elif len(claim_verbs_found) < 2:
            comments.append(
                f"claim verb が {len(claim_verbs_found)} 件のみ ("
                f"{', '.join(claim_verbs_found)})。"
                "3 項目あれば 3 動詞 (propose → demonstrate → show 等) で"
                "貢献を階層的に宣言したい。"
            )

        # 全て OK → 肯定
        if (
            len(claim_verbs_found) >= 3
            and 3 <= item_count <= 5
            and is_parallel_pos
            and hedge_count == 0
            and has_quantitative
        ):
            comments.append(
                "Contribution block 検出、claim verbs・parallelism・"
                "定量要素すべて揃う。"
            )

    # 2-6 個担保
    if len(comments) < 2:
        comments.append(
            "score は粗い指標。skill.md §C / T1 を直接参照して"
            "人間判断で補完すること。"
        )
    comments = comments[:6]

    hint = (
        "score は Introduction 末尾の Contribution 段落を前提とした粗い温度計。"
        "Introduction の構造自体 (CARS: Establish territory → Niche → Fill niche) は"
        "別途 check。9 以上は誤差範囲、6 以下は改善サイン。"
    )

    return {
        "score": float(score),
        "score_max": 10,
        "comments": comments,
        "observations": observations,
        "hint": hint,
        "source": (
            "skill.md T1 §C §P / Swales CARS (Create A Research Space) / "
            "Wallwork §14"
        ),
    }
