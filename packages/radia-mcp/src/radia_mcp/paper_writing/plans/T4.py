"""Plan B T4: Discussion 内 Limitation 段落の有無・位置・充実度診断.

paper_writing skill.md §O / NG #4 / Wallwork §9.12 に基づく粗い温度計。
Discussion セクションを検出し、その内部に limitation / caveat 段落が
存在するか、位置 (early / middle / end) と長さ (文字数) を観察する。

設計方針:
- score は 0-10 の粗い温度計。9 以上は誤差範囲、6 以下は改善サイン。
- PRIMARY OUTPUT は comments (人間レビュアー風の助言文)。
- observations は raw measurement。自動化 numeric 判定には使わない。
- skill.md §O: Limitation は Discussion 末尾ではなく中間配置が推奨
  (末尾は防衛的印象、中盤配置は future work への橋渡しに使える)。
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# パターン
# ------------------------------------------------------------------
# Discussion 相当 section 名 (英日、markdown 含む)
_DISCUSSION_SECTION_NAMES = (
    r"Discussion"
    r"|Discussions"
    r"|DISCUSSION"
    r"|DISCUSSIONS"
    r"|考察"
    r"|議論"
    r"|考 察"
)

# \section{Discussion} or \section*{Discussion}
_DISCUSSION_LATEX_PATTERN = re.compile(
    rf"\\section\*?\{{\s*(?:{_DISCUSSION_SECTION_NAMES})\s*\}}",
    re.IGNORECASE,
)

# markdown heading: "## Discussion" or "# Discussion"
_DISCUSSION_MD_PATTERN = re.compile(
    rf"^\s*#{{1,6}}\s+(?:{_DISCUSSION_SECTION_NAMES})\s*$",
    re.IGNORECASE | re.MULTILINE,
)

# 任意の次 section (LaTeX or md)
_ANY_SECTION_LATEX = re.compile(r"\\section\*?\{[^}]+\}")
_ANY_SECTION_MD = re.compile(r"^\s*#{1,6}\s+\S.*$", re.MULTILINE)

# Limitation キーワード (英 + 日)
_LIMITATION_KEYWORDS_EN = [
    r"\blimitations?\b",
    r"\bcaveats?\b",
    r"\bshortcomings?\b",
    r"\bweakness(?:es)?\b",
    r"\bdrawbacks?\b",
]
_LIMITATION_KEYWORDS_JA = [
    r"制約",
    r"限界",
    r"課題",
]

_LIMITATION_PATTERN = re.compile(
    r"(?:"
    + "|".join(_LIMITATION_KEYWORDS_EN + _LIMITATION_KEYWORDS_JA)
    + r")",
    re.IGNORECASE,
)


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def _find_discussion_range(text: str) -> tuple[int, int] | None:
    """Discussion セクションの [start, end] char range を返す。

    LaTeX `\\section{Discussion}` または markdown `## Discussion` を検出。
    範囲は header 直後 (body 先頭) から次の section 見出し直前まで。
    見つからなければ None。
    """
    # LaTeX 優先
    m = _DISCUSSION_LATEX_PATTERN.search(text)
    if m:
        start = m.end()
        next_m = _ANY_SECTION_LATEX.search(text, pos=start)
        end = next_m.start() if next_m else len(text)
        return (start, end)

    # markdown fallback
    m_md = _DISCUSSION_MD_PATTERN.search(text)
    if m_md:
        start = m_md.end()
        # 次の md heading (same or higher level) を探す
        # heuristic: 次の ^#{1,6} heading を探す
        next_md = _ANY_SECTION_MD.search(text, pos=start)
        end = next_md.start() if next_md else len(text)
        return (start, end)

    return None


def _split_paragraphs(body: str) -> list[tuple[int, str]]:
    """空行区切りで段落を分け、(offset_in_body, paragraph_text) のリストを返す。"""
    paragraphs: list[tuple[int, str]] = []
    # 空行で分割しつつ offset を追跡
    cursor = 0
    for m in re.finditer(r"\n\s*\n", body):
        chunk = body[cursor:m.start()]
        if chunk.strip():
            paragraphs.append((cursor, chunk))
        cursor = m.end()
    # tail
    tail = body[cursor:]
    if tail.strip():
        paragraphs.append((cursor, tail))
    return paragraphs


def _classify_position(
    para_offset: int,
    body_length: int,
) -> str:
    """early / middle / end を判定。

    - early: pos < 0.33 * body_length
    - middle: 0.33 <= pos < 0.67
    - end: pos >= 0.67
    """
    if body_length <= 0:
        return "middle"
    ratio = para_offset / body_length
    if ratio < 0.33:
        return "early"
    if ratio < 0.67:
        return "middle"
    return "end"


def _text_head(s: str, n: int = 60) -> str:
    """段落先頭 n 字 (改行は空白に)。"""
    s = re.sub(r"\s+", " ", s).strip()
    if len(s) <= n:
        return s
    return s[:n] + "..."


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def paper_writing_limitation_statement_presence(
    tex_or_text: str,
) -> dict:
    """Discussion 内で limitation/caveat 段落の有無・位置・充実度を検査。

    skill.md §O / Wallwork §9.12: Limitation は Discussion 末尾でなく中間に置く。

    引数:
        tex_or_text: 論文原稿 (LaTeX またはプレーンテキスト / markdown)

    返り値:
        dict: score / score_max / comments / observations / hint / source
    """
    if tex_or_text is None:
        tex_or_text = ""

    # ---- Step 1: Discussion 範囲検出 ----
    disc_range = _find_discussion_range(tex_or_text)
    discussion_found = disc_range is not None

    if discussion_found:
        disc_start, disc_end = disc_range  # type: ignore[misc]
        discussion_body = tex_or_text[disc_start:disc_end]
        discussion_char_range: list[int] | None = [disc_start, disc_end]
    else:
        disc_start, disc_end = 0, 0
        discussion_body = ""
        discussion_char_range = None

    body_length = len(discussion_body)

    # ---- Step 2: limitation 段落検出 ----
    limitation_paragraphs: list[dict] = []
    if discussion_found and body_length > 0:
        for offset_in_body, para in _split_paragraphs(discussion_body):
            if _LIMITATION_PATTERN.search(para):
                position = _classify_position(offset_in_body, body_length)
                # "pos" は document 全体内の絶対位置
                abs_pos = disc_start + offset_in_body
                limitation_paragraphs.append(
                    {
                        "pos": abs_pos,
                        "text_head": _text_head(para, 60),
                        "length_chars": len(para.strip()),
                        "position_in_discussion": position,
                    }
                )

    # ---- Step 3: 集計 ----
    mid_limitation_present = any(
        p["position_in_discussion"] == "middle" for p in limitation_paragraphs
    )
    has_substantial_limitation = any(
        p["length_chars"] >= 150 for p in limitation_paragraphs
    )
    end_only = (
        bool(limitation_paragraphs)
        and all(
            p["position_in_discussion"] == "end" for p in limitation_paragraphs
        )
    )

    observations = {
        "discussion_found": discussion_found,
        "discussion_char_range": discussion_char_range,
        "limitation_paragraphs": limitation_paragraphs,
        "mid_limitation_present": mid_limitation_present,
        "has_substantial_limitation": has_substantial_limitation,
    }

    # ---- Step 4: score ----
    if not discussion_found:
        score = 3.0
    elif not limitation_paragraphs:
        score = 0.0
    else:
        score = 5.0  # base
        if mid_limitation_present:
            score += 3.0
        if has_substantial_limitation:
            score += 2.0
        score = max(0.0, min(10.0, score))

    # ---- Step 5: comments ----
    comments: list[str] = []

    if not discussion_found:
        comments.append(
            "Discussion 見出しが検出されない。"
            "\\section{Discussion} or \\section{考察} を明示。"
            "limitation check はその前提。"
        )
        comments.append(
            "markdown 原稿なら `## Discussion` / `## 考察` を使う。"
            "Discussion が無いと査読者は『結果の解釈が欠けている』と評価する "
            "(skill.md §O)。"
        )
    elif not limitation_paragraphs:
        comments.append(
            "Discussion に limitation/caveat 段落が無い。"
            "査読者は『著者は自分の弱点を理解しているか』を見る。"
            "正直に書いた方が採択率が高い (skill.md NG #4)。"
        )
        comments.append(
            "limitation / caveat / shortcoming / weakness / drawback / "
            "制約 / 限界 / 課題 のいずれかで明示的に段落を起こす。"
        )
    else:
        # limitation 段落あり
        if end_only:
            comments.append(
                "Limitation が Discussion 末尾のみ。末尾は印象が弱い。"
                "Discussion 中盤に置き、その後で対応策・反証可能性を議論する構成 "
                "(Wallwork §9.12)。"
            )

        # 短い limitation 指摘
        short_lengths = [
            p["length_chars"]
            for p in limitation_paragraphs
            if p["length_chars"] < 150
        ]
        if short_lengths and not has_substantial_limitation:
            n = min(short_lengths)
            comments.append(
                f"Limitation 段落が短い ({n} 字)。"
                "トークン制限レベルで留めず、具体的に X / Y / Z の制約を明示。"
            )

        if mid_limitation_present and has_substantial_limitation:
            comments.append(
                "Limitation が Discussion 中盤に配置、長さも十分。"
                "skill.md §O 構成を満たす。"
            )
        elif mid_limitation_present and not has_substantial_limitation:
            comments.append(
                "Limitation は Discussion 中盤にあるが、"
                "150 字未満と薄い。具体的な制約 (データ規模 / 仮定 / "
                "適用範囲) を 1 段落で正直に記述 (skill.md §O)。"
            )
        elif not mid_limitation_present and has_substantial_limitation and not end_only:
            comments.append(
                "Limitation 段落はあるが Discussion 序盤 (early) 配置。"
                "序盤で弱点を出しすぎると主張が埋もれる。"
                "主解釈段落の後、結論段落の前 (中盤) に移す (Wallwork §9.12)。"
            )

    # 2-5 個に調整
    if len(comments) < 2:
        comments.append(
            "score は粗い指標。skill.md §O / Wallwork §9.12 を直接参照して "
            "人間判断で補完すること。"
        )
    comments = comments[:5]

    hint = (
        "\\section 検出は LaTeX 前提。markdown の ## Discussion も "
        "heuristic 対応。Discussion を複数 section に分割している場合 "
        "(Principal findings / Strengths and limitations / ...) は最初の "
        "Discussion 範囲のみ検査する点に注意。"
    )

    return {
        "score": float(score),
        "score_max": 10,
        "comments": comments,
        "observations": observations,
        "hint": hint,
        "source": (
            "skill.md §O NG #4 / Wallwork §9.12 / Reviewer 2 archetype"
        ),
    }
