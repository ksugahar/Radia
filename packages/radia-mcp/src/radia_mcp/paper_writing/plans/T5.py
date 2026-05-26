"""Plan B T5: Introduction 内 Related Work 密度 / 自己引用 / 年度分布 診断.

paper_writing skill.md NG #5 (Related work が自己引用ばかり) / Wallwork §14 に
基づく粗い温度計。Introduction 節の `\\cite` 密度・自己引用比率・recency (5 年
以内比率) をまとめて観測し、competitor 研究の引用量が薄くないかを測る。

設計方針:
- score は 0-10 の粗い温度計。9 以上は誤差範囲、6 以下は改善サイン。
- PRIMARY OUTPUT は comments (人間レビュアー風の助言文)。
- observations は raw measurement。自動 numeric 判定には使わない。
- tools.py には触れない。bib parser は brace-balanced な local 実装。
"""

from __future__ import annotations

import re


# ------------------------------------------------------------------
# Intro section detection
# ------------------------------------------------------------------
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

_CITE_PATTERN = re.compile(
    r"\\(?:cite|citep|citet|autocite|parencite|citealp|citealt|nocite)"
    r"\*?(?:\[[^\]]*\])?(?:\[[^\]]*\])?\{([^}]+)\}"
)


# Current year anchor (skill.md / spec: current = 2026)
CURRENT_YEAR = 2026
RECENT_YEARS = 5


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def _find_intro_body(text: str) -> str | None:
    """Introduction 本文 (次の \\section まで) を返す。見つからなければ None。"""
    m = _INTRO_HEAD_PATTERN.search(text)
    if not m:
        return None
    start = m.end()
    next_m = _SECTION_PATTERN.search(text, pos=start)
    end = next_m.start() if next_m else len(text)
    return text[start:end]


def _extract_cite_keys(text: str) -> list[str]:
    """text 中の \\cite{...} キーを出現順に抽出。重複含む。"""
    keys: list[str] = []
    for m in _CITE_PATTERN.finditer(text):
        for k in m.group(1).split(","):
            k = k.strip()
            if k:
                keys.append(k)
    return keys


def _extract_balanced(text: str, start: int) -> tuple[int, str]:
    """text[start] == '{' から対応する '}' までを返す。

    Returns (end_index_after_closing_brace, inner).
    ミスマッチは (-1, "")。
    """
    if start >= len(text) or text[start] != "{":
        return -1, ""
    depth = 0
    i = start
    n = len(text)
    while i < n:
        c = text[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i + 1, text[start + 1 : i]
        i += 1
    return -1, ""


def _parse_bib_entries(bib: str) -> dict[str, dict[str, str]]:
    """brace-balanced .bib parser。

    返り値: {cite_key: {"author": "...", "year": "...", "type": "article", ...}}
    """
    entries: dict[str, dict[str, str]] = {}
    if not bib:
        return entries

    for hdr in re.finditer(r"@(\w+)\s*\{", bib):
        brace_start = hdr.end() - 1  # position of '{'
        entry_type = hdr.group(1).lower()
        entry_end, entry_body = _extract_balanced(bib, brace_start)
        if entry_end < 0:
            continue

        key_match = re.match(r"\s*([^,\s]+)\s*,", entry_body)
        if not key_match:
            continue
        key = key_match.group(1)
        fields_region = entry_body[key_match.end():]

        fields: dict[str, str] = {"type": entry_type}

        # field = value 抽出 (author / year / title etc)
        for fm in re.finditer(
            r"([A-Za-z_]+)\s*=\s*", fields_region
        ):
            fname = fm.group(1).lower()
            vstart = fm.end()
            if vstart >= len(fields_region):
                continue
            vc = fields_region[vstart]
            value = ""
            if vc == "{":
                _, inner = _extract_balanced(fields_region, vstart)
                value = inner
            elif vc == '"':
                end_q = fields_region.find('"', vstart + 1)
                if end_q > 0:
                    value = fields_region[vstart + 1 : end_q]
            else:
                # bare value (number / abbrev)
                end_bare = re.search(r"[,\n}]", fields_region[vstart:])
                if end_bare:
                    value = fields_region[vstart : vstart + end_bare.start()]
                else:
                    value = fields_region[vstart:]
            # 既出の field を上書きしない (最初の出現を優先)
            if fname not in fields:
                fields[fname] = value.strip()

        entries[key] = fields

    return entries


def _extract_year(year_field: str) -> int | None:
    """bib の year field から西暦 4 桁を抽出。失敗したら None。"""
    if not year_field:
        return None
    m = re.search(r"\b(1[89]\d{2}|20\d{2}|21\d{2})\b", year_field)
    if not m:
        return None
    return int(m.group(1))


def _author_matches(author_field: str, authors: set[str]) -> bool:
    """author field に author_last_names のいずれかが含まれるかを判定。

    BibTeX の author field は "Last, First and Last2, First2 and ..." 形式または
    "First Last and First2 Last2 and ..." 形式。粗く姓を拾うため case-insensitive
    な word-boundary 照合を行う。
    """
    if not author_field or not authors:
        return False
    # 波括弧 / backslash コマンド除去
    cleaned = re.sub(r"[{}\\]", " ", author_field)
    low = cleaned.lower()
    for last in authors:
        last_low = last.strip().lower()
        if not last_low:
            continue
        # word boundary 照合 (日本語姓も半角 ASCII 扱いで fallback)
        pat = r"(?:^|[^A-Za-z0-9_])" + re.escape(last_low) + r"(?:[^A-Za-z0-9_]|$)"
        if re.search(pat, low):
            return True
        # ASCII で bound が取れない場合 (純 CJK) は substring 検索
        if not re.match(r"^[A-Za-z0-9_].*[A-Za-z0-9_]$", last_low):
            if last_low in low:
                return True
    return False


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def paper_writing_related_work_density(
    tex: str,
    bib: str = "",
    author_last_names: str = "",
) -> dict:
    """Introduction 内の \\cite 密度・自己引用比率・年度分布を診断。

    skill.md NG #5 (related work が自己引用だらけ) / Wallwork §14 に基づく。

    Args:
        tex: 論文原稿 (LaTeX 全文)。
        bib: .bib ファイル内容 (文字列)。year / author 抽出に使用。
        author_last_names: 自分+共著者の姓をカンマ区切りで指定。
            例: "Sugahara, Hane, Matsuo"

    Returns:
        dict: score / score_max / comments / observations / hint / source
    """
    if tex is None:
        tex = ""
    if bib is None:
        bib = ""
    if author_last_names is None:
        author_last_names = ""

    authors_set = {
        a.strip() for a in author_last_names.split(",") if a.strip()
    }

    # ---- Step 1: Introduction 抽出 ----
    intro_body = _find_intro_body(tex)
    introduction_found = intro_body is not None

    # ---- Step 2: cite 抽出 ----
    all_cite_keys = _extract_cite_keys(tex)
    total_cites_in_paper = len(set(all_cite_keys))

    if introduction_found:
        intro_cite_keys_ordered = _extract_cite_keys(intro_body or "")
    else:
        intro_cite_keys_ordered = []

    # unique, preserve order
    seen: set[str] = set()
    introduction_cite_keys: list[str] = []
    for k in intro_cite_keys_ordered:
        if k not in seen:
            seen.add(k)
            introduction_cite_keys.append(k)
    introduction_cite_total = len(introduction_cite_keys)

    # ---- Step 3: bib parse ----
    bib_entries = _parse_bib_entries(bib)

    # ---- Step 4: self cite ----
    self_cite_keys: list[str] = []
    if authors_set:
        for key in introduction_cite_keys:
            entry = bib_entries.get(key)
            if not entry:
                continue
            if _author_matches(entry.get("author", ""), authors_set):
                self_cite_keys.append(key)

    self_cite_ratio = (
        len(self_cite_keys) / introduction_cite_total
        if introduction_cite_total
        else 0.0
    )
    competitor_cite_count = introduction_cite_total - len(self_cite_keys)

    # ---- Step 5: year histogram (intro citations) ----
    year_histogram: dict[int, int] = {}
    recent_count = 0
    year_known_count = 0
    for key in introduction_cite_keys:
        entry = bib_entries.get(key)
        if not entry:
            continue
        y = _extract_year(entry.get("year", ""))
        if y is None:
            continue
        year_histogram[y] = year_histogram.get(y, 0) + 1
        year_known_count += 1
        if y >= CURRENT_YEAR - RECENT_YEARS:
            recent_count += 1

    recent_ratio = (
        recent_count / introduction_cite_total
        if introduction_cite_total
        else 0.0
    )

    observations = {
        "introduction_found": introduction_found,
        "introduction_cite_keys": introduction_cite_keys,
        "introduction_cite_total": introduction_cite_total,
        "self_cite_keys": self_cite_keys,
        "self_cite_ratio": round(self_cite_ratio, 3),
        "competitor_cite_count": competitor_cite_count,
        "year_histogram": year_histogram,
        "recent_ratio": round(recent_ratio, 3),
        "total_cites_in_paper": total_cites_in_paper,
    }

    # ---- Step 6: score ----
    if not introduction_found or introduction_cite_total == 0:
        score = 2.0
    else:
        score = 0.0
        if 20 <= introduction_cite_total <= 40:
            score += 3
        if self_cite_ratio <= 0.20:
            score += 3
        if recent_ratio >= 0.40:
            score += 2
        if competitor_cite_count >= 5:
            score += 2
        score = max(0.0, min(10.0, score))

    # ---- Step 7: comments ----
    comments: list[str] = []

    if not introduction_found:
        comments.append(
            "Introduction 見出しが見つからない。"
            "\\section{Introduction} 明示後に再実行。"
        )
    elif introduction_cite_total == 0:
        comments.append(
            "Introduction 内に \\cite が 0 件。"
            "関連研究引用が皆無では reviewer-2 の「何も読んでいない」懸念を招く "
            "(skill.md NG #5 / Wallwork §14)。"
        )

    if introduction_found and introduction_cite_total > 0:
        if introduction_cite_total < 20:
            comments.append(
                f"Introduction の引用が {introduction_cite_total} 件と少ない "
                "(target: 20-40)。関連研究への配置が薄いと reviewer-2 の"
                "「何も読んでいない」懸念を招く。"
            )
        if introduction_cite_total > 40:
            comments.append(
                f"Introduction に引用 {introduction_cite_total} 件は多すぎ "
                "(target 20-40)。Related Work 節に切り分けて "
                "Intro を compact にする。"
            )
        if self_cite_ratio > 0.20:
            pct = int(round(self_cite_ratio * 100))
            comments.append(
                f"自己引用比率 {pct}% が高め (20% 以下推奨)。"
                "競合研究の引用を増やす (skill.md NG #5)。"
            )
        if recent_ratio < 0.40:
            pct = int(round(recent_ratio * 100))
            comments.append(
                f"5 年以内の引用が {pct}% と少ない (40% 以上推奨)。"
                "最新動向を追えていない印象を与える。"
            )
        if competitor_cite_count < 5:
            comments.append(
                f"競合研究の引用が少ない (競合 {competitor_cite_count} 件)。"
                "「why this vs [17]」系の review comments を誘発 "
                "(archetype Q)。"
            )

        # all good
        if (
            20 <= introduction_cite_total <= 40
            and self_cite_ratio <= 0.20
            and recent_ratio >= 0.40
            and competitor_cite_count >= 5
        ):
            comments.append(
                "Intro 引用密度・自己引用比率・recency バランス、"
                "いずれも良好。"
            )

    # 2-5 comments 担保
    if len(comments) < 2:
        comments.append(
            "score は粗い指標。skill.md NG #5 / Wallwork §14 を直接参照して"
            "人間判断で補完すること。"
        )
    comments = comments[:5]

    hint = (
        "bib parser は brace-balanced (v0.11.1 fix 済み)。"
        "author_last_names = 自分+共著者の姓をカンマ区切り。"
    )

    source = (
        "skill.md NG #5 / Wallwork §14 / "
        "reviewer archetype Q '[17] と比べてなぜ?'"
    )

    return {
        "score": float(score),
        "score_max": 10,
        "comments": comments,
        "observations": observations,
        "hint": hint,
        "source": source,
    }
