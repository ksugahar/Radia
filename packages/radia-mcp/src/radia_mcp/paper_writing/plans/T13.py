"""Plan B T13: Citation health 4 axes
(paper_writing_citation_health_4_axes).

論文の reference list を 4 軸で診断:
  1. **recency** — 引用の median age が古すぎないか (out-of-date literature)
  2. **author_concentration** — 同一著者への引用集中 (echo chamber)
  3. **geographic_diversity** — 著者の地理分散 (英語圏 only / 日本のみ)
  4. **self_citation_ratio** — 自己引用比率 (既存 T 番号と相補的に細分化)

入力: bib 文字列 (BibTeX) または references セクション抜粋。
パース: BibTeX エントリ (year / author / address / journal) を粗く抽出。

設計方針:
- score は 4 軸の coverage 率の温度計 (0-10)。
- 各軸は warning threshold を持ち、外れると flag。
- bib 形式が混在しても regex で粗く取れる範囲で動作。
- 既存 paper_writing_check_self_citation_ratio と組み合わせて使う。
"""

from __future__ import annotations

import datetime as _dt
import re
import statistics
from collections import Counter


def _parse_bib_entries(bib: str) -> list[dict]:
    """Parse BibTeX via the package-wide balanced-brace parser."""
    from ...bibliography._bibparse import parse_bib

    entries: list[dict] = []
    for parsed in parse_bib(bib):
        if not parsed.key:
            continue
        fields = parsed.fields
        entry: dict = {}
        # year
        year_str = fields.get("year", "")
        ym = re.search(r"\d{4}", year_str)
        if ym:
            try:
                entry["year"] = int(ym.group(0))
            except ValueError:
                pass
        # author
        if "author" in fields:
            authors = re.split(r"\s+and\s+", fields["author"])
            entry["authors"] = [a.strip() for a in authors if a.strip()]
        # journal / booktitle
        entry["journal"] = fields.get("journal", "") or fields.get(
            "booktitle", ""
        )
        # address (城址 detection)
        entry["address"] = fields.get("address", "") or fields.get(
            "publisher", ""
        )
        entries.append(entry)
    return entries


# Geographic mapping: very rough, based on country/city tokens in journal names + address
_REGION_TOKENS = {
    "USA": ["USA", "United States", "California", "Boston", "MIT",
             "Stanford", "Harvard", "Princeton", "Berkeley", "MA",
             "CA", "NY", "American"],
    "Europe": ["UK", "Britain", "London", "Cambridge", "Oxford",
                "Germany", "France", "Italy", "Switzerland",
                "Netherlands", "Spain", "Sweden", "Denmark",
                "European", "EPL"],
    "Japan": ["Japan", "Tokyo", "Osaka", "Kyoto", "Nagoya", "Sendai",
               "Sapporo", "Tsukuba", "RIKEN", "JST", "JSPS", "IEEJ",
               "Japanese", "東京", "大阪", "京都"],
    "China": ["China", "Beijing", "Shanghai", "Tsinghua", "Peking",
                "Chinese", "Hong Kong"],
    "Other_Asia": ["Korea", "Seoul", "Taiwan", "Taipei", "Singapore",
                    "India", "Mumbai", "Delhi"],
}


def _region_of(entry: dict) -> str:
    haystack = " ".join([
        entry.get("journal", ""),
        entry.get("address", ""),
        " ".join(entry.get("authors", [])),
    ])
    for region, tokens in _REGION_TOKENS.items():
        for t in tokens:
            if re.search(r"\b" + re.escape(t) + r"\b", haystack,
                         re.IGNORECASE):
                return region
    return "Unknown"


def _last_name(author: str) -> str:
    """Extract surname. BibTeX format: 'Lastname, Firstname' or 'Firstname Lastname'."""
    if "," in author:
        return author.split(",")[0].strip()
    parts = author.strip().split()
    return parts[-1] if parts else ""


def paper_writing_citation_health_4_axes(
    bib: str,
    author_last_names: str = "",
    current_year: int = 0,
    text_for_inline_check: str = "",
) -> dict:
    """論文 reference の 4 軸 health 診断。

    Args:
        bib: .bib ファイル内容 (BibTeX)
        author_last_names: 自分+共著者の姓 (カンマ区切り、self-citation 用)
        current_year: 「現在年」、0 なら今日の year を使用
        text_for_inline_check: 本文 (\\cite{} カウントが必要なら)

    Returns:
        dict: score / per_axis (recency / author_concentration /
              geographic_diversity / self_citation) / comments
    """
    if bib is None:
        bib = ""
    if current_year == 0:
        current_year = _dt.datetime.now().year

    entries = _parse_bib_entries(bib)
    n_entries = len(entries)

    if n_entries == 0:
        return {
            "score": None,
            "comments": [
                "BibTeX エントリが 1 件も parse できなかった。"
                "正しい .bib 形式 (@article{key, year={2024}, ...}) で渡されているか確認。"
            ],
            "hint": "input が空/不正。",
            "source": "T13 citation health (paper-writing 2026-04)",
        }

    # ---- 1. recency ----
    years = [e["year"] for e in entries if "year" in e]
    if years:
        median_year = statistics.median(years)
        median_age = current_year - median_year
        recency_status = (
            "FRESH" if median_age <= 5
            else ("MODERATE" if median_age <= 10 else "OUT_OF_DATE")
        )
    else:
        median_year = None
        median_age = None
        recency_status = "UNKNOWN"

    recency_axis = {
        "median_year": median_year,
        "median_age_years": median_age,
        "status": recency_status,
        "warning": median_age is not None and median_age > 10,
        "year_distribution": {
            "n_within_5yr": sum(
                1 for y in years if current_year - y <= 5
            ),
            "n_5_to_10yr": sum(
                1 for y in years if 5 < current_year - y <= 10
            ),
            "n_over_10yr": sum(
                1 for y in years if current_year - y > 10
            ),
        },
    }

    # ---- 2. author_concentration ----
    all_first_authors: list[str] = []
    for e in entries:
        if "authors" in e and e["authors"]:
            all_first_authors.append(_last_name(e["authors"][0]))
    author_counter = Counter(all_first_authors)
    most_common_5 = author_counter.most_common(5)
    top_share = (
        most_common_5[0][1] / n_entries if most_common_5 else 0
    )
    concentration_status = (
        "DIVERSE" if top_share <= 0.10
        else ("MODERATE" if top_share <= 0.20 else "ECHO_CHAMBER")
    )
    concentration_axis = {
        "top5_first_authors": most_common_5,
        "top1_share": round(top_share, 3),
        "unique_first_authors": len(author_counter),
        "status": concentration_status,
        "warning": top_share > 0.20,
    }

    # ---- 3. geographic_diversity ----
    regions = [_region_of(e) for e in entries]
    region_counter = Counter(regions)
    n_regions_known = sum(
        v for k, v in region_counter.items() if k != "Unknown"
    )
    region_status = (
        "DIVERSE" if len(
            [k for k, v in region_counter.items()
             if k != "Unknown" and v >= 2]
        ) >= 3
        else ("MODERATE" if n_regions_known >= 2 else "REGION_BIAS")
    )
    geographic_axis = {
        "region_distribution": dict(region_counter),
        "n_regions_with_2plus": len(
            [k for k, v in region_counter.items()
             if k != "Unknown" and v >= 2]
        ),
        "status": region_status,
        "warning": region_status == "REGION_BIAS",
    }

    # ---- 4. self_citation ----
    self_authors = {
        s.strip().lower() for s in author_last_names.split(",")
        if s.strip()
    }
    n_self = 0
    if self_authors:
        for e in entries:
            for a in e.get("authors", []):
                ln = _last_name(a).lower()
                if ln in self_authors:
                    n_self += 1
                    break
    self_ratio = n_self / n_entries if n_entries else 0
    self_status = (
        "OK" if self_ratio <= 0.20
        else ("HIGH" if self_ratio <= 0.30 else "EXCESSIVE")
    )
    self_citation_axis = {
        "n_self_citations": n_self,
        "n_total": n_entries,
        "ratio": round(self_ratio, 3),
        "self_authors_input": sorted(self_authors),
        "status": self_status,
        "warning": self_ratio > 0.30,
    }

    # ---- score ----
    n_warnings = sum([
        recency_axis["warning"],
        concentration_axis["warning"],
        geographic_axis["warning"],
        self_citation_axis["warning"],
    ])
    score = round((4 - n_warnings) / 4 * 10, 1)

    # ---- comments ----
    comments: list[str] = []

    comments.append(
        f"References parsed: {n_entries} 件。warnings: {n_warnings}/4 軸。"
    )

    if recency_axis["warning"]:
        comments.append(
            f"⚠ Recency: 中央年 {median_year} (median age {median_age} 年)。"
            "10 年超は『最新研究のフォロー不足』と評価される。"
            "直近 5 年の論文を 1-2 本追加推奨。"
        )

    if concentration_axis["warning"]:
        top1, top1_count = most_common_5[0]
        comments.append(
            f"⚠ Author concentration: 第1著者 '{top1}' が "
            f"{top1_count}/{n_entries} 件 ({top_share*100:.0f}%) を占有。"
            "Echo chamber 印象。多様な研究者からの引用を増やす。"
        )

    if geographic_axis["warning"]:
        comments.append(
            f"⚠ Geographic diversity: regions={dict(region_counter)}。"
            "1 地域の引用に偏っている可能性。国際 review 視点で多様化推奨。"
        )

    if self_citation_axis["warning"]:
        comments.append(
            f"⚠ Self-citation: {n_self}/{n_entries} 件 "
            f"({self_ratio*100:.0f}%) が自己引用。30% 超は reviewer #2 "
            "トリガー (引用操作の疑い)。"
        )

    if n_warnings == 0:
        comments.append(
            f"4 軸すべて健全。citation health score {score}/10。"
        )

    if not self_authors:
        comments.append(
            "ℹ author_last_names が指定されていないので self-citation は "
            "0 計上。本人の姓をカンマ区切りで渡すと精度向上。"
        )

    comments = comments[:7]

    hint = (
        "BibTeX 形式の bib を渡す前提。reference list の text 抜粋では "
        "不正確になることがある。region 判定は journal 名 / address / "
        "author 名から rough に推定 — Unknown が多い場合は信頼性低。"
        "self-citation は author_last_names の姓 lowercase 比較なので、"
        "common な姓 (Tanaka, Smith) は false positive リスクあり。"
    )

    return {
        "score": score,
        "score_max": 10,
        "n_entries": n_entries,
        "n_warnings": n_warnings,
        "per_axis": {
            "recency": recency_axis,
            "author_concentration": concentration_axis,
            "geographic_diversity": geographic_axis,
            "self_citation": self_citation_axis,
        },
        "comments": comments,
        "hint": hint,
        "source": (
            "Bornmann & Daniel (2008) 'What do citation counts measure?' "
            "Journal of Documentation 64(1):45-80 + "
            "ICMJE Recommendations (citation manipulation guidance) + "
            "DORA (San Francisco Declaration on Research Assessment 2012)。"
            "(paper-writing 2026-04 採用)。"
        ),
    }
