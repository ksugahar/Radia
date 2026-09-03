"""Build the applicant's own publication list from the canonical bibliography.

A grant form asks for past achievements, and the list is nearly always the
applicant's own papers. Retyping it per application is how a wrong volume number
gets into three proposals at once, so it comes out of the one bibliography the
lab maintains.

Two things this deliberately does not do:

  It does not claim peer review. Nothing in a BibTeX entry says whether a paper
  was refereed, so entries are grouped by kind and the applicant marks the
  review status. Asserting it would be inventing evidence on a funding
  application.

  It does not quietly include software repositories, in-preparation papers, or
  entries missing a year. They are reported in a separate section, because a
  reviewer who finds a GitHub URL listed among journal papers stops trusting
  the rest of the list.
"""
from __future__ import annotations

import re
from datetime import UTC, datetime
from pathlib import Path

from ..bibliography._bibparse import read_bib_file

# how each BibTeX kind is presented on a Japanese grant form
GROUPS = [
    ("article", "学術論文"),
    ("inproceedings", "国際会議・研究発表"),
    ("incollection", "国際会議・研究発表"),
    ("conference", "国際会議・研究発表"),
    ("techreport", "研究会資料・技術報告"),
    ("book", "著書"),
    ("inbook", "著書"),
]
EXCLUDED_KINDS = {"misc", "unpublished", "online", "software"}


def _resolve_bib_path(bib_path: str) -> Path:
    """Return an explicit bibliography path; grant-writing bundles no private data."""
    if not bib_path:
        raise ValueError(
            "bib_path is required; supply the canonical BibTeX file used for this "
            "application"
        )
    path = Path(bib_path)
    if not path.is_file():
        raise FileNotFoundError(f"bibliography not found: {path}")
    return path


def _clean(s: str) -> str:
    """Strip the brace protection that belongs in .bib, not in a form."""
    s = re.sub(r"\\[a-zA-Z]+\s*", "", s or "")
    return re.sub(r"\s+", " ", s.replace("{", "").replace("}", "")).strip()


def _authors(field: str, me: str, mark: str) -> str:
    """Render authors in order, marking the applicant's own name."""
    out = []
    for a in re.split(r"\s+and\s+", field or ""):
        a = _clean(a)
        if not a:
            continue
        if "," in a:                       # "Surname, Given" -> "Given Surname"
            last, first = (x.strip() for x in a.split(",", 1))
            a = f"{first} {last}".strip()
        if re.search(me, a, re.IGNORECASE):
            a = mark.replace("{}", a)
        out.append(a)
    return ", ".join(out)


def _venue(e) -> str:
    f = e.fields
    parts = [_clean(f.get("journal") or f.get("booktitle") or "")]
    if f.get("volume"):
        parts.append(f"vol. {_clean(f['volume'])}")
    if f.get("number"):
        parts.append(f"no. {_clean(f['number'])}")
    if f.get("pages"):
        parts.append(f"pp. {_clean(f['pages']).replace('--', '-')}")
    if f.get("year"):
        parts.append(_clean(f["year"]))
    return ", ".join(p for p in parts if p)


def grant_writing_publication_list(bib_path: str,
                                   author: str = "Sugahara|菅原",
                                   since: str = "",
                                   mark: str = "**{}**") -> str:
    """List the applicant's own publications for a grant achievement section.

    bib_path: canonical BibTeX file used for this application.
    author : regex matched against each author name; the default catches both
             the roman and the Japanese spelling.
    since  : four-digit year; omit for everything.
    mark   : how to highlight the applicant's own name, ``{}`` being the name.
             Grant forms usually want it underlined; ``**{}**`` suits Markdown.
    """
    path = _resolve_bib_path(bib_path)
    entries = [e for e in read_bib_file(path) if not e.kind.startswith("@")]
    mine = [e for e in entries
            if re.search(author, e.fields.get("author", ""), re.IGNORECASE)]
    if since:
        mine = [e for e in mine
                if e.fields.get("year", "").strip()[:4].isdigit()
                and e.fields["year"].strip()[:4] >= since]

    out = [f"grant_writing_publication_list: author /{author}/"
           + (f", {since} 年以降" if since else "")]
    out.append(f"  {len(mine)} 件が該当（正典 {len(entries)} 件中）")
    out.append("  査読の有無は本ツールでは判定しない。申請者が付すこと。")

    seen = set()
    for label in dict.fromkeys(label for _, label in GROUPS):
        kinds = {kind for kind, group_label in GROUPS if group_label == label}
        group = [e for e in mine if e.kind.lower() in kinds and id(e) not in seen]
        if not group:
            continue
        for e in group:
            seen.add(id(e))
        group.sort(key=lambda e: e.fields.get("year", "0"), reverse=True)
        out.append(f"\n## {label}（{len(group)} 件）")
        for i, e in enumerate(group, 1):
            line = (f"{i}. {_authors(e.fields.get('author',''), author, mark)}, "
                    f"「{_clean(e.fields.get('title',''))}」, {_venue(e)}.")
            if e.fields.get("doi"):
                line += f" DOI: {_clean(e.fields['doi'])}"
            out.append(line)
            gaps = [f for f in ("year", "pages", "doi")
                    if not e.fields.get(f)
                    and e.kind.lower() in {"article", "inproceedings"}]
            if gaps:
                out.append(f"   [要確認] 未記入: {', '.join(gaps)}")
            placeholders = [
                field
                for field, value in e.fields.items()
                if re.search(r"(?i)(?:\bTBD\b|\bTODO\b|0xx|x{3,})", value or "")
            ]
            if placeholders:
                out.append(
                    "   [要確認] プレースホルダを含む書誌項目: "
                    + ", ".join(sorted(placeholders))
                )
            year = e.fields.get("year", "")[:4]
            if year.isdigit() and int(year) > datetime.now(UTC).astimezone().year:
                out.append(
                    "   [要確認] 将来年の業績である。採択・受理・発表予定の状態を明記する。"
                )

    dropped = [e for e in mine if id(e) not in seen]
    if dropped:
        out.append(f"\n## 業績欄に載せなかったもの（{len(dropped)} 件）")
        out.append("   ソフトウェア・未公表・種別不明。必要なら別欄に記載する。")
        for e in dropped:
            out.append(f"   - [{e.kind}] {_clean(e.fields.get('title',''))[:64]}"
                       f"  ({e.fields.get('year','年不明')})")
    return "\n".join(out)


# ---------------------------------------------------------------- count check

# What a proposal claims, and which bibliography entry kinds could support it.
# "査読付" is deliberately mapped to the article count and then reported as
# unverifiable: nothing in a BibTeX entry records refereeing.
_CLAIM_PATTERNS = [
    (r"査読\s*(?:付き?|有り?|あり)\s*(?:の)?(?:学術)?論文", "article", "peer_review"),
    (r"(?:原著|学術|学術雑誌|ジャーナル)\s*論文", "article", ""),
    (r"(?:国際)\s*会議\s*(?:論文|発表|予稿)?", "inproceedings", "international"),
    (r"(?:国内)?\s*研究会\s*(?:資料|発表)?", "techreport", ""),
    (r"peer[- ]reviewed\s+(?:journal\s+)?(?:papers?|articles?)", "article", "peer_review"),
    (r"journal\s+(?:papers?|articles?)", "article", ""),
    (r"(?:international\s+)?conference\s+(?:papers?|presentations?)",
     "inproceedings", "international"),
]
# The English patterns take the English number-first form. They are the entries
# written in ASCII, which is what distinguishes them from the Japanese ones.
_EN_PATTERNS = {p for p, _, _ in _CLAIM_PATTERNS if p.isascii()}
# A particle may sit between the noun and its number. 、 is excluded on purpose:
# it separates one claim from the next, and including it made the number from
# "査読付き論文を12件、国際会議発表を8件" attach to the conference claim too.
_PARTICLE = r"(?:[をはがもでの]|など|about|of)*"
_COUNTER = (r"[\s]*" + _PARTICLE
            + r"[\s]*(\d+)\s*(?:件|編|報|本|篇|papers?|articles?)?")
# number-first: Japanese needs its counter word, English only a space, so a bare
# digit next to a noun is not mistaken for a count in either language
_BEFORE_JA = r"(\d+)\s*(?:件|編|報|本|篇)\s*(?:の)?\s*"
_BEFORE_EN = r"(\d+)\s+"
_KIND_LABEL = {"article": "学術論文", "inproceedings": "国際会議・研究発表",
               "techreport": "研究会資料・技術報告"}


def _scope_year(text: str) -> tuple[str, str]:
    """The period the text says it is counting over, if it says one."""
    m = re.search(r"(\d{4})\s*年\s*以降", text)
    if m:
        return m.group(1), f"{m.group(1)}年以降"
    m = re.search(r"(?:過去|最近|直近)\s*(\d+)\s*年", text)
    if m:
        return "", f"過去{m.group(1)}年"
    return "", "全期間"


def grant_writing_achievement_count_check(text: str,
                                          bib_path: str,
                                          author: str = "Sugahara|菅原") -> str:
    """Check the publication counts a proposal claims against the bibliography.

    A number written by hand into an achievement section drifts as soon as one
    more paper appears, and a reviewer who counts the list and gets a different
    answer has found a reason to doubt everything else on the page.

    A claim about peer review cannot be settled here -- BibTeX does not record
    it -- so those are reported as needing the applicant's confirmation rather
    than passed or failed.

    bib_path is required because the public package does not bundle a private
    laboratory bibliography.
    """
    path = _resolve_bib_path(bib_path)
    entries = [e for e in read_bib_file(path) if not e.kind.startswith("@")]
    mine = [e for e in entries
            if re.search(author, e.fields.get("author", ""), re.IGNORECASE)]

    since, scope_label = _scope_year(text)
    if since:
        mine = [e for e in mine if e.fields.get("year", "")[:4] >= since]
    elif scope_label.startswith("過去"):
        n = int(re.search(r"\d+", scope_label).group())
        cutoff = datetime.now(UTC).astimezone().year - n + 1
        mine = [e for e in mine
                if e.fields.get("year", "")[:4].isdigit()
                and int(e.fields["year"][:4]) >= cutoff]

    counts = {}
    for kind in ("article", "inproceedings", "techreport"):
        counts[kind] = sum(1 for e in mine if e.kind.lower() == kind)

    out = [
        f"grant_writing_achievement_count_check  (集計範囲: {scope_label})",
        (
            f"  書誌中の該当: 学術論文 {counts['article']} / "
            f"会議・研究発表 {counts['inproceedings']} / "
            f"研究会 {counts['techreport']}"
        ),
    ]

    found = 0
    claimed_spans: list[tuple[int, int]] = []
    for pat, kind, uncertain_kind in _CLAIM_PATTERNS:
        before = _BEFORE_EN if pat in _EN_PATTERNS else _BEFORE_JA
        hits = list(re.finditer(pat + _COUNTER, text))
        hits += list(re.finditer(before + pat, text))
        for m in hits:
            if any(m.start() < end and m.end() > start for start, end in claimed_spans):
                continue
            claimed = next((g for g in m.groups() if g and g.isdigit()), "")
            if not claimed:
                continue
            claimed_spans.append(m.span())
            found += 1
            claimed_n, actual = int(claimed), counts[kind]
            phrase = m.group(0).strip().replace("\n", " ")
            if uncertain_kind == "peer_review":
                out.append(f"  [要確認] 「{phrase}」")
                out.append(f"      査読の有無は書誌からは判定できない。"
                           f"{_KIND_LABEL[kind]}の総数は {actual} 件。")
                out.append("      査読付がこの数と一致するか申請者が確認すること。")
            elif uncertain_kind == "international":
                out.append(f"  [要確認] 「{phrase}」")
                out.append(
                    f"      BibTeX種別だけでは国際・国内を判定できない。"
                    f"会議・研究発表の総数は {actual} 件。"
                )
                out.append(
                    "      国際会議に該当する内訳と件数を申請者が確認すること。"
                )
            elif claimed_n == actual:
                out.append(f"  [一致]   「{phrase}」= {actual} 件")
            else:
                out.append(f"  [不一致] 「{phrase}」 と書かれているが、"
                           f"{_KIND_LABEL[kind]}は {actual} 件")
                out.append("      本文か書誌のどちらかが古い。"
                           "grant_writing_publication_list で内訳を確認すること。")

    if not found:
        out.append("  本文に業績件数の記述が見当たらない。"
                   "件数を書くなら書誌と突き合わせること。")
    return "\n".join(out)
