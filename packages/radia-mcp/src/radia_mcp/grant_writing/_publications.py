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

from ..bibliography._bibparse import read_bib_file
from ..bibliography.plans.T14_canonical import CANONICAL

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
        if re.search(me, a, re.I):
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


def grant_writing_publication_list(author: str = "Sugahara|菅原",
                                   since: str = "",
                                   mark: str = "**{}**",
                                   bib_path: str = "") -> str:
    """List the applicant's own publications for a grant achievement section.

    author : regex matched against each author name; the default catches both
             the roman and the Japanese spelling.
    since  : four-digit year; omit for everything.
    mark   : how to highlight the applicant's own name, ``{}`` being the name.
             Grant forms usually want it underlined; ``**{}**`` suits Markdown.
    """
    path = bib_path or str(CANONICAL)
    entries = [e for e in read_bib_file(path) if not e.kind.startswith("@")]
    mine = [e for e in entries
            if re.search(author, e.fields.get("author", "")
                         + " " + e.fields.get("editor", ""), re.I)]
    if since:
        mine = [e for e in mine
                if not e.fields.get("year", "").strip()[:4].isdigit()
                or e.fields["year"].strip()[:4] >= since]

    out = [f"grant_writing_publication_list: author /{author}/"
           + (f", {since} 年以降" if since else "")]
    out.append(f"  {len(mine)} 件が該当（正典 {len(entries)} 件中）")
    out.append("  査読の有無は本ツールでは判定しない。申請者が付すこと。")

    seen = set()
    for kind, label in GROUPS:
        group = [e for e in mine if e.kind.lower() == kind and id(e) not in seen]
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
                    and kind in {"article", "inproceedings"}]
            if gaps:
                out.append(f"   [要確認] 未記入: {', '.join(gaps)}")

    dropped = [e for e in mine if id(e) not in seen]
    if dropped:
        out.append(f"\n## 業績欄に載せなかったもの（{len(dropped)} 件）")
        out.append("   ソフトウェア・未公表・種別不明。必要なら別欄に記載する。")
        for e in dropped:
            out.append(f"   - [{e.kind}] {_clean(e.fields.get('title',''))[:64]}"
                       f"  ({e.fields.get('year','年不明')})")
    return "\n".join(out)
