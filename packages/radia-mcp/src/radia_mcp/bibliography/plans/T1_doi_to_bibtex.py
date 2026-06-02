"""Tier 1 — bibliography_doi_to_bibtex.

Look up a DOI via the Crossref REST API and emit a BibTeX entry.

Crossref REST endpoint: https://api.crossref.org/works/{doi}
- No API key required
- Polite usage: include mailto in User-Agent
- Output: structured JSON which we map to BibTeX fields
"""
from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request

from .._bibparse import BibEntry, make_cite_key, write_bib


_CROSSREF_BASE = "https://api.crossref.org/works/"
_USER_AGENT = ("mcp-server-document/3.0 (mailto:ksugahar@ele.kindai.ac.jp) "
               "Crossref-Plus")


def _fetch_crossref(doi: str, timeout: float = 15.0) -> dict | None:
    doi = doi.strip()
    if doi.lower().startswith("https://doi.org/"):
        doi = doi[16:]
    elif doi.lower().startswith("doi:"):
        doi = doi[4:].strip()
    url = _CROSSREF_BASE + urllib.parse.quote(doi, safe="/.")
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT,
                                                "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data.get("message")
    except (urllib.error.URLError, urllib.error.HTTPError,
            TimeoutError, OSError, ValueError):
        return None


_TYPE_MAP = {
    "journal-article": "article",
    "proceedings-article": "inproceedings",
    "book": "book",
    "book-chapter": "incollection",
    "monograph": "book",
    "edited-book": "book",
    "report": "techreport",
    "posted-content": "misc",
    "dissertation": "phdthesis",
}


def _author_field(authors: list[dict]) -> str:
    parts = []
    for a in authors or []:
        family = a.get("family", "").strip()
        given = a.get("given", "").strip()
        if family and given:
            parts.append(f"{family}, {given}")
        elif family:
            parts.append(family)
        elif a.get("name"):
            parts.append(a["name"])
    return " and ".join(parts)


def _year_from_date_parts(date_parts) -> str:
    try:
        return str(date_parts[0][0])
    except (TypeError, IndexError, KeyError):
        return ""


def _crossref_to_bibentry(msg: dict) -> BibEntry:
    kind = _TYPE_MAP.get(msg.get("type", ""), "misc")
    fields: dict[str, str] = {}
    title = msg.get("title", [""])
    if isinstance(title, list):
        title = title[0] if title else ""
    fields["title"] = title or ""
    fields["author"] = _author_field(msg.get("author", []))
    issued = msg.get("issued", {})
    year = _year_from_date_parts(issued.get("date-parts"))
    if year:
        fields["year"] = year
    journal = msg.get("container-title", [""])
    if isinstance(journal, list):
        journal = journal[0] if journal else ""
    if kind == "article" and journal:
        fields["journal"] = journal
    elif kind == "inproceedings" and journal:
        fields["booktitle"] = journal
    elif journal:
        fields["booktitle"] = journal
    vol = msg.get("volume")
    if vol:
        fields["volume"] = str(vol)
    issue = msg.get("issue")
    if issue:
        fields["number"] = str(issue)
    pages = msg.get("page")
    if pages:
        fields["pages"] = pages.replace("-", "--")
    doi = msg.get("DOI") or msg.get("doi")
    if doi:
        fields["doi"] = doi
    pub = msg.get("publisher")
    if pub:
        fields["publisher"] = pub
    fields = {k: v for k, v in fields.items() if v}
    entry = BibEntry(kind=kind, key="", fields=fields)
    entry.key = make_cite_key(entry)
    return entry


def bibliography_doi_to_bibtex(doi: str) -> str:
    """Look up a DOI via Crossref and return a BibTeX entry as text.

    Args:
        doi: The DOI. Accepts bare DOI, ``doi:...`` prefix, or full
             ``https://doi.org/...`` URL.

    Returns:
        BibTeX entry string, or an error message when the lookup fails.
    """
    msg = _fetch_crossref(doi)
    if msg is None:
        return f"Error: Crossref lookup failed for {doi!r}"
    entry = _crossref_to_bibentry(msg)
    return write_bib([entry]).strip()
