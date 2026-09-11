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
from .._doi import normalize_doi
from .._metadata_text import bibtex_text


_CROSSREF_BASE = "https://api.crossref.org/works/"
_USER_AGENT = ("mcp-server-document/3.0 (mailto:ksugahar@ele.kindai.ac.jp) "
               "Crossref-Plus")


def _fetch_crossref(doi: str, timeout: float = 15.0) -> dict | None:
    doi = normalize_doi(doi)
    url = _CROSSREF_BASE + urllib.parse.quote(doi, safe="/.")
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT,
                                                "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            message = data.get("message") if isinstance(data, dict) else None
            return message if isinstance(message, dict) else None
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
    if not isinstance(authors, list) or not authors:
        raise ValueError("missing authors")
    parts = []
    for a in authors:
        if not isinstance(a, dict):
            raise ValueError("malformed author")
        if any(a.get(k) is not None and not isinstance(a[k], str)
               for k in ("family", "given", "name")):
            raise ValueError("malformed author name")
        family = bibtex_text(a.get("family") or "")
        given = bibtex_text(a.get("given") or "")
        if family and given:
            parts.append(f"{family}, {given}")
        elif family:
            parts.append(family)
        elif a.get("name"):
            name = bibtex_text(a["name"])
            if not name:
                raise ValueError("empty author name")
            parts.append("{" + name + "}")
        else:
            raise ValueError("incomplete author name")
    return " and ".join(parts)


def _year_from_date_parts(date_parts) -> str:
    try:
        year = date_parts[0][0]
        return str(year) if type(year) is int and 1 <= year <= 9999 else ""
    except (TypeError, IndexError, KeyError):
        return ""


def _crossref_to_bibentry(msg: dict) -> BibEntry:
    if not isinstance(msg, dict):
        raise ValueError("malformed Crossref metadata")
    kind = _TYPE_MAP.get(msg.get("type", ""), "misc")
    fields: dict[str, str] = {}
    title = msg.get("title", [""])
    if isinstance(title, list):
        title = title[0] if title else ""
    if not isinstance(title, str) or not title.strip():
        raise ValueError("missing title")
    fields["title"] = bibtex_text(title)
    if not fields["title"]:
        raise ValueError("empty title")
    fields["author"] = _author_field(msg.get("author", []))
    year = ""
    for key in ("published", "issued"):
        record = msg.get(key)
        if isinstance(record, dict):
            year = _year_from_date_parts(record.get("date-parts"))
            if year:
                break
    if not year:
        raise ValueError("missing publication year")
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
        fields["pages"] = re.sub(r"-+", "--", str(pages))
    doi = msg.get("DOI") or msg.get("doi")
    if doi:
        fields["doi"] = doi
    pub = msg.get("publisher")
    if pub:
        fields["publisher"] = pub
    for key in ("journal", "booktitle", "volume", "number", "pages", "publisher"):
        if key in fields:
            if not isinstance(fields[key], str):
                raise ValueError(f"malformed {key}")
            fields[key] = bibtex_text(fields[key])
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
    doi = normalize_doi(doi)
    if not re.fullmatch(r"10\.\d{4,9}/\S+", doi):
        return "Error: invalid DOI identifier"
    msg = _fetch_crossref(doi)
    if msg is None:
        return f"Error: Crossref lookup failed for {doi!r}"
    try:
        returned = msg.get("DOI") or msg.get("doi")
        if not isinstance(returned, str) or normalize_doi(returned).casefold() != doi.casefold():
            raise ValueError("returned DOI does not match requested DOI")
        entry = _crossref_to_bibentry(msg)
    except (ValueError, TypeError, AttributeError) as exc:
        return f"Error: verify Crossref metadata manually: {exc}"
    return write_bib([entry]).strip()
