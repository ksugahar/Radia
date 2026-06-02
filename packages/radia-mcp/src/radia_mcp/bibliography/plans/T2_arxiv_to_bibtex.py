"""Tier 1 — bibliography_arxiv_to_bibtex.

Query the arXiv API for a paper ID and emit a BibTeX entry.

arXiv API endpoint: https://export.arxiv.org/api/query?id_list=...
- Atom XML response
- Polite usage: 1 req / 3 sec, include mailto if scripting
"""
from __future__ import annotations

import re
import ssl
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET


def _ssl_context() -> ssl.SSLContext | None:
    """Build an SSL context using certifi's bundle when available.

    Some Python installs lack a complete trust store for hosts like
    export.arxiv.org. certifi (bundled with pip / requests) covers
    that gap. Returns ``None`` to let urllib use its default.
    """
    try:
        import certifi  # type: ignore
    except ImportError:
        return None
    ctx = ssl.create_default_context(cafile=certifi.where())
    return ctx

from .._bibparse import BibEntry, make_cite_key, write_bib


_ARXIV_API = "https://export.arxiv.org/api/query"
_NS = {"a": "http://www.w3.org/2005/Atom",
       "arxiv": "http://arxiv.org/schemas/atom"}
_USER_AGENT = ("mcp-server-document/3.0 "
               "(mailto:ksugahar@ele.kindai.ac.jp)")


def _normalize_arxiv_id(s: str) -> str:
    s = s.strip()
    if s.lower().startswith("arxiv:"):
        s = s[6:]
    if s.lower().startswith("https://arxiv.org/abs/"):
        s = s[22:]
    # Strip version suffix like v2 for the lookup, the API accepts both.
    return s


def _fetch_arxiv(arxiv_id: str, timeout: float = 15.0):
    qs = urllib.parse.urlencode({"id_list": arxiv_id})
    req = urllib.request.Request(f"{_ARXIV_API}?{qs}",
                                  headers={"User-Agent": _USER_AGENT})
    ctx = _ssl_context()
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
            return resp.read().decode("utf-8")
    except (urllib.error.URLError, urllib.error.HTTPError,
            TimeoutError, OSError):
        return None


def _atom_to_bibentry(atom_xml: str) -> BibEntry | None:
    try:
        root = ET.fromstring(atom_xml)
    except ET.ParseError:
        return None
    entry = root.find("a:entry", _NS)
    if entry is None:
        return None
    title = (entry.findtext("a:title", default="", namespaces=_NS) or "").strip()
    summary = (entry.findtext("a:summary", default="", namespaces=_NS) or "").strip()
    published = (entry.findtext("a:published", default="", namespaces=_NS) or "").strip()
    year = published[:4] if published else ""
    authors = []
    for a in entry.findall("a:author", _NS):
        nm = (a.findtext("a:name", default="", namespaces=_NS) or "").strip()
        if nm:
            authors.append(nm)
    # Convert "First Last" → "Last, First" if possible.
    formatted_authors = []
    for nm in authors:
        if "," in nm:
            formatted_authors.append(nm)
        else:
            parts = nm.rsplit(" ", 1)
            if len(parts) == 2:
                formatted_authors.append(f"{parts[1]}, {parts[0]}")
            else:
                formatted_authors.append(nm)
    fields = {
        "title": " ".join(title.split()),
        "author": " and ".join(formatted_authors),
        "year": year,
        "eprint": entry.findtext("a:id", default="", namespaces=_NS).rsplit("/", 1)[-1],
        "archivePrefix": "arXiv",
        "primaryClass": (entry.find("arxiv:primary_category", _NS) or ET.Element("x")).get("term", ""),
        "abstract": " ".join(summary.split()),
    }
    fields = {k: v for k, v in fields.items() if v}
    out = BibEntry(kind="misc", key="", fields=fields)
    out.key = make_cite_key(out)
    return out


def bibliography_arxiv_to_bibtex(arxiv_id: str) -> str:
    """Query arXiv for a paper id and return a BibTeX entry.

    Args:
        arxiv_id: e.g. ``2603.17339`` or ``arxiv:2603.17339`` or
                  ``https://arxiv.org/abs/2603.17339``.

    Returns:
        BibTeX entry text, or an error message.
    """
    aid = _normalize_arxiv_id(arxiv_id)
    body = _fetch_arxiv(aid)
    if body is None:
        return f"Error: arXiv lookup failed for {aid!r}"
    entry = _atom_to_bibentry(body)
    if entry is None or "title" not in entry.fields:
        return f"Error: no entry returned for {aid!r}"
    return write_bib([entry]).strip()
