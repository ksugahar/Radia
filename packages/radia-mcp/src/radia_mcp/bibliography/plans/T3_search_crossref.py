"""Tier 1 — bibliography_search_crossref.

Free-text search the Crossref REST API and return ranked candidates.
Useful when the DOI isn't known but the title/author is.
"""
from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
from .._doi import normalize_doi


_USER_AGENT = ("mcp-server-document/3.0 (mailto:ksugahar@ele.kindai.ac.jp) "
               "Crossref-Plus")


def bibliography_search_crossref(query: str, limit: int = 5) -> str:
    """Search Crossref for the top ``limit`` matches.

    Returns a human-readable list (title / authors / year / DOI) so the
    user can pick the right DOI and call ``bibliography_doi_to_bibtex``.

    Args:
        query: Free-text search (title, author, etc.).
        limit: Max number of results to return (default 5).
    """
    if not isinstance(query, str) or not query.strip():
        return "Error: query must be a nonempty string"
    if type(limit) is not int:
        return "Error: limit must be an integer"
    if limit < 1 or limit > 25:
        limit = 5
    qs = urllib.parse.urlencode({
        "query": query,
        "rows": limit,
        "select": "DOI,title,author,issued,container-title,type",
    })
    url = f"https://api.crossref.org/works?{qs}"
    req = urllib.request.Request(url, headers={
        "User-Agent": _USER_AGENT,
        "Accept": "application/json",
    })
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, urllib.error.HTTPError,
            TimeoutError, OSError, ValueError) as e:
        return f"Error: Crossref query failed: {type(e).__name__}"

    message = data.get("message") if isinstance(data, dict) else None
    items = message.get("items") if isinstance(message, dict) else None
    if (not isinstance(items, list) or data.get("error")
            or data.get("status") not in (None, "ok")
            or not all(isinstance(item, dict) for item in items)):
        return "Error: Crossref returned an invalid search response"
    if not items:
        return f"No results for {query!r}."

    def text(value):
        return " ".join(value.split()) if isinstance(value, str) else ""

    def first_text(value):
        return text(value[0] if value else "") if isinstance(value, list) else text(value)

    lines = [f"bibliography_search_crossref: {query!r} ({len(items)} hit(s))"]
    for i, item in enumerate(items, 1):
        title = first_text(item.get("title"))
        raw_doi = item.get("DOI")
        doi = normalize_doi(raw_doi) if isinstance(raw_doi, str) else ""
        if not title or not re.fullmatch(r"10\.\d{4,9}/\S+", doi):
            return f"Error: Crossref candidate {i} lacks a valid title or DOI"
        title = title[:80] + ("…" if len(title) > 80 else "")
        authors = item.get("author", [])
        first_au = ""
        if isinstance(authors, list) and authors and isinstance(authors[0], dict):
            au = authors[0]
            first_au = text(au.get("family")) or text(au.get("name"))
        year = ""
        issued = item.get("issued", {})
        dp = issued.get("date-parts") if isinstance(issued, dict) else None
        if (isinstance(dp, list) and dp and isinstance(dp[0], list) and dp[0]
                and type(dp[0][0]) is int and 1 <= dp[0][0] <= 9999):
            year = str(dp[0][0])
        container = first_text(item.get("container-title"))
        lines.append(f"  {i}. [{year or '?'}] {first_au or '?'} — {title}")
        if container:
            lines.append(f"      {container[:80]}")
        lines.append(f"      DOI: {doi}")
    return "\n".join(lines)
