"""Tier 1 — bibliography_search_crossref.

Free-text search the Crossref REST API and return ranked candidates.
Useful when the DOI isn't known but the title/author is.
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request


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

    items = data.get("message", {}).get("items", [])
    if not items:
        return f"No results for {query!r}."
    lines = [f"bibliography_search_crossref: {query!r} ({len(items)} hit(s))"]
    for i, item in enumerate(items, 1):
        title = item.get("title", [""])
        if isinstance(title, list):
            title = title[0] if title else ""
        title = title[:80] + ("…" if len(title) > 80 else "")
        authors = item.get("author", [])
        first_au = ""
        if authors:
            au = authors[0]
            first_au = au.get("family") or au.get("name") or ""
        year = ""
        issued = item.get("issued", {})
        dp = issued.get("date-parts") if issued else None
        if dp and isinstance(dp, list) and dp and isinstance(dp[0], list) and dp[0]:
            year = str(dp[0][0])
        container = item.get("container-title", [""])
        if isinstance(container, list):
            container = container[0] if container else ""
        doi = item.get("DOI", "")
        lines.append(f"  {i}. [{year or '?'}] {first_au or '?'} — {title}")
        if container:
            lines.append(f"      {container[:80]}")
        lines.append(f"      DOI: {doi}")
    return "\n".join(lines)
