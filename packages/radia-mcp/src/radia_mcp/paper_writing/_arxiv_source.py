"""arXiv LaTeX-source fetch + EM-category search for radia_mcp.paper_writing.

Added 2026-05-26 (radia-mcp v0.89.0).  Absorbed from
``takashiishida/arxiv-latex-mcp`` GitHub survey 2026-05-26.

Why this is useful for EM papers (lab perspective):

  * IEEE TAP / TMag / TMTT papers are math-heavy; reviewers care
    about exact equation forms (e.g. dx/ds = B vs dx/dt = v;
    grad x B vs nabla x B; complex conjugate convention).  When
    citing a related paper, fetching the SOURCE LaTeX lets you
    quote / interpret the math exactly without OCR errors from
    the PDF.

  * arXiv's e-print URL exposes the tarball of the original
    submission's TeX sources.  This module pulls just the main
    .tex file and decoded into a single string for inline review.

  * Companion: the existing ``paper_writing_ieee_download_pdf``,
    ``paper_writing_sciencedirect_download_pdf``,
    ``paper_writing_emerald_download_pdf`` tools handle publisher
    PDFs.  This module is the COMPLEMENTARY arXiv source-LaTeX
    pathway -- when the paper has an arXiv preprint AND the
    reviewer/author wants to read the math (not the typeset PDF).

Sources absorbed:
  * https://github.com/takashiishida/arxiv-latex-mcp
    -- "arxiv-to-prompt under the hood for math interpretation"
  * https://github.com/openags/paper-search-mcp
    -- multi-source arXiv search pattern
  * https://github.com/blazickjp/arxiv-mcp-server
    -- arXiv ID -> Semantic Scholar references/citations link
"""

from __future__ import annotations

import io
import re
import tarfile
import urllib.parse
from typing import Optional


# ============================================================
# Constants
# ============================================================

# EM-relevant arXiv categories (lab default).  Override with categories=""
# to search ALL arXiv (default), or pass a different list.
EM_DEFAULT_CATEGORIES = (
    "eess.SP",           # Signal Processing (often EM-adjacent)
    "physics.class-ph",  # Classical Physics
    "physics.app-ph",    # Applied Physics
    "physics.comp-ph",   # Computational Physics
    "physics.plasm-ph",  # Plasma Physics
    "physics.space-ph",  # Space Physics
    "cs.CE",             # Computational Engineering
    "math.NA",           # Numerical Analysis
)


# ============================================================
# Helpers
# ============================================================


def _require_requests():
    """Lazy import of requests with clear ImportError."""
    try:
        import requests  # type: ignore
        return requests
    except ImportError as e:
        raise ImportError(
            "paper_writing arXiv tools require requests: pip install requests"
        ) from e


def _normalize_arxiv_id(s: str) -> str:
    """Normalize various arXiv ID forms to canonical 'YYMM.NNNNN' or 'arch/yymmNNN'.

    Accepts:
      "2603.17339"
      "arXiv:2603.17339"
      "https://arxiv.org/abs/2603.17339"
      "https://arxiv.org/pdf/2603.17339"
      "https://arxiv.org/abs/2603.17339v2"   (version stripped)
      "physics/0501123"                       (old-style)
    """
    s = s.strip()
    # Strip URL prefix
    for prefix in ("https://arxiv.org/abs/", "http://arxiv.org/abs/",
                    "https://arxiv.org/pdf/", "http://arxiv.org/pdf/"):
        if s.startswith(prefix):
            s = s[len(prefix):]
            break
    # Strip "arXiv:" prefix
    if s.lower().startswith("arxiv:"):
        s = s[6:]
    # Strip ".pdf" suffix
    if s.endswith(".pdf"):
        s = s[:-4]
    # Strip version suffix (vN)
    s = re.sub(r"v\d+$", "", s)
    return s


# ============================================================
# paper_writing_arxiv_fetch_latex_source
# ============================================================


def paper_writing_arxiv_fetch_latex_source(
    arxiv_id: str,
    max_files: int = 5,
    max_chars_per_file: int = 200_000,
) -> dict:
    """Fetch the LaTeX source of an arXiv preprint.

    Downloads the e-print tarball from arxiv.org/e-print/{id}, extracts
    all .tex files, and returns them as a dict.  Identifies the "main"
    .tex file by looking for ``\\documentclass`` declaration.

    POLICY: This complements the publisher-PDF tools
    (``paper_writing_ieee_download_pdf`` etc.) -- use THIS when you
    want to read the AUTHOR'S original equations directly, without
    publisher typesetting or OCR noise.

    Args:
        arxiv_id: arXiv identifier in any common form.  Examples:
            "2603.17339", "arXiv:2603.17339",
            "https://arxiv.org/abs/2603.17339",
            "physics/0501123" (old-style).
        max_files: cap on number of .tex files returned (a large paper
            with appendices may have 10+; default 5 protects against
            blowing the context).
        max_chars_per_file: cap on chars per .tex file.

    Returns:
        dict with:
          "arxiv_id": canonicalized ID
          "main_tex": name of the main .tex file (with \\documentclass)
          "files": list of {"name": str, "size_bytes": int, "content": str}
          "n_files_total": total .tex files in the tarball (may exceed
                            max_files)
          "advice": text hint for next steps

    Sources:
      * https://github.com/takashiishida/arxiv-latex-mcp
      * arxiv-to-prompt: https://github.com/em-bzh/arxiv-to-prompt
    """
    requests = _require_requests()
    aid = _normalize_arxiv_id(arxiv_id)
    url = f"https://arxiv.org/e-print/{aid}"
    try:
        r = requests.get(
            url,
            headers={"User-Agent": "radia-mcp (mailto:ksugahar@ele.kindai.ac.jp)"},
            timeout=60,
        )
        r.raise_for_status()
    except Exception as e:  # noqa: BLE001 - re-raise as clean MCP-JSON
        return {"error": f"arXiv e-print fetch failed for {aid}: {e}"}

    # arXiv returns gzipped tar (occasionally a single .gz of a .tex)
    blob = r.content
    if len(blob) > 50 * 1024 * 1024:
        return {"error": "arXiv source archive exceeds the 50 MiB safety limit"}
    files = []
    main_tex = None
    n_total = 0

    try:
        tar = tarfile.open(fileobj=io.BytesIO(blob), mode="r:*")
    except tarfile.ReadError:
        # Maybe it's a single gzipped .tex (rare but happens)
        import gzip
        try:
            txt = gzip.decompress(blob).decode("utf-8", errors="replace")
            files.append({
                "name": f"{aid}.tex",
                "size_bytes": len(txt.encode("utf-8")),
                "content": txt[:max_chars_per_file],
            })
            return {
                "arxiv_id": aid,
                "main_tex": f"{aid}.tex",
                "files": files,
                "n_files_total": 1,
                "advice": "Single-file .gz arXiv submission (rare path).",
            }
        except Exception as e:  # noqa: BLE001
            return {
                "error": f"arXiv tarball decode failed for {aid}: not tar, "
                         f"not gz, exception={e}",
            }

    for member in tar.getmembers():
        if not member.isfile():
            continue
        if not member.name.endswith(".tex"):
            continue
        n_total += 1
        if len(files) >= max_files:
            continue
        if member.size > max_chars_per_file * 8:
            return {
                "error": (
                    f"arXiv member {member.name!r} is unexpectedly large "
                    f"({member.size} bytes)"
                )
            }
        fh = tar.extractfile(member)
        if fh is None:
            continue
        raw = fh.read(max_chars_per_file * 8 + 1)
        text = raw.decode("utf-8", errors="replace")
        is_main = "\\documentclass" in text and main_tex is None
        if is_main:
            main_tex = member.name
        files.append({
            "name": member.name,
            "size_bytes": member.size,
            "content": text[:max_chars_per_file],
        })

    tar.close()

    if not files:
        return {"error": f"No .tex files found in arXiv tarball for {aid}"}

    if main_tex is None:
        main_tex = files[0]["name"]  # fallback

    return {
        "arxiv_id": aid,
        "main_tex": main_tex,
        "files": files,
        "n_files_total": n_total,
        "advice": (
            f"Main .tex = {main_tex}.  For equation-level inspection use "
            f"paper_writing_arxiv_extract_equations(content).  For citation "
            f"graph use paper_writing_semantic_scholar_lookup('arXiv:{aid}')."
        ),
    }


# ============================================================
# paper_writing_arxiv_extract_equations
# ============================================================


_EQ_PATTERNS = [
    # \begin{equation*?} ... \end{equation*?}
    re.compile(r"\\begin\{equation\*?\}([\s\S]*?)\\end\{equation\*?\}"),
    # \begin{align*?} ... \end{align*?}
    re.compile(r"\\begin\{align\*?\}([\s\S]*?)\\end\{align\*?\}"),
    # \begin{eqnarray*?} ... \end{eqnarray*?}
    re.compile(r"\\begin\{eqnarray\*?\}([\s\S]*?)\\end\{eqnarray\*?\}"),
    # \begin{gather*?} ... \end{gather*?}
    re.compile(r"\\begin\{gather\*?\}([\s\S]*?)\\end\{gather\*?\}"),
    # $$...$$ displayed
    re.compile(r"\$\$([\s\S]+?)\$\$"),
    # \[...\] displayed
    re.compile(r"\\\[([\s\S]+?)\\\]"),
]


def paper_writing_arxiv_extract_equations(
    latex_source: str,
    max_equations: int = 100,
    snippet_chars: int = 600,
) -> dict:
    """Extract all displayed equations from a LaTeX source.

    Useful when reviewing related-work equations or quoting prior-art
    derivations.  Catches the 6 common displayed-math environments:
    equation / align / eqnarray / gather / $$..$$ / \\[...\\].
    Skips inline $..$ math (too noisy for equation-by-equation review).

    Args:
        latex_source: the .tex content (typically from
            paper_writing_arxiv_fetch_latex_source returns).
        max_equations: cap the returned list size.  Default 100.
        snippet_chars: cap per-equation character length.  Default 600
            (long enough for most theorem statements; truncate longer).

    Returns:
        dict with:
          "n_equations": total count (may exceed max_equations)
          "equations": list of {"env": str, "body": str, "char_offset": int}
            ordered by appearance in source.

    Companion to paper_writing_arxiv_fetch_latex_source.
    """
    if not isinstance(latex_source, str):
        return {"error": "latex_source must be a string"}

    eqs = []
    env_names = ["equation", "align", "eqnarray", "gather",
                  "displaymath_dollar", "displaymath_bracket"]
    for env_name, pat in zip(env_names, _EQ_PATTERNS):
        for m in pat.finditer(latex_source):
            body = m.group(1).strip()
            if len(body) > snippet_chars:
                body = body[:snippet_chars] + " ... [truncated]"
            eqs.append({
                "env": env_name,
                "body": body,
                "char_offset": m.start(),
            })

    eqs.sort(key=lambda e: e["char_offset"])
    n_total = len(eqs)
    return {
        "n_equations": n_total,
        "equations": eqs[:max_equations],
        "truncated": n_total > max_equations,
    }


# ============================================================
# paper_writing_arxiv_search
# ============================================================


def paper_writing_arxiv_search(
    query: str,
    max_results: int = 10,
    categories: str = "",
    sort: str = "relevance",
) -> dict:
    """Search arXiv via the official Atom XML API.

    EM-tuned: by default filters to the 8 most likely EM-relevant
    arXiv categories.  Pass an empty ``categories`` string to search
    ALL of arXiv.

    Args:
        query: free-text search query (title / abstract / author).
            Examples:
              "topology optimization magnetostatic"
              "sommerfeld layered green function"
              "method of moments rwg low frequency".
        max_results: cap at 100 (arXiv API limit per request).
            Default 10.
        categories: comma-separated arXiv categories.  Empty string
            (default) means: use EM_DEFAULT_CATEGORIES.  Pass "all"
            to remove the category filter entirely.
            Examples:
              ""           -> default lab EM filter
              "all"        -> no filter (search ALL arXiv)
              "eess.SP"    -> only signal processing
              "physics.class-ph,physics.comp-ph"
        sort: 'relevance' (default), 'lastUpdatedDate',
            'submittedDate'.

    Returns:
        dict with:
          "query": echoed
          "categories_filter": list of arXiv categories actually used
          "n_results": count returned
          "papers": list of {"arxiv_id", "title", "authors", "summary",
                              "published", "categories", "pdf_url",
                              "abs_url"}

    Sources:
      * https://github.com/openags/paper-search-mcp
        (multi-source arXiv search pattern -- here trimmed to just
        arXiv because lab's other publisher tools already handle
        IEEE / ScienceDirect / Emerald).
    """
    requests = _require_requests()
    max_results = max(1, min(int(max_results), 100))

    # Build category filter
    if categories.strip().lower() == "all":
        cat_list: list[str] = []
    elif categories.strip() == "":
        cat_list = list(EM_DEFAULT_CATEGORIES)
    else:
        cat_list = [c.strip() for c in categories.split(",") if c.strip()]

    # arXiv API uses "all:" prefix for full-text search; category
    # filter uses cat:CATNAME.  Combine with AND.
    q_parts = [f"all:{query}"]
    if cat_list:
        cat_clause = " OR ".join(f"cat:{c}" for c in cat_list)
        q_parts.append(f"({cat_clause})")
    search_query = " AND ".join(q_parts)

    sort_map = {
        "relevance": "relevance",
        "lastUpdatedDate": "lastUpdatedDate",
        "submittedDate": "submittedDate",
    }
    sort_by = sort_map.get(sort, "relevance")

    url = "https://export.arxiv.org/api/query"
    params = {
        "search_query": search_query,
        "start": 0,
        "max_results": max_results,
        "sortBy": sort_by,
        "sortOrder": "descending",
    }
    try:
        r = requests.get(
            url,
            params=params,
            headers={"User-Agent": "radia-mcp (mailto:ksugahar@ele.kindai.ac.jp)"},
            timeout=30,
        )
        r.raise_for_status()
    except Exception as e:  # noqa: BLE001
        return {"error": f"arXiv search failed: {e}"}

    # Parse Atom XML
    import xml.etree.ElementTree as ET
    ns = {
        "atom": "http://www.w3.org/2005/Atom",
        "arxiv": "http://arxiv.org/schemas/atom",
    }
    try:
        root = ET.fromstring(r.text)
    except ET.ParseError as e:
        return {"error": f"arXiv search XML parse error: {e}"}

    papers = []
    for entry in root.findall("atom:entry", ns):
        title = (entry.findtext("atom:title", "", ns) or "").strip()
        summary = (entry.findtext("atom:summary", "", ns) or "").strip()
        published = entry.findtext("atom:published", "", ns) or ""
        authors = [
            (a.findtext("atom:name", "", ns) or "").strip()
            for a in entry.findall("atom:author", ns)
        ]
        # arXiv id from <id> URL
        id_url = entry.findtext("atom:id", "", ns) or ""
        if not id_url or title.casefold() == "error":
            summary_text = re.sub(r"\s+", " ", summary).strip()
            return {"error": f"arXiv API error: {summary_text or title}"}
        aid = _normalize_arxiv_id(id_url) if id_url else ""
        # PDF link
        pdf_url = ""
        abs_url = ""
        for link in entry.findall("atom:link", ns):
            if link.get("title") == "pdf":
                pdf_url = link.get("href", "")
            elif link.get("rel") == "alternate":
                abs_url = link.get("href", "")
        # Categories
        cats = [
            cat.get("term", "")
            for cat in entry.findall("atom:category", ns)
            if cat.get("term")
        ]
        papers.append({
            "arxiv_id": aid,
            "title": title,
            "authors": authors,
            "summary": summary[:1000],
            "published": published,
            "categories": cats,
            "pdf_url": pdf_url,
            "abs_url": abs_url,
        })

    return {
        "query": query,
        "categories_filter": cat_list,
        "n_results": len(papers),
        "papers": papers,
    }


# ============================================================
# paper_writing_semantic_scholar_lookup
# ============================================================


def paper_writing_semantic_scholar_lookup(
    paper_id: str,
    fields: str = "title,abstract,authors,year,venue,referenceCount,citationCount,externalIds",
) -> dict:
    """Look up a paper via Semantic Scholar API.

    Accepts DOI, arXiv ID, Semantic Scholar Corpus ID, or URL.

    Args:
        paper_id: any of:
          DOI:           "10.1109/TAP.2006.880747"
          arXiv:         "arXiv:2603.17339" or "2603.17339"
          S2 Corpus ID:  "CorpusID:12345"
          PubMed ID:     "PMID:12345"
          ACL ID:        "ACL:2020.acl-main.1"
          URL:           "https://www.semanticscholar.org/paper/abc123"
        fields: comma-separated fields per Semantic Scholar API.
            Default covers title / abstract / authors / year / venue /
            referenceCount / citationCount / externalIds.

    Returns:
        dict with the paper metadata (forwarded from S2 API), or
        {"error": ...} on failure.

    Sources:
      * https://github.com/blazickjp/arxiv-mcp-server
        (S2 reference / citation traversal pattern)
      * https://api.semanticscholar.org/graph/v1
    """
    requests = _require_requests()

    # Build S2-style ID
    pid = paper_id.strip()
    if pid.startswith("https://www.semanticscholar.org/paper/"):
        pid = pid.split("/paper/")[1].split("/")[0]
    elif pid.lower().startswith("arxiv:"):
        pid = f"ARXIV:{pid[6:]}"
    elif re.match(r"^(?:doi:|https?://(?:dx\.)?doi\.org/|10\.)", pid,
                  re.IGNORECASE):
        pid = re.sub(r"^doi\s*:\s*", "", pid, flags=re.IGNORECASE)
        pid = re.sub(
            r"^https?://(?:dx\.)?doi\.org/", "", pid, flags=re.IGNORECASE
        )
        pid = f"DOI:{pid}"
    elif "/" in pid and re.match(r"^[a-z\-]+/\d{7}$", pid.lower()):
        # old-style arXiv
        pid = f"ARXIV:{pid}"
    # else: pass through as-is (S2 will resolve)

    url = (
        "https://api.semanticscholar.org/graph/v1/paper/"
        + urllib.parse.quote(pid, safe=":")
    )
    params = {"fields": fields}
    try:
        r = requests.get(url, params=params, timeout=30)
        if r.status_code == 404:
            return {"error": f"paper not found in Semantic Scholar: {paper_id}"}
        r.raise_for_status()
    except Exception as e:  # noqa: BLE001
        return {"error": f"Semantic Scholar lookup failed: {e}"}

    data = r.json()
    return {
        "paper_id_input": paper_id,
        "paper_id_resolved": pid,
        "metadata": data,
        "hint": (
            "For the reference graph: pass the returned 'paperId' to "
            "paper_writing_semantic_scholar_references(paper_id).  For "
            "citing-papers: paper_writing_semantic_scholar_citations(paper_id)."
        ),
    }


# ============================================================
# paper_writing_semantic_scholar_references
# ============================================================


def paper_writing_semantic_scholar_references(
    paper_id: str,
    limit: int = 50,
    fields: str = "title,year,venue,externalIds,authors",
) -> dict:
    """List the references CITED BY a given paper.

    Useful for verifying that a related-work section's citations are
    grounded in the lineage that the paper itself uses.

    Args:
        paper_id: as in paper_writing_semantic_scholar_lookup.
        limit: max references to return (S2 max 1000).  Default 50.
        fields: comma-separated S2 fields for each reference.

    Returns:
        dict with "references": list of paper dicts.
    """
    requests = _require_requests()
    pid = paper_id.strip()
    if pid.lower().startswith("arxiv:"):
        pid = f"ARXIV:{pid[6:]}"
    elif pid.startswith("10."):
        pid = f"DOI:{pid}"

    limit = max(1, min(int(limit), 1000))
    url = f"https://api.semanticscholar.org/graph/v1/paper/{pid}/references"
    params = {"fields": fields, "limit": limit}
    try:
        r = requests.get(url, params=params, timeout=30)
        if r.status_code == 404:
            return {"error": f"paper not found: {paper_id}"}
        r.raise_for_status()
    except Exception as e:  # noqa: BLE001
        return {"error": f"S2 references lookup failed: {e}"}
    data = r.json()
    return {
        "paper_id": pid,
        "n_references": len(data.get("data", [])),
        "references": [d.get("citedPaper", {}) for d in data.get("data", [])],
    }


# ============================================================
# paper_writing_semantic_scholar_citations
# ============================================================


def paper_writing_semantic_scholar_citations(
    paper_id: str,
    limit: int = 50,
    fields: str = "title,year,venue,externalIds,authors",
) -> dict:
    """List the papers CITING a given paper.

    Useful for finding more-recent work that built on a foundational
    paper -- e.g. "given Koh-Yook 2006 Sommerfeld impedance plane,
    who has cited it since?".

    Args:
        paper_id: as in paper_writing_semantic_scholar_lookup.
        limit: max citations to return (S2 max 1000).  Default 50.
        fields: comma-separated S2 fields for each citing paper.

    Returns:
        dict with "citations": list of paper dicts.
    """
    requests = _require_requests()
    pid = paper_id.strip()
    if pid.lower().startswith("arxiv:"):
        pid = f"ARXIV:{pid[6:]}"
    elif pid.startswith("10."):
        pid = f"DOI:{pid}"

    limit = max(1, min(int(limit), 1000))
    url = f"https://api.semanticscholar.org/graph/v1/paper/{pid}/citations"
    params = {"fields": fields, "limit": limit}
    try:
        r = requests.get(url, params=params, timeout=30)
        if r.status_code == 404:
            return {"error": f"paper not found: {paper_id}"}
        r.raise_for_status()
    except Exception as e:  # noqa: BLE001
        return {"error": f"S2 citations lookup failed: {e}"}
    data = r.json()
    return {
        "paper_id": pid,
        "n_citations": len(data.get("data", [])),
        "citations": [d.get("citingPaper", {}) for d in data.get("data", [])],
    }


# ============================================================
# Knowledge string -- external paper sources survey
# ============================================================


_EXTERNAL_SOURCES_KNOWLEDGE = r"""
# External paper sources for radia_mcp.paper_writing
# (GitHub MCP survey 2026-05-26)

The lab already has publisher-PDF download tools
(``paper_writing_ieee_download_pdf`` etc.) and a Crossref-based
``paper_writing_fetch_and_cite`` orchestrator.  The 2026-05-26
GitHub survey identified 7 academic-writing MCP servers; the lab
absorbed the 2 most useful PATTERNS:

  1. arXiv LaTeX-source fetch (Ishida -- arxiv-latex-mcp)
  2. Semantic Scholar reference / citation graph
     (Blazick -- arxiv-mcp-server)

## Decision: which tool when?

| Have | Want | Use |
|---|---|---|
| DOI from IEEE / Elsevier / Emerald | full PDF | ``paper_writing_ieee_download_pdf`` etc. |
| arXiv ID | full author LaTeX source (equations exact) | ``paper_writing_arxiv_fetch_latex_source`` |
| arXiv ID | search arXiv for related work in EM | ``paper_writing_arxiv_search`` |
| arXiv .tex source | extract all displayed equations | ``paper_writing_arxiv_extract_equations`` |
| DOI or arXiv ID | paper metadata + abstract + counts | ``paper_writing_semantic_scholar_lookup`` |
| Paper S2 ID | list of REFERENCES (what it cites) | ``paper_writing_semantic_scholar_references`` |
| Paper S2 ID | list of CITATIONS (who cites it) | ``paper_writing_semantic_scholar_citations`` |
| DOI -> BibTeX | citation entry | ``paper_writing_doi_to_bibtex`` (existing) |
| Reference list verification | bibliographic lint | ``paper_writing_lint_reference_format`` |

## EM-tuned defaults

``paper_writing_arxiv_search`` filters by default to 8 categories
likely to contain EM / signal-processing / computational-physics
work:

  eess.SP, physics.class-ph, physics.app-ph, physics.comp-ph,
  physics.plasm-ph, physics.space-ph, cs.CE, math.NA

Pass ``categories="all"`` to remove the filter and search ALL of
arXiv (rarely needed for EM work, but useful for cross-discipline
surveys like "topology optimization in structural mech AND magnetics").

## Surveyed but NOT absorbed (over-scope for lab)

  * paper-search-mcp (openags): 11-source unified search.  The
    lab already has dedicated tools for the 3 publishers it
    actually uses; the other 8 (PubMed, bioRxiv, DBLP, OpenAIRE,
    Springer, Scholar, etc.) are out of scope.
  * Academic-MCP-Server (nanyang12138): 6 databases + citation
    analysis.  Replicates Semantic Scholar functionality.
  * PaperDebugger: multi-agent Research->Critique->Revision.  The
    lab already has ``paper_writing_run_full_workflow`` and
    ``paper_writing_root_cause_diagnosis`` which serve the same
    purpose.
  * citecheck (arxiv:2603.17339): workspace-aware bibliographic
    repair.  Partially overlaps with the lab's existing
    ``paper_writing_lint_reference_format`` +
    ``paper_writing_fetch_and_cite``.

## API key notes

  * arXiv API: no key needed.  Polite rate limit 1 request / 3 seconds.
  * Semantic Scholar API: no key needed for the endpoints used here.
    Rate-limited to 100 requests / 5 minutes for unauthenticated.
    For heavier use, register at
    https://www.semanticscholar.org/product/api for a free API key
    and set ``S2_API_KEY`` env var (NOT YET WIRED -- future
    enhancement).

## Citations / credit

Pattern absorbed FROM:
  * I. Takashi Ishida, ``arxiv-latex-mcp``,
    https://github.com/takashiishida/arxiv-latex-mcp
  * J. Blazick, ``arxiv-mcp-server``,
    https://github.com/blazickjp/arxiv-mcp-server
  * openags, ``paper-search-mcp``,
    https://github.com/openags/paper-search-mcp

(All MIT-licensed; this absorption is by API IDEA only, no source
code copy.)
"""


def paper_writing_external_sources_recipe() -> str:
    """Return the GitHub-survey writeup + decision tree + credits."""
    return _EXTERNAL_SOURCES_KNOWLEDGE
