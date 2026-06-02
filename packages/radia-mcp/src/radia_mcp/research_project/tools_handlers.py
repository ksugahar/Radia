"""Map artifact kind -> health-report function.

Lazy-imported so the research_project package doesn't pull every
sub-skill at module load.
"""
from __future__ import annotations


def _paper_health(path: str) -> str:
    from radia_mcp.paper_writing import tools as t
    if hasattr(t, "paper_writing_health_report"):
        return t.paper_writing_health_report(path)
    return f"(no paper_writing_health_report for {path})"


def _poster_health(path: str) -> str:
    from radia_mcp.poster import tools as t
    return t.poster_health_report(path)


def _presentation_pptx_health(path: str) -> str:
    from radia_mcp.doc_convert import tools as t
    return t.doc_convert_health_report(path, kind="presentation")


def _figure_pptx_health(path: str) -> str:
    from radia_mcp.doc_convert import tools as t
    return t.doc_convert_health_report(path, kind="figure", target_width_cm=8.0)


def _bib_health(path: str) -> str:
    from radia_mcp.bibliography import tools as t
    return t.bibliography_health_report(path)


def _pdf_health(path: str) -> str:
    from radia_mcp.pdf import tools as t
    return t.pdf_health_report(path)


def _grant_health(path: str) -> str:
    # grant_writing lint is LAB-private (mcp-server-document) and is NOT part of
    # radia-mcp.  Use it opportunistically if that private package happens to be
    # installed (e.g. on the lab machine); otherwise report it unavailable.
    try:
        from mcp_server_document.grant_writing import tools as t  # LAB-private, optional
    except ImportError:
        return ("(grant health is LAB-private to mcp-server-document; "
                "not available in this radia-mcp install)")
    if hasattr(t, "grant_writing_health_report"):
        # Many grant-writing health tools accept a string of text;
        # we pass the file contents to be safe.
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                text = f.read()
            return t.grant_writing_health_report(text)
        except OSError:
            return f"(could not read {path})"
    return f"(no grant_writing_health_report for {path})"


handlers = {
    "paper":              _paper_health,
    "poster":             _poster_health,
    "presentation_pptx":  _presentation_pptx_health,
    "figure_pptx":        _figure_pptx_health,
    "bib":                _bib_health,
    "output_pdf":         _pdf_health,
    "grant":              _grant_health,
}
