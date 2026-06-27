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
    from radia_mcp.grant_writing import tools as t
    if hasattr(t, "grant_writing_health_report"):
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                text = f.read()
            return t.grant_writing_health_report(text, program="generic")
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
