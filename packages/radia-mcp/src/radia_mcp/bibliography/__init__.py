"""bibliography sub-module — BibTeX lookup, parse, lint, normalize, score.

Promoted from 'paper_writing has 54 tools but assumes .bib exists' gap.
Lab standard: one bundled canonical ``references.bib``; manuscripts cite its
keys and ship a generated ``.bbl`` rather than a private bibliography copy.
Keys follow ``<author><year><word>`` convention (e.g. ``sugahara2026binput``).

Tools (prefix ``bibliography_``):

    Tier 1 — Online lookup:
        bibliography_doi_to_bibtex(doi)
        bibliography_arxiv_to_bibtex(arxiv_id)
        bibliography_search_crossref(query, limit=5)

    Tier 2 — Parse + Edit:
        bibliography_parse(bib_path)
        bibliography_dedupe(bib_path)
        bibliography_canonicalize_keys(bib_path, dry_run=True)

    Tier 3 — Lint:
        bibliography_lint(bib_path)
        bibliography_cite_validation(tex_path, bib_path)
        bibliography_self_citation_ratio(bib_path, author_lastname)
        bibliography_year_distribution(bib_path)
        bibliography_check_surname_braces(bib_path, fix=False)

    Tier 4 — Format conversion:
        bibliography_normalize_journal_names(bib_path, style="ieee", dry_run=True)

    Tier 5 — Intelligence Layer:
        bibliography_health_report(bib_path)

    Tier 6 — Canonical bibliography:
        bibliography_canonical_path()
        bibliography_make_bbl(tex_path, style="", out_path=None)
"""
from . import tools as _tools


def register(mcp) -> int:
    """Register bibliography_* tools with the given FastMCP instance."""
    count = 0
    for name in dir(_tools):
        if name.startswith("bibliography_") and callable(getattr(_tools, name)):
            mcp.tool()(getattr(_tools, name))
            count += 1
    return count
