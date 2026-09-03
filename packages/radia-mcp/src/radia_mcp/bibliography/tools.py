"""Tool functions for the bibliography sub-skill.

Each function is a plain callable; ``__init__.py:register`` wraps them
with ``@mcp.tool()`` at registration time.
"""
from __future__ import annotations

# Tier 1 — online lookup
from .plans.T1_doi_to_bibtex import bibliography_doi_to_bibtex  # noqa: F401
from .plans.T2_arxiv_to_bibtex import bibliography_arxiv_to_bibtex  # noqa: F401
from .plans.T3_search_crossref import bibliography_search_crossref  # noqa: F401

# Tier 2 — parse + edit
from .plans.T4_parse import bibliography_parse  # noqa: F401
from .plans.T5_dedupe import bibliography_dedupe  # noqa: F401
from .plans.T6_canonicalize_keys import bibliography_canonicalize_keys  # noqa: F401

# Tier 3 — lint
from .plans.T7_lint import bibliography_lint  # noqa: F401
from .plans.T8_cite_validation import bibliography_cite_validation  # noqa: F401
from .plans.T9_self_citation_ratio import bibliography_self_citation_ratio  # noqa: F401
from .plans.T12_year_distribution import bibliography_year_distribution  # noqa: F401
from .plans.T13_check_surname_braces import bibliography_check_surname_braces  # noqa: F401
from .plans.T14_canonical import (  # noqa: F401
    bibliography_canonical_path,
    bibliography_make_bbl,
)

# Tier 4 — format conversion
from .plans.T10_normalize_journals import bibliography_normalize_journal_names  # noqa: F401

# Tier 5 — Intelligence Layer
from .plans.T11_health_report import bibliography_health_report  # noqa: F401
