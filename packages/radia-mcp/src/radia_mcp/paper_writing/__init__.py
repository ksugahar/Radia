"""radia_mcp.paper_writing -- journal-paper writing skill suite.

Promoted on 2026-05-26 from
  legacy private source tree
into radia-mcp following the standard ``radia_mcp.<topic>.server``
pattern (statusable, --selftest-able, discoverable via the meta
catalog).

What's inside:

  * tools.py (67 paper_writing_* tools)
       imrad balance, abstract length, hedge counting, figure forward
       reference, equation numbering, citation usage, passive voice
       ratio, paragraph length, sentence variety, reviewer comment
       classification, cover/response letter generation, PDF edge
       overflow check, page-limit validation, ... plus 20 Plan-B
       skill modules in plans/T1-T20 (Tier 1-4 composite scoring,
       reproducibility, statistical reporting, IMRaD discussion
       structure, journal fit, root-cause diagnosis, etc.)
  * cross_lint.py (8 JA-lint wrappers, calling _ja_lint)
  * _ja_lint.py (inline copy of 8 grant_writing JA-lint helpers --
       lab-language linters reused so paper authors get the SAME
       tool names as grant authors)
  * _shared/hedges.py (10 JA hedge patterns shared with the EN
       hedge set inside count_weak_expressions)
  * paper_download.py (IEEE Xplore / ScienceDirect / Emerald PDF
       download with cookie-seeded sessions + fetch-and-cite
       orchestration)
  * _pdf_layout_visual.py (NEW v0.88.0: pymupdf-based image
       rendering for vision-agent layout inspection)
  * _tex_figure_placement.py (NEW v0.88.0: knowledge module on
       LaTeX float specifiers / placeins / widths / anti-patterns /
       per-journal profiles)
  * skill.md (712-line writing-skill knowledge file)
  * plans/T1-T20 (20 Plan-B composite-score helpers)
"""

# Re-export the 67 paper_writing_* tools from tools.py for convenient
# direct import.  Server-side registration happens via server.py.
from . import tools as _tools  # noqa: F401
from ._pdf_layout_visual import (  # noqa: F401
    paper_writing_render_pages_to_png,
    paper_writing_detect_page_whitespace_anomalies,
    paper_writing_layout_thumbnail_strip,
    paper_writing_check_floats_far_from_reference,
    paper_writing_layout_visual_recipe,
)
from ._tex_figure_placement import (  # noqa: F401
    paper_writing_tex_figure_placement,
)

# v0.89.0: external paper-source tools (GitHub survey 2026-05-26):
#   arxiv-latex-mcp + arxiv-mcp-server + paper-search-mcp absorption
from ._arxiv_source import (  # noqa: F401
    paper_writing_arxiv_fetch_latex_source,
    paper_writing_arxiv_extract_equations,
    paper_writing_arxiv_search,
    paper_writing_semantic_scholar_lookup,
    paper_writing_semantic_scholar_references,
    paper_writing_semantic_scholar_citations,
    paper_writing_external_sources_recipe,
)

# v0.90.0: EM-paper-specific style + pre-submission gate orchestrator
from ._em_paper_style import (  # noqa: F401
    paper_writing_em_paper_style,
    paper_writing_em_submission_gate,
)

# v0.91.0: citation verification workflow (references.bib enforcement)
from ._citation_verify import (  # noqa: F401
    paper_writing_verify_citation,
    paper_writing_citation_workflow_recipe,
)

# v0.92.0: PDF overlap/overflow detection + undefined-variable check
from ._pdf_overlap_detection import (  # noqa: F401
    paper_writing_detect_text_image_overlap,
    paper_writing_detect_text_overflow_page,
    paper_writing_detect_overlapping_text_blocks,
    paper_writing_pdf_overlap_recipe,
)
from ._undefined_variables import (  # noqa: F401
    paper_writing_check_undefined_variables,
    UNIVERSAL_WHITELIST,
)

# v0.93.0: multi-file .tex resolver + abstract auto-extract
from ._tex_resolver import (  # noqa: F401
    resolve_input_chain,
    extract_abstract_from_tex,
    paper_writing_resolve_input_chain,
    paper_writing_extract_abstract,
)
