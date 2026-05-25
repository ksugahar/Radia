"""radia_mcp.graph — Sugahara Lab publication-figure style guide.

Promoted on 2026-05-26 from
  s:/mcp-server/src/mcp_server_document/graph/
into radia-mcp as a standalone subpackage following the
`radia_mcp.<topic>.server` pattern (statusable, --selftest-able,
discoverable via the meta catalog).

Two MCP tools:
  - graph_style_guide(target='all'|'paper_single_column'|...)
        Lab-standard rules (Times New Roman, units in parentheses,
        no in-figure title, font sizes per embed width).
  - graph_size_for_target(target, embed_width_cm)
        Compute Matplotlib/MATLAB figure size and font sizes that
        will display correctly when embedded at the given column width.

Python callable helpers (importable; NOT MCP tools):
  - lab_figsize / apply_lab_style / lab_savefig
  - tighten_margins, label_curve_endpoints, add_slope_guide
  - check_legend_overlap, find_best_legend_loc
  - plot_asymptote_ratio_sweep, plot_basis_size_convergence

The helpers require matplotlib at runtime, but only when called -- the
MCP server itself loads without matplotlib so `mcp-server-graph
--selftest` can run on machines without a plotting stack.
"""

from . import tools          # noqa: F401
# `_paper_figure` is the private module that owns the PROFILES + helpers;
# we re-export the public symbols below.  The module is intentionally
# underscore-prefixed so that `from radia_mcp.graph import paper_figure`
# unambiguously gives the FUNCTION (the user's main entry point), not
# the module (which would shadow it under a flat `from .x import x`
# pattern).
from . import _paper_figure   # noqa: F401

# Public helper functions for direct import:
#   from radia_mcp.graph import lab_figsize, apply_lab_style
from .tools import (  # noqa: F401
    lab_figsize,
    apply_lab_style,
    lab_savefig,
    tighten_margins,
    label_curve_endpoints,
    add_slope_guide,
    check_legend_overlap,
    find_best_legend_loc,
    plot_asymptote_ratio_sweep,
    plot_basis_size_convergence,
)

# Paper-quality figure scaffolds (2026-05-26, v0.78.0):
#   from radia_mcp.graph import paper_figure, emit_paper_figure
from ._paper_figure import (  # noqa: F401
    PaperProfile,
    PROFILES,
    get_profile,
    IEEE_SINGLE_COLUMN,
    IEEE_DOUBLE_COLUMN,
    IEEJ_SINGLE_COLUMN,
    IEEJ_DOUBLE_COLUMN,
    IGTE_DIGEST_DOUBLE,
    IGTE_DIGEST_SINGLE,
    paper_figure,
    measure_figure_efficiency,
    auto_tighten,
    add_panel_labels,
    emit_paper_figure,
    # v0.80.0 additions (GitHub MCP plotting servers + tueplots + Wong 2011):
    OKABE_ITO,                       # CVD-safe palette
)
