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

from . import tools  # noqa: F401  -- public re-export for `from radia_mcp.graph import tools`

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
