"""radia_mcp.figure — Sugahara Lab publication-figure toolkit.

Promoted on 2026-05-26 from LAB-private figure tooling into radia-mcp as
a standalone subpackage following the `radia_mcp.<topic>.server` pattern.

2026-07-18: figure MERGED into mcp-server-paper-writing (Sugahara).
figure was a shared middle-layer server for paper-writing + presentation;
since presentation merged into paper-writing (2026-07-17), figure follows.
The standalone mcp-server-figure entry point is retired -- all figure_* /
paper_figure_* tools are now served by mcp-server-paper-writing via the
register() function at the end of this module.

The figure MCP tools (registered by register()):
  - figure_style_guide(target='all'|'paper_single_column'|...)
        Lab-standard rules (Times New Roman, units in parentheses,
        no in-figure title, font sizes per embed width).
  - figure_size_for_target(target, embed_width_cm)
        Compute Matplotlib/MATLAB figure size and font sizes that
        will display correctly when embedded at the given column width.
  - figure_tikz_recipe(query, target)
        TikZ / PGFPlots drawing templates and review rules for
        paper-ready schematics, data axes, and MATLAB export.

Python callable helpers (importable; NOT MCP tools):
  - lab_figsize / apply_lab_style / lab_savefig
  - tighten_margins, label_curve_endpoints, add_slope_guide
  - check_legend_overlap, find_best_legend_loc
  - audit_text_overflow
  - plot_asymptote_ratio_sweep, plot_basis_size_convergence

The helpers require matplotlib at runtime, but only when called -- the
MCP server itself loads without matplotlib so `mcp-server-figure
--selftest` can run on machines without a plotting stack.
"""

# `_paper_figure` is the private module that owns the PROFILES + helpers;
# we re-export the public symbols below.  The module is intentionally
# underscore-prefixed so that `from radia_mcp.figure import paper_figure`
# unambiguously gives the FUNCTION (the user's main entry point), not
# the module (which would shadow it under a flat `from .x import x`
# pattern).
#
# 2026-05-26: removed redundant `from . import tools` (line 28) and
# `from . import _paper_figure` (line 35) statements that registered
# the submodules as attributes early.  Under certain .pyc-cache states
# on the self-hosted Windows CI runner, those bare `from . import X`
# lines fired before the submodule was attribute-resolvable in the
# partially-initialized package, producing a spurious "circular import"
# even though the actual import graph is acyclic.  The `from .X import
# symbol` lines below register the submodules as a side effect of
# Python's submodule-resolution logic, so removing the explicit lines
# loses nothing and avoids the race.

# Public helper functions for direct import:
#   from radia_mcp.figure import lab_figsize, apply_lab_style
from .tools import (  # noqa: F401
    lab_figsize,
    apply_lab_style,
    lab_savefig,
    check_min_font,
    minimum_visible_font_pt,
    PAPER_MIN_VISIBLE_FONT_PT,
    PRESENTATION_MIN_VISIBLE_FONT_PT,
    tighten_margins,
    label_curve_endpoints,
    add_slope_guide,
    check_legend_overlap,
    find_best_legend_loc,
    plot_asymptote_ratio_sweep,
    plot_basis_size_convergence,
)

# Paper-quality figure scaffolds (2026-05-26, v0.78.0):
#   from radia_mcp.figure import paper_figure, emit_paper_figure
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
    BEAMER_169_FULL,
    BEAMER_169_HALF,
    paper_figure,
    measure_figure_efficiency,
    auto_tighten,
    add_panel_labels,
    audit_text_overflow,
    emit_paper_figure,
    # v0.80.0 additions (GitHub MCP plotting servers + tueplots + Wong 2011):
    OKABE_ITO,                       # CVD-safe palette
)

# Misuse-proof slide/paper API + audit (2026-06-03): author AT the embed
# width, fail-loud gates on save (no in-figure title, Times New Roman
# actually used, font embedding, no CJK), and emit the exact
# \includegraphics[width=<cm>] snippet.
#   from radia_mcp.figure import lab_figure, save_lab_figure, audit_tex_figures
from ._lab_api import (  # noqa: F401
    lab_figure,
    save_lab_figure,
    legend_no_overlap,
    audit_tex_figures,
    audit_label_overflow,
    _assert_times_new_roman,
)

# The same "author AT the embed width" gate on the SLIDE side: recover the
# authored width of every pasted picture and report the paste scale
# (2026-08-16).
from ._pptx_audit import (  # noqa: F401
    audit_pptx_figures,
    slide_png_from_pdf,
)

# One-call publication-figure builders (return fig+axes via lab_figure):
#   from radia_mcp.figure import scaling_loglog, grouped_bars, convergence, ...
from ._builders import (  # noqa: F401
    scaling_loglog,
    grouped_bars,
    convergence,
    quiver_pair,
    bh_curve,
)


def register(mcp) -> int:
    """Register all figure_* / paper_figure_* tools onto the given FastMCP.

    2026-07-18: figure MERGED into mcp-server-paper-writing (Sugahara):
    figure was a shared middle-layer server for paper-writing + presentation;
    since presentation merged into paper-writing, figure follows. The
    standalone mcp-server-figure entry point is retired.
    """
    from . import tools as _tools
    from . import server as _srv
    count = 0
    for name in sorted(dir(_tools)):
        if name.startswith("figure_") and callable(getattr(_tools, name)):
            mcp.tool()(getattr(_tools, name)); count += 1
    for fn in (_srv.paper_figure_profiles, _srv.paper_figure_recipe,
               _srv.paper_figure_quality_rules, _srv.figure_design_principles,
               _srv.figure_diagram_recipes, _srv.figure_audit_embeds):
        mcp.tool()(fn); count += 1
    return count
