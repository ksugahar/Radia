"""
Graph MCP Server (radia_mcp.graph)

Surfaces the Sugahara Lab publication-figure style guide via two MCP
tools:

    graph_style_guide(target='all'|'paper_single_column'|...)
        Returns the lab-standard rules (Times New Roman, units in
        parentheses, no in-figure title, font sizes per embed width,
        IEEE/IEEJ conventions) plus the numeric profile for the
        selected target.

    graph_size_for_target(target, embed_width_cm=None)
        Computes the recommended figure size + font sizes + a
        ready-to-paste MATLAB or Matplotlib snippet for the chosen
        embedding width.

Promoted on 2026-05-26 from
  s:/mcp-server/src/mcp_server_document/graph/
into radia-mcp following the standard `radia_mcp.<topic>.server`
pattern (statusable, --selftest-able, discoverable via the meta
catalog).

The companion Python helpers (apply_lab_style, lab_figsize,
lab_savefig, tighten_margins, label_curve_endpoints, add_slope_guide,
check_legend_overlap, find_best_legend_loc, plot_asymptote_ratio_sweep,
plot_basis_size_convergence) are NOT registered as MCP tools — they
have richer signatures (matplotlib objects in/out) than MCP-JSON can
carry.  They live in `radia_mcp.graph.tools` for direct import by
analysis scripts.

Usage:
    mcp-server-graph              # start MCP server (stdio transport)
    mcp-server-graph --selftest   # smoke check every profile in both tools
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool

from . import tools as _tools
from . import _paper_figure as _paper

mcp = FastMCP("mcp-server-graph")


# ============================================================
# Tool registration
# ============================================================
# Auto-register the MCP-callable subset (graph_*).  The Python helpers
# (apply_lab_style etc.) have matplotlib-typed signatures that don't
# round-trip through MCP-JSON, so we skip them deliberately.
_REGISTERED: list[str] = []
for _name in dir(_tools):
    if _name.startswith("graph_"):
        _fn = getattr(_tools, _name)
        if callable(_fn):
            mcp.tool()(_fn)
            _REGISTERED.append(_name)


# ============================================================
# Paper-quality figure MCP tools (2026-05-26, v0.78.0)
# ============================================================
# These do NOT execute matplotlib — they return text recipes / catalog
# data that the AI passes to a separate Python execution.  This keeps
# the MCP server purely informational (no matplotlib import at server
# load time, no figure rendering in MCP RPC, no PIL/PDF dependencies).
# The actual rendering happens when a user runs the recipe locally:
#   from radia_mcp.graph import paper_figure, emit_paper_figure
#   fig, axes = paper_figure('ieee_double_column', nrows=1, ncols=2)
#   ... plot ...
#   emit_paper_figure(fig, 'out', 'ieee_double_column')


@mcp.tool()
def paper_figure_profiles(query: str = "all") -> str:
    """List paper-quality figure profiles + their exact journal geometry.

    Args:
        query: 'all' (default) for every profile, or a profile key
            (e.g. 'ieee_double_column') for the single-profile spec.

    Returns multi-line text with the width in mm, font/legend size,
    default aspect, default subplots_adjust margins, and the upstream
    journal spec URL.

    Use this BEFORE calling `paper_figure_recipe` so you know which
    profile name to pass.
    """
    q = (query or "all").strip().lower()
    if q == "all":
        lines = ["radia_mcp.graph paper-figure profiles:", ""]
        for name, prof in _paper.PROFILES.items():
            lines.append(f"  {name}")
            lines.append(f"    {prof.full_name}")
            lines.append(f"    width = {prof.width_mm:.2f} mm "
                          f"({prof.width_in:.3f} in)")
            lines.append(f"    font  = {prof.font_pt:.1f} pt  "
                          f"legend = {prof.legend_pt:.1f} pt  "
                          f"tick = {prof.tick_pt:.1f} pt")
            lines.append(f"    default_aspect (h/w) = {prof.default_aspect:.2f}")
            lines.append(f"    margins L/R/T/B = "
                          f"{prof.margin_left:.3f}/{prof.margin_right:.3f}/"
                          f"{prof.margin_top:.3f}/{prof.margin_bottom:.3f}")
            lines.append(f"    wspace/hspace = "
                          f"{prof.wspace:.2f}/{prof.hspace:.2f}")
            if prof.spec_url:
                lines.append(f"    spec: {prof.spec_url}")
            lines.append("")
        return "\n".join(lines)
    # Single profile drill-down
    try:
        prof = _paper.get_profile(q)
    except ValueError as e:
        return str(e) + "\nUse paper_figure_profiles('all') to list."
    out = [
        f"Profile: {prof.name}",
        f"  full_name = {prof.full_name}",
        f"  width_mm  = {prof.width_mm:.2f} ({prof.width_in:.3f} in)",
        f"  column    = {prof.column}",
        f"  font_pt   = {prof.font_pt}",
        f"  legend_pt = {prof.legend_pt}",
        f"  tick_pt   = {prof.tick_pt}",
        f"  linewidth_pt      = {prof.linewidth_pt}",
        f"  axes_linewidth_pt = {prof.axes_linewidth_pt}",
        f"  marker_size_pt    = {prof.marker_size_pt}",
        f"  margin_left   = {prof.margin_left}",
        f"  margin_right  = {prof.margin_right}",
        f"  margin_top    = {prof.margin_top}",
        f"  margin_bottom = {prof.margin_bottom}",
        f"  wspace = {prof.wspace}    hspace = {prof.hspace}",
        f"  default_aspect = {prof.default_aspect}",
        f"  spec: {prof.spec_url}",
    ]
    return "\n".join(out)


@mcp.tool()
def paper_figure_recipe(
    profile: str = "ieee_double_column",
    nrows: int = 1,
    ncols: int = 1,
    panel_labels: bool = False,
    aspect: float | None = None,
) -> str:
    """Generate a self-contained Python recipe for a paper-quality figure.

    Output is ready-to-paste code that:
      1. imports paper_figure + emit_paper_figure from radia_mcp.graph
      2. calls paper_figure(profile, nrows, ncols, ...) to create
         (fig, axes_2d) at the journal's EXACT width with the
         pre-tuned subplots_adjust margins for that layout
      3. has a placeholder `# plot here` block per axis
      4. ends with emit_paper_figure(fig, 'out', profile) which acts
         as the GATE: raises ValueError if axes_area / total_area <
         0.78, with a per-margin suggestion of which margin is the
         biggest waste

    Args:
        profile: A key from paper_figure_profiles (e.g.
            'ieee_double_column', 'ieej_single_column',
            'igte_digest_double').
        nrows, ncols: Subplot grid.  Profiles' subplots_adjust deltas
            for common layouts (1x1, 1x2, 1x3, 2x1, 2x2) are pre-baked
            into paper_figure() so you do NOT need to tune them.
        panel_labels: True to auto-place (a), (b), (c)... in each panel
            (skipped for 1x1).
        aspect: Override the figure aspect (h/w).  Default None uses
            the profile's recommended aspect.

    Returns:
        A Python recipe as a string, ready to paste into a script.
    """
    try:
        prof = _paper.get_profile(profile)
    except ValueError as e:
        return f"# ERROR: {e}\n# Use paper_figure_profiles() to list profiles."

    nax = nrows * ncols
    panel_str = ""
    if panel_labels and nax > 1:
        panel_str = ", panel_labels=True"
    aspect_str = ""
    if aspect is not None:
        aspect_str = f", aspect={aspect}"

    plot_loop = []
    if nax == 1:
        plot_loop.append("ax = axes[0, 0]")
        plot_loop.append("# --- plot here ---")
        plot_loop.append("# ax.plot(x, y, label='...')")
        plot_loop.append(r"# ax.set_xlabel(r'$f$ (Hz)')  # units in PARENS,"
                          " not [Hz] (IEEE/IEEJ convention)")
        plot_loop.append(r"# ax.set_ylabel(r'$|Z|$ ($\Omega$)')")
        plot_loop.append("# ax.legend(loc='best', frameon=False)")
    else:
        plot_loop.append("for i, ax in enumerate(axes.flat):")
        plot_loop.append("    # --- plot here ---")
        plot_loop.append("    # ax.plot(x, y, label='...')")
        plot_loop.append(r"    # ax.set_xlabel(r'$f$ (Hz)')")
        plot_loop.append(r"    # ax.set_ylabel(r'$|Z|$ ($\Omega$)')")
        plot_loop.append("    # ax.legend(loc='best', frameon=False)")

    recipe = f"""# Paper-quality figure: {prof.full_name}
# width = {prof.width_mm:.1f} mm exactly ({prof.width_in:.3f} in)
# layout = {nrows} x {ncols} = {nax} panel(s)
#
# Pre-tuned for: font {prof.font_pt} pt, axes_linewidth {prof.axes_linewidth_pt}
#                pt, marker {prof.marker_size_pt} pt, wspace
#                {prof.wspace:.2f}, hspace {prof.hspace:.2f}
# Margins set via per-layout deltas for {nrows}x{ncols} (see paper_figure.py
# _MARGIN_DELTAS).  axes-area / total-area should land near 0.80-0.88.

from radia_mcp.graph import paper_figure, emit_paper_figure

fig, axes = paper_figure(
    profile={prof.name!r},
    nrows={nrows}, ncols={ncols}{aspect_str}{panel_str},
)

{chr(10).join(plot_loop)}

# Save with validation gate.  on_fail='raise' will refuse the save if
# axes_area / total_area < 0.78 and tell you which margin to cut.
# Use on_fail='auto_tighten' to let the gate fix it instead of you.
emit_paper_figure(
    fig,
    path='out',                     # writes out.pdf + out.png at 600 DPI
    profile={prof.name!r},
    min_axes_fraction=0.72,         # paper-quality floor (achievable
                                    # after auto_tighten); raise to
                                    # 0.78-0.82 for ultra-tight cases
    on_fail='raise',                # 'raise' | 'warn' | 'auto_tighten'
)
"""
    return recipe


@mcp.tool()
def paper_figure_quality_rules(query: str = "all") -> str:
    """Why paper-quality figures need a margin-efficiency gate.

    Returns text on what 'axes-area / total-area fraction' is, why the
    paper-quality floor is around 0.78, what aspects to watch (units in
    parentheses not brackets, no in-figure title, TrueType font embed
    pdf.fonttype=42, 8-9 pt font for IEEE/IEEJ figure text), and how
    auto_tighten + emit_paper_figure compose to make the workflow
    refuse to ship a wasteful figure.

    Topics:
        'all'             - full text
        'efficiency'      - how axes_area_fraction is computed + thresholds
        'margins'         - per-margin breakdown reading guide
        'units'           - units in parentheses convention (IEEE/IEEJ)
        'font_embedding'  - Type 42 requirement
        'multipanel'      - 1x2 / 2x1 / 2x2 layout tactics
        'side_by_side'    - two figures in 8 cm -> each <= 4 cm,
                            font 10 pt, legend 8-9 pt, no overlap
    """
    rules = {
        "efficiency": """\
[efficiency]

`axes_area_fraction = sum(axes.bbox_in2) / fig_area_in2`.

WHAT THIS COUNTS AS WASTE (= 余白):
  the white region BETWEEN the axes box outer edges and the FIGURE
  BOUNDING BOX outer edges.  This is the area that carries NO data,
  no labels, no tick numbers -- pure white pixels around the axes.
  The lab principle is "情報がなく無駄はやめる": every mm of the
  figure bbox should either be axes interior, axis label, tick label,
  tick mark, or legend.  Idle whitespace = waste.

WHAT IS *NOT* COUNTED AS WASTE:
  - space INSIDE the axes between curves and the axis frame (that is
    the data window; it carries information about magnitudes)
  - x/y axis labels in the bottom/left margins (they carry the units
    and quantity name)
  - tick labels (carry the numerical scale)

Empirical thresholds (after auto_tighten):
  >= 0.85   excellent
  0.72-0.85 paper-quality, no obvious waste
  0.65-0.72 first-draft acceptable
  < 0.65    reviewer-visible waste -- fix via auto_tighten() or
            tighten subplots_adjust manually
""",
        "font_rule": """\
[font_rule]

SUGAHARA LAB ABSOLUTE FONT RULE:

  Bounding-box width = 8 cm  ->  figure-text font MUST be 10 pt.

The 10 pt size matches IEEE / IEEJ body text.  When the figure prints
at 100% scale alongside body text, the figure text is the same size
as the surrounding paragraphs -- the reader doesn't have to squint.

WIDER COLUMNS (16-18 cm \\figure* / digest full-width):
  The font stays 10 pt.  Wider columns get a BIGGER AXES BOX, NOT
  bigger text.  The text is ABSOLUTE, not relative to figure size.

  Wrong intuition: "the figure is twice as wide, so font should be
                    twice as big (20 pt)" -- this is WRONG.  At 18 cm
                    a 20-pt font reads as a billboard, not a paper.
  Right intuition: a 10-pt font on 18 cm = 10 pt on 8 cm.  Same
                    readability; bigger axes for the data.

LEGEND / TICK FONT:
  legend = 10 pt (same as body, never shrunk -- legends are
                  data, not afterthoughts)
  tick   = 9  pt (one pt below body; IEEE / IEEJ convention)

This rule pins every paper_figure profile's font_pt to 10.0 regardless
of column width.  If you find yourself wanting smaller text "to make
the axes fit", the answer is to use a wider column or simplify the
plot -- never to shrink font below 10 pt.

EXCEPTION -- two panels side-by-side in 8 cm (each ~4 cm wide):
  When you place TWO graphs in an 8 cm width (1 row x 2 cols, each
  sub-panel <= 4 cm), the BODY / AXIS font stays 10 pt, but the
  LEGEND drops to 8-9 pt.  A 10 pt legend physically crowds a 4 cm
  panel; 8-9 pt keeps it readable without stealing data area.  This
  is the ONE sanctioned place to shrink the legend below body size.
  The legend must still NOT overlap the graph (see no_legend_overlap).
  See the `side_by_side` topic and the
  'digest_double_column_side_by_side' profile.
""",
        "no_title_in_figure": """\
[no_title_in_figure]

LAB RULE: NO TITLES INSIDE THE FIGURE.

  WRONG:  ax.set_title('AC impedance vs frequency')
  WRONG:  fig.suptitle('Figure 3')
  RIGHT:  \\caption{AC impedance vs frequency.  Sample frequency
          1-100 kHz, ...}

Reason:
  - Titles in figures are duplicated in the LaTeX caption -> visual
    noise.
  - The title text steals 10-15% of figure height that could go to
    the axes.
  - Authors edit captions in LaTeX; updating figure titles requires
    re-rendering.  Sole source of truth = the caption.

emit_paper_figure() defaults to `check_title_in_figure=True` and
RAISES ValueError when a title is found.  Pass `False` only if you
are intentionally producing a slide-deck figure (not a paper figure).
""",
        "colorblind_safe": """\
[colorblind_safe]

LAB RULE: every line color must be either GREYSCALE or in the
Okabe-Ito 8-color colorblind-safe palette.

Okabe-Ito palette (Wong 2011, Nature Methods 8:441):
  black        #000000
  orange       #E69F00
  sky-blue     #56B4E9
  bluish-green #009E73
  yellow       #F0E442   (use sparingly -- low contrast on white)
  blue         #0072B2
  vermillion   #D55E00
  reddish-purple #CC79A7

Why:
  ~8% of men are red-green colorblind (deuteranopia / protanopia).
  Matplotlib's default `tab10` cycle has confusable red+orange that
  these readers cannot distinguish.  Okabe-Ito is ALSO distinguishable
  in greyscale print, important for low-cost mono printers.

paper_figure() sets the lab default rcParams['axes.prop_cycle'] to
Okabe-Ito.  emit_paper_figure() additionally LINTS every Line2D color
against the palette + greyscale exception, and raises on violation.

Override with `check_colorblind_safe=False` ONLY for project-mandated
brand palettes you've verified externally.
""",
        "font_embedding": """\
[font_embedding]

LAB RULE: PDFs MUST embed TrueType (Type-42) fonts.

Matplotlib defaults: `pdf.fonttype = 3` (raster glyphs).  Type-3 PDFs
blur when zoomed AND fail IEEE / Elsevier pre-flight.  paper_figure()
sets rcParams['pdf.fonttype'] = 42 (TrueType) at scaffold time.

emit_paper_figure() then verifies AFTER save by scanning the PDF
binary for /Subtype /Type3.  Catches the case where a downstream
import or user code reset rcParams between scaffold and save.

Verify manually:
  pdffonts out.pdf | grep -i type
  # All entries should be 'Type 42' (TrueType).
""",
        "no_legend_overlap": """\
[no_legend_overlap]

LAB RULE: LEGENDS MUST NOT OVERLAP DATA LINES.

A legend that overlaps even ONE curve is the single most common
reviewer-visible obvious flaw -- programmatic detection is reliable,
eyeball-check at small embed size is not.

How emit_paper_figure() detects it:
  For each axis-legend pair, sample 200 equally-spaced points along
  every Line2D in the axes and test each against the legend bbox.
  Any point INSIDE the bbox = overlap = fail the gate.

How to fix (in order of lab preference):

  1. Direct labels in the right margin (BEST for time-series, sweeps):
       from radia_mcp.graph import label_curve_endpoints
       label_curve_endpoints(ax, [
           {"y_data": 141.0, "text": "Schur",  "color": "C1"},
           {"y_data":  27.0, "text": "Exact",  "color": "k"},
           {"y_data":   8.5, "text": "CLN-N",  "color": "C0"},
       ])
       fig.subplots_adjust(right=0.78)   # reserve right-margin room

  2. Programmatic best location:
       from radia_mcp.graph import find_best_legend_loc
       best, summary = find_best_legend_loc(ax)
       ax.legend(loc=best, frameon=False)

  3. Place outside the axes (often LESS efficient than 1):
       ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1),
                 frameon=False)
       fig.subplots_adjust(right=0.80)

Override with `check_legend_overlap=False` only when you have
manually inspected the figure at the final embed scale.
""",
        "no_legend_frame": """\
[no_legend_frame]

LAB RULE: NO BOX / FRAME AROUND THE LEGEND.

  WRONG:  ax.legend()                                # default frameon=True
  WRONG:  ax.legend(frameon=True)
  WRONG:  ax.legend(framealpha=0.9)                  # frame still rendered
  RIGHT:  ax.legend(frameon=False)
  RIGHT:  ax.legend(loc=find_best_legend_loc(ax)[0], frameon=False)

Why:
  - A box around the legend competes visually with the axis frame and
    the data lines.  In small embed (4-8 cm), the legend box itself
    becomes a noticeable rectangle that the eye reads as "another panel".
  - IEEE / IEEJ / Nature / Science figure conventions all expect a
    frameless legend.  A boxed legend reads as "PowerPoint slide" to
    journal reviewers.
  - Removing the frame increases the apparent axes-area fraction (less
    visual clutter in the plot area).
  - With direct-endpoint labels (`label_curve_endpoints`) you skip the
    legend entirely, which is the lab-preferred alternative for
    time-series / sweep plots (see `no_legend_overlap` rule).

  matplotlib default is `frameon=True` -- you MUST explicitly pass
  `frameon=False` every time, OR set it once in rcParams:
      matplotlib.rcParams['legend.frameon'] = False

  `paper_figure_recipe` and `paper_figure(...)` set this in rcParams for
  you; manual `ax.legend(...)` calls must still pass `frameon=False`.

How emit_paper_figure() detects it:
  After the figure is built, walk every axes' get_legend() and check
  legend.get_frame().get_visible() == False.  If ANY legend renders a
  frame, fail the gate with "legend N on axes M has frameon=True; pass
  frameon=False" as the message.

Override with `check_legend_frame=False` only for the rare case where
the legend NEEDS a frame to be readable (e.g. legend placed over a
busy heatmap and the box prevents data behind it from bleeding through).
""",
        "margins": """\
[margins]

measure_figure_efficiency(fig)['margin_breakdown'] = {
    'left':   <fraction of fig width left of leftmost axis>,
    'right':  <fraction of fig width right of rightmost axis>,
    'top':    <fraction of fig height above topmost axis>,
    'bottom': <fraction of fig height below bottommost axis>,
}

Reading guide (figure fractions):
  - any margin > 0.20  -> major waste, shave it
  - bottom typically the biggest because of x-tick labels + xlabel;
    a target of 0.15-0.18 is realistic for 8 pt font
  - left  similarly, target 0.10-0.16 for the ylabel + y-tick labels
  - top   target 0.02-0.05 (no title in figure, lab convention)
  - right target 0.01-0.04 (legend should NOT live in the right margin;
    use direct labels via label_curve_endpoints instead)
""",
        "units": """\
[units]

LAB STYLE (per IEEE Editorial Style Manual + IEEJ practice):
    units in PARENTHESES, NOT square brackets.

    correct:  'f (Hz)'  'B (T)'  'Temperature (K)'  'M (kA/m)'
    wrong:    'f [Hz]'  'B [T]'  'Temperature [K]'  'M [kA/m]'

ISO 80000 'f / Hz' is a distinct convention NOT used in IEEE/IEEJ
engineering papers.  Don't mix.

Math symbols italic, units upright:
    r'${\\it B}$ (T)'  r'${\\it H}$ (kA/m)'  r'$|{\\it Z}|$ ($\\Omega$)'
""",
        "font_embedding": """\
[font_embedding]

CRITICAL: matplotlib defaults to Type 3 PostScript fonts for PDF/EPS
output.  Type 3 fonts are RASTER glyphs that blur when zoomed and
fail many publisher pre-flight checks (IEEE, Elsevier).

Required rcParams (paper_figure() sets these for you):
    pdf.fonttype = 42     # Type 42 = TrueType, vector glyphs
    ps.fonttype  = 42

Verify on the output PDF:
    pdfinfo out.pdf | grep -i font   # any "Type 3" line is a failure
""",
        "multipanel": """\
[multipanel]

paper_figure(profile, nrows=R, ncols=C) auto-tunes the subplots_adjust
margins per (R, C):

  1x1  no delta -- single axes uses the profile's full margins
  1x2  -0.020 left, -0.02 wspace -- shared horizontal extent allows
       a slightly looser left margin to host the leftmost axis's
       ylabel without overlap, and tighter wspace
  1x3  -0.030 left, -0.04 wspace
  2x1  -0.020 bottom, -0.05 hspace -- shared x-axis amortizes the
       xlabel; can tighten interior gap
  2x2  -0.020 left + bottom, -0.02 wspace, -0.05 hspace

If you need non-standard layouts (e.g. 3x2, 1x4), pass them anyway
and tune via auto_tighten(fig, target_axes_fraction=0.80) after
plotting.  The auto_tighten loop will iteratively shrink each safe
margin in 0.005-step increments until labels would clip or target
is reached.

Direct labeling pattern for 1xC layouts:
  fig, axes = paper_figure('ieee_double_column', nrows=1, ncols=2,
                            panel_labels=True)
  # places (a), (b) in top-left of each panel automatically
""",
        "side_by_side": """\
[side_by_side]

TWO FIGURES IN 8 cm -> SIDE BY SIDE, each <= 4 cm.

When two graphs must share an 8 cm width, lay them out HORIZONTALLY
(1 row x 2 columns).  Do NOT stack them or shrink one; split the width.

GEOMETRY:
  - total embed width      : 8 cm
  - each sub-panel width    : <= 4 cm  (8 / 2; leave a small inter-panel
                              gap so each lands ~3.5-4 cm)
  - layout                  : 1 row x 2 cols (横並び, horizontal)

FONTS:
  - body / axis label / tick: 10 pt  (the absolute lab font rule -- do
                              NOT shrink it for the narrow 4 cm panel)
  - legend                  : 8-9 pt (REDUCED from 10 pt; a 10 pt
                              legend crowds a 4 cm panel).  This is the
                              ONE sanctioned exception to "legend = body".

LEGEND PLACEMENT:
  - the legend MUST NOT overlap the graph (curves / markers).
  - put it in an empty corner with frameon=False, OR use direct
    endpoint labels (label_curve_endpoints).
  - emit_paper_figure() rejects an overlapping legend (see
    no_legend_overlap topic).

TICKS:
  - keep sparse (~4-5 per axis); a 4 cm panel crowds easily.

HOW TO BUILD IT:
  fig, axes = paper_figure('ieee_double_column', nrows=1, ncols=2,
                            panel_labels=True)
  for ax in axes.ravel():
      ax.legend(loc=find_best_legend_loc(ax)[0], frameon=False,
                fontsize=8)            # 8-9 pt, not 10
  emit_paper_figure(fig, 'out', 'ieee_double_column', on_fail='raise')

Or use the size/font recipe directly:
  graph_size_for_target('digest_double_column_side_by_side')
  # -> 8 cm wide, font 10 pt, legend 8-9 pt
""",
    }
    q = (query or "all").strip().lower()
    if q == "all":
        return "\n\n".join(rules.values()) + """\

Recipe pipeline:
  1.  paper_figure_profiles('all')     -- pick profile + layout
  2.  paper_figure_recipe(profile=...) -- copy/paste recipe
  3.  Plot your data (NO ax.set_title!)
  4.  ax.legend(loc=find_best_legend_loc(ax)[0]) or
      label_curve_endpoints(ax, [...])  -- avoid legend-overlap
  5.  emit_paper_figure(fig, 'out', profile=..., on_fail='raise')
      -- gate refuses to save if any of:
           (a) ax.set_title() / fig.suptitle() set
           (b) legend overlaps a data line
           (c) axes_area / fig_area < 0.72
"""
    if q in rules:
        return rules[q]
    return (
        f"Unknown topic {query!r}. Available: "
        f"{', '.join(['all'] + list(rules))}"
    )


# ============================================================
# Self-introspection (uniform with other radia_mcp servers)
# ============================================================

register_status_tool(
    mcp,
    server_name="mcp-server-graph",
    description="Sugahara Lab publication-figure style guide: "
                "IEEE / IEEJ font/size profiles, MATLAB + Matplotlib "
                "snippets, lab style rules (units in parentheses, no "
                "in-figure title, Times New Roman serif).",
    subpackage="radia_mcp.graph",
    related_servers=["mathematica", "literature-index"],
    # matplotlib is needed for the python helper functions
    # (apply_lab_style / lab_figsize / etc.), NOT for the MCP tools.
    # Server loads without it.
    optional_deps=["matplotlib"],
)


# ============================================================
# Entry point
# ============================================================

def main():
    if "--selftest" in sys.argv:
        print("graph MCP server self-test:")
        print(f"  Registered tools ({len(_REGISTERED)}):")
        for name in sorted(_REGISTERED):
            print(f"    - {name}")
        # Smoke-test both legacy tools across every profile so a typo
        # in the _PROFILES dict immediately surfaces in CI.
        profiles = sorted(_tools._PROFILES.keys())
        print(f"  Lab-style profiles ({len(profiles)}):")
        for prof in profiles:
            guide = _tools.graph_style_guide(prof)
            size = _tools.graph_size_for_target(prof)
            assert len(guide) > 200, (
                f"graph_style_guide({prof!r}) returned only "
                f"{len(guide)} chars (suspiciously short)"
            )
            assert "Unknown target" not in guide, (
                f"graph_style_guide({prof!r}) reported unknown target"
            )
            assert "Unknown target" not in size, (
                f"graph_size_for_target({prof!r}) reported unknown target"
            )
            print(f"    {prof:42s} guide={len(guide):4d} ch, "
                  f"size={len(size):4d} ch")
        # Sanity: 'all' returns everything
        full = _tools.graph_style_guide("all")
        assert all(p in full for p in profiles), \
            "graph_style_guide('all') missing some profiles"
        # Sanity: unknown target returns help text (not crash)
        unk = _tools.graph_style_guide("not-a-real-profile")
        assert "Unknown target" in unk and "Valid:" in unk
        print(f"  graph_style_guide('all')               -> "
              f"{len(full):5d} chars")
        print(f"  graph_style_guide('not-a-real-profile') -> "
              f"unknown-target help text emitted")

        # --- paper-figure tier (v0.78.0) ---
        print(f"  Paper-quality profiles ({len(_paper.PROFILES)}):")
        for pname, p in _paper.PROFILES.items():
            print(f"    {pname:24s} {p.width_mm:6.2f} mm  "
                  f"font={p.font_pt:4.1f} pt  L/R/T/B="
                  f"{p.margin_left:.3f}/{p.margin_right:.3f}/"
                  f"{p.margin_top:.3f}/{p.margin_bottom:.3f}")
        prof_list = paper_figure_profiles("all")
        assert "ieee_double_column" in prof_list
        assert "ieej_single_column" in prof_list
        assert "igte_digest_double" in prof_list
        print(f"  paper_figure_profiles('all')      -> "
              f"{len(prof_list):5d} chars")

        # Recipe smoke-test: every layout shape generates a non-empty
        # recipe that mentions the profile name + the gate.
        layouts = [(1, 1), (1, 2), (2, 1), (2, 2), (1, 3)]
        for nrows, ncols in layouts:
            rec = paper_figure_recipe(
                profile="ieee_double_column",
                nrows=nrows, ncols=ncols,
                panel_labels=(nrows * ncols > 1),
            )
            assert "paper_figure(" in rec
            assert "emit_paper_figure(" in rec
            assert "ieee_double_column" in rec
            assert "min_axes_fraction=0.72" in rec
            print(f"    recipe[{nrows}x{ncols}]       -> "
                  f"{len(rec):5d} chars")
        # Quality-rules tool
        rules = paper_figure_quality_rules("all")
        for sec in ("efficiency", "margins", "units",
                    "font_embedding", "multipanel"):
            assert f"[{sec}]" in rules, f"section {sec!r} missing"
        print(f"  paper_figure_quality_rules('all') -> "
              f"{len(rules):5d} chars")

        # Try importing matplotlib + smoke-test the runtime helpers
        # (paper_figure / measure_figure_efficiency).  Skip if mpl is
        # absent (the MCP tools above already loaded without mpl).
        try:
            import matplotlib
            matplotlib.use("Agg")
            from radia_mcp.graph import (
                paper_figure as _pf,
                measure_figure_efficiency as _meas,
                emit_paper_figure as _emit,
                auto_tighten as _at,
            )
            fig, axes = _pf("ieee_double_column", nrows=1, ncols=2)
            for ax in axes.flat:
                ax.plot([0, 1], [0, 1])
                ax.set_xlabel(r"$f$ (Hz)")
                ax.set_ylabel(r"$|Z|$ ($\Omega$)")
            m = _meas(fig)
            print(f"  paper_figure(ieee_double_column, 1, 2):")
            print(f"    fig {m['fig_size_inches'][0]:.2f} x "
                  f"{m['fig_size_inches'][1]:.2f} in, "
                  f"axes_area_fraction = {m['axes_area_fraction']:.3f}")
            assert m["axes_area_fraction"] > 0.55, (
                f"baseline IEEE 1x2 fell below 0.55 "
                f"({m['axes_area_fraction']:.3f}) — profile margins drifted"
            )
            import matplotlib.pyplot as plt
            plt.close(fig)
        except ImportError:
            print("  [skipped runtime smoke-test: matplotlib not installed]")

        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
