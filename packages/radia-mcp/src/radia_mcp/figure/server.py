"""
Graph MCP Server (radia_mcp.figure)

Surfaces the Sugahara Lab publication-figure style guide via two MCP
tools:

    figure_style_guide(target='all'|'paper_single_column'|...)
        Returns the lab-standard rules (Times New Roman, units in
        parentheses, no in-figure title, font sizes per embed width,
        IEEE/IEEJ conventions) plus the numeric profile for the
        selected target.

    figure_size_for_target(target, embed_width_cm=None)
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
carry.  They live in `radia_mcp.figure.tools` for direct import by
analysis scripts.

Usage:
    mcp-server-figure              # start MCP server (stdio transport)
    mcp-server-figure --selftest   # smoke check every profile in both tools
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool

from . import tools as _tools
from . import _paper_figure as _paper

mcp = FastMCP("mcp-server-figure")


# ============================================================
# Tool registration
# ============================================================
# Auto-register the MCP-callable subset (graph_*).  The Python helpers
# (apply_lab_style etc.) have matplotlib-typed signatures that don't
# round-trip through MCP-JSON, so we skip them deliberately.
_REGISTERED: list[str] = []
for _name in dir(_tools):
    if _name.startswith("figure_"):
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
#   from radia_mcp.figure import paper_figure, emit_paper_figure
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
        lines = ["radia_mcp.figure paper-figure profiles:", ""]
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
      1. imports paper_figure + emit_paper_figure from radia_mcp.figure
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

from radia_mcp.figure import paper_figure, emit_paper_figure

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
        'tikz_export'     - MATLAB -> matlab2tikz -> LaTeX TikZ:
                            preferred over PDF includegraphics for
                            LaTeX papers (font / math matches body)
        'export_targets'  - format matrix: vector PDF (paper) / TikZ
                            (LaTeX) / EMF (Word/PowerPoint, MATLAB
                            -dmeta) / PNG 400 dpi (draft/web) -- which
                            format for which venue, from real lab scripts
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

THE ON-PAGE RATIO both cases obey: 10 pt at 8 cm width = 1.25 pt per
cm of figure width.  What must be 10 pt is the text ON THE PRINTED PAGE.

WIDER COLUMNS embedded at 100% (16-18 cm \\figure* / digest full-width,
authored AT that width in matplotlib):
  The font stays 10 pt.  Wider columns get a BIGGER AXES BOX, NOT
  bigger text.  The text is ABSOLUTE for a directly-embedded figure.

  Wrong intuition: "the figure is twice as wide, so font should be
                    twice as big (20 pt)" -- WRONG for a 100%-embedded
                    matplotlib figure.  At 18 cm a 20-pt font reads as
                    a billboard, not a paper.
  Right intuition: a 10-pt font on 18 cm = 10 pt on 8 cm.  Same
                    readability; bigger axes for the data.

MATLAB (authored OVERSIZED, then \\includegraphics-DOWNSCALED):
  The lab authors MATLAB plots at ~2x (16 cm) and embeds them at the
  8 cm column -- a 0.5x downscale that HALVES the on-page font.  So
  author at 20 pt @ 16 cm -> it lands at 10 pt @ 8 cm on the page.
  This is the `matlab_oversized_for_8cm_embed` profile.  The
  discriminator is the EMBED SCALE FACTOR, not the tool: author-at-
  embed-size -> 10 pt @ 8 cm; author-at-2x-then-downscale -> 20 pt
  @ 16 cm.  Both deliver 10 pt on the printed page.

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
       from radia_mcp.figure import label_curve_endpoints
       label_curve_endpoints(ax, [
           {"y_data": 141.0, "text": "Schur",  "color": "C1"},
           {"y_data":  27.0, "text": "Exact",  "color": "k"},
           {"y_data":   8.5, "text": "CLN-N",  "color": "C0"},
       ])
       fig.subplots_adjust(right=0.78)   # reserve right-margin room

  2. Programmatic best location:
       from radia_mcp.figure import find_best_legend_loc
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
  figure_size_for_target('digest_double_column_side_by_side')
  # -> 8 cm wide, font 10 pt, legend 8-9 pt
""",
        "tikz_export": """\
[tikz_export]

MATLAB FIGURE -> TikZ (via matlab2tikz) -> LaTeX paper.

For figures rendered in MATLAB and embedded in an IEEE / IEEJ / IGTE
LaTeX paper, the lab-preferred export path is TikZ via matlab2tikz,
NOT saveas('fig.pdf') / exportgraphics('fig.pdf').

WHY TikZ beats PDF here:
  - axis / tick / legend text INHERITS the paper's LaTeX font
    (Times New Roman for IEEE/IEEJ).  PDF includegraphics bakes in
    whatever font MATLAB happened to render with -- never an exact
    match for body text.
  - inline math like $\\sigma_{xy}$ in axis labels renders in the
    paper's MATH font, not MATLAB's LaTeX-interpreter approximation.
  - fully vector, editable in .tex after export (tweak labels,
    colors, ticks without re-running MATLAB).
  - pgfplots scales the figure to \\columnwidth / \\textwidth -- one
    .tikz file works for both single-column and double-column
    layouts.

WHEN TO STAY ON PDF (TikZ exceptions):
  - heatmaps / pcolor / large image overlays (raster, slow in TikZ).
  - >10000 plot points without cleanfigure pre-processing
    (LaTeX compile time blows up).
  - complex 3-D scenes that pgfplots struggles to reproduce.
  - photographs.
  -> for these, use exportgraphics PDF or a hybrid (raster the heavy
     layer, TikZ-overlay the axes/labels).

THE LAB RECIPE (parameterise with the lab profile):

    cleanfigure('targetResolution', 300);   % decimate dense data
    matlab2tikz('fig/result.tikz', ...
        'width', '\\figureWidth', ...        % LaTeX-side \\setlength
        'height', '\\figureHeight', ...      %   to \\columnwidth
        'parseStrings',  false, ...          % keep your $\\LaTeX$
        'showInfo',      false, ...
        'showWarnings',  false, ...
        'standalone',    false);             % embed; not stand-alone

LaTeX-side prelude (one-time):

    \\usepackage{pgfplots}
    \\pgfplotsset{compat=1.18}
    \\newlength\\figureWidth   \\setlength\\figureWidth{\\columnwidth}
    \\newlength\\figureHeight  \\setlength\\figureHeight{6cm}
    \\input{fig/result.tikz}

GETTING THE EXACT LAB-PROFILE RECIPE:

    figure_matlab2tikz_recipe(target='paper_single_column')
        -> ready-to-paste MATLAB recipe sized for IEEE single column
           (88.9 mm) with the lab's Times New Roman / 10 pt / 0.7 pt
           axis-linewidth / sparse-ticks defaults baked in.

    figure_matlab2tikz_recipe(target='paper_double_column')
        -> double-column (181 mm) wide figure.

    figure_matlab2tikz_recipe(target='digest_double_column_side_by_side')
        -> two-panel-in-8cm digest layout; remember to clip the
           legend font down to 8-9 pt for 4 cm sub-panels
           (see `side_by_side` topic).

PRE-FLIGHT one-time install:

    % After git clone https://github.com/matlab2tikz/matlab2tikz
    addpath(genpath('<install-dir>/matlab2tikz/src'));

CAVEATS:
  - matlab2tikz currently supports up to MATLAB R2024a-ish; very new
    graphics objects (e.g. some R2024b chart types) may export as
    rasterised fallbacks.
  - For very wide / dense plots, increase
    cleanfigure(..., 'targetResolution', 600) and tune the
    'minimumPointsDistance' option.
""",
        "export_targets": """\
[export_targets]

LAB EXPORT-FORMAT MATRIX (extracted 2026-06 from the lab's real MATLAB
+ matplotlib scripts on S:, including the FEMM folder).  The correct
format depends on WHERE the figure is embedded -- the lab uses three
in practice, not just paper PDF:

  Venue / embed              Format          How
  -------------------------  --------------  --------------------------
  IEEE/IEEJ/IGTE paper       vector PDF      emit_paper_figure(...) ->
   (LaTeX includegraphics)    (Type-42)       .pdf  (this server default)
  LaTeX paper, exact font    TikZ            matlab2tikz / see the
                                             tikz_export topic
  Word / PowerPoint          EMF (vector)    MATLAB exportgraphics(gcf,
   (Office embed)                            'f.emf','ContentType',
                                             'vector') -- stays vector in
                                             Office + prints crisp.
                                             (figure_office_export_recipe)
  Draft / web / slide        PNG 400-600 dpi matplotlib savefig(dpi=400)
   quick-look                 (raster)        / MATLAB exportgraphics(gcf,
                                             'f.png','Resolution',400).

MATLAB EXPORT: use exportgraphics (R2020a+), NOT print.  The old
print('-dmeta','f.emf') / print('-dpng') idiom (seen in legacy FEMM
scripts) is DEPRECATED -- exportgraphics is the modern, supported path
and keeps EMF vector.  Ready-to-paste recipes:
  - figure_office_export_recipe()  -> MATLAB exportgraphics EMF + PNG
  - figure_everyday_recipe()       -> matplotlib everyday analysis figure

OBSERVED LAB HABIT (from S: scripts):
  - FEMM MATLAB (legacy): set(gca,'FontName','Times'); xlabel('{\\it X}
                   (m)'); print('-dmeta','f.emf');  % print is DEPRECATED
                   -> modern: exportgraphics(...,'ContentType','vector')
  - COMSOL/CoreformCubit matplotlib: figsize=(3,4), dpi=400, Times New
                   Roman 10 pt, inward ticks, savefig PNG.

RULE: NEVER embed a raster PNG in a CAMERA-READY paper -- re-render to
vector PDF (or TikZ).  PNG/EMF are for drafts + Office.  EMF (not PNG)
is the right Office format because it keeps the text vector inside Word.

VERIFICATION (2026-06): the lab's real S: scripts AGREE with this
server's core rules -- Times/TNR, italic variable + unit in PARENTHESES
('{\\it X} (m)'), box-on + inward ticks, frameless legend, NO in-figure
title.  The everyday matplotlib analysis style additionally uses a
subtle two-level grid (major dotted, minor dashed, light-gray
~gainsboro) with minor ticks on; the paper profiles deliberately keep
the lighter single dotted grid (grid.alpha 0.4) to stay reviewer-clean.
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


@mcp.tool()
def figure_design_principles(topic: str = "all") -> str:
    """The figure-MAKING (作図, *sakuzu*) DESIGN canon, distilled from the
    authoritative external scientific-visualization literature (GitHub repos +
    the canonical papers), each topic attributed to its source.

    This is the DESIGN layer that sits ABOVE two other layers:
      * the lab MECHANICS / save-time gates -> paper_figure_quality_rules
        (no-title, frameless legend, units-in-parens, axes-efficiency,
        no-overlap, Type-42, Okabe-Ito);
      * the data-PLOTTING (グラフ, the line/scatter/bar itself) ->
        radia_mcp.chart2d + paper_figure()/emit_paper_figure().

    The user's distinction: 作図 (figure DESIGN -- message, encoding, colour,
    labelling, composition, medium) decided FIRST, vs グラフ (data PLOTTING)
    drawn after.  This tool carries the 作図 knowledge.

    Topics:
        'all'                - everything
        'ten_rules'          - Rougier, Droettboom & Bourne (2014), Ten Simple Rules
        'perception'         - Cleveland & McGill (1984) graphical-perception ranking
        'color'              - perceptually-uniform maps (Crameri 2020) + Okabe-Ito
        'chartjunk'          - Tufte data-ink ratio; strip non-data ink
        'direct_labeling'    - label curves directly > legend (adjustText, ...)
        'defaults'           - "do not trust the defaults"; journal styles
        'external_resources' - where to learn 作図 (books, repos, papers)
        'sakuzu_vs_graph'    - the 作図 vs グラフ split + the radia_mcp layering
    """
    principles = {
        "ten_rules": """\
[ten_rules]  -- Rougier, Droettboom & Bourne (2014), "Ten Simple Rules for Better
Figures", PLOS Comput. Biol. 10(9):e1003833 (the canonical principle checklist;
lead author wrote matplotlib's scientific-visualization reference book).

  1.  Know your audience       - peers vs students vs public set the detail level.
  2.  Identify your message    - one figure = ONE message; design backward from it.
  3.  Adapt to the medium      - paper / slide / poster need DIFFERENT figures
                                 (font, detail, size); never reuse a paper figure
                                 verbatim on a slide.
  4.  Captions are not optional - the caption carries what the graphic cannot;
                                 this is WHY the lab forbids in-figure titles
                                 (paper_figure_quality_rules('no_title_in_figure')).
  5.  Do not trust the defaults - matplotlib / MATLAB defaults are NOT publication
                                 ready (see the 'defaults' topic).
  6.  Use colour effectively   - colour is a DATA channel, not decoration ('color').
  7.  Do not mislead           - honest axes/scales; pick the encoding the eye reads
                                 accurately (see 'perception').
  8.  Avoid chartjunk          - every drop of ink should be data ('chartjunk').
  9.  Message trumps beauty    - clarity / readability over prettiness.
  10. Get the right tool       - match the tool to the job (vector for line art,
                                 raster for images).
""",
        "perception": """\
[perception]  -- Cleveland & McGill (1984), "Graphical Perception: Theory,
Experimentation, and Application...", J. Amer. Statist. Assoc. 79(387):531
(the experimental ranking of how ACCURATELY the eye decodes each channel).

  Most -> least accurate quantitative decoding:
    1. position along a COMMON scale   (scatter, dot plot, aligned bars)
    2. position on non-aligned scales
    3. length                          (unaligned bars)
    4. angle / slope                   (line slope, pie wedges)
    5. area                            (bubble size)
    6. volume / colour-hue / saturation

  CONSEQUENCE for chart choice (the evidence behind rule 7 "do not mislead"):
    - to COMPARE magnitudes, encode as POSITION or LENGTH, never area or colour.
      Prefer a dot / bar / line over a pie (angle), a bubble (area) or a heat-tile
      (colour) when the reader must read values accurately.
    - reserve colour-HUE for CATEGORIES and colour-LUMINANCE for ORDERED fields
      (colour is a weak quantitative channel -- good for "which", poor for
      "how much").
    - a convergence / defect plot (e.g. the lab act2_1x figures) is rightly a
      log-y POSITION encoding -- the most accurate channel.
""",
        "color": """\
[color]  -- Crameri, Shephard & Heron (2020), "The misuse of colour in science
communication", Nature Commun. 11:5444 + Okabe & Ito (2008) / Wong (2011),
Nature Methods 8:441.

  THREE colormap CLASSES -- choose by DATA type, never by taste:
    - QUALITATIVE (categories; lines/markers): the Okabe-Ito 8-colour colourblind
      -safe set (lab line default; paper_figure_quality_rules('colorblind_safe')).
    - SEQUENTIAL (ordered 0->max field, e.g. |B|): a PERCEPTUALLY-UNIFORM map --
      viridis / cividis / magma (matplotlib) or Crameri batlow / lajolla.
    - DIVERGING (signed field about 0, e.g. +/- error): vik / roma / coolwarm
      (perceptually-uniform, symmetric about a neutral midpoint).

  NEVER use jet / rainbow / hsv for data:
    Crameri 2020 shows a physically-built (not perceptually-built) map ADDS
    artificial boundaries where its luminance jumps and HIDES variation where
    luminance is flat -- it DISTORTS the data -- AND is unreadable to ~4-8% of
    readers with colour-vision deficiency AND fails greyscale print.  Perceptually
    -uniform maps are perceptually ordered, colourblind-friendly, greyscale-safe.

  In radia_mcp: lines already default to Okabe-Ito (paper_figure); for FIELD plots
  (chart2d pcolormesh / contourf) pass cmap='viridis' (or a Crameri map via the
  `cmcrameri` package), never matplotlib's historical 'jet'.
""",
        "chartjunk": """\
[chartjunk]  -- Tufte (1983), "The Visual Display of Quantitative Information"
(the data-ink ratio) + Rougier rule 8.

  data-ink ratio = ink-encoding-data / total-ink.  Maximise it: every stroke
  should carry information; erase the rest.

  REMOVE the non-data ink (the lab gates already enforce most of this):
    - 3-D bars / pie / extruded charts (3-D adds area+volume distortion, the
      worst-decoded channels -- see 'perception').
    - the LEGEND BOX (frameon=False; paper_figure_quality_rules('no_legend_frame')).
    - the in-figure TITLE (-> the caption; 'no_title_in_figure').
    - heavy / redundant gridlines, drop shadows, gradient fills, dense ticks.
    - idle whitespace margins (the axes-area gate, 'efficiency').

  KEEP: the data marks, the axis frame + sparse ticks, axis labels with
  units-in-parentheses, and direct curve labels.  A lab paper figure is
  deliberately spare -- the lab principle "情報がなく無駄はやめる"
  (drop anything carrying no information).
""",
        "direct_labeling": """\
[direct_labeling]  -- label the curve WHERE IT IS, not in a legend the eye must
round-trip to (Rougier rules 8/9; the R `ggrepel` lineage).

  WHY: a legend forces the reader to match colour -> name -> curve repeatedly; an
  inline label at the curve's end is read in place.  For time-series / sweeps /
  convergence plots with a handful of curves, direct labels beat a legend.

  EXTERNAL tools (the GitHub lineage to learn from):
    - adjustText  (github.com/Phlya/adjustText) -- `ggrepel`-for-Python; iteratively
      repositions many text labels to remove overlap with points and each other.
      `pip install adjustText`; `from adjustText import adjust_text`.
    - matplotlib-label-lines (github.com/cphyc/matplotlib-label-lines) -- places
      each line's label ON the line: `labelLines(ax.get_lines())`.

  LAB tools (radia_mcp.figure.tools -- the gate-aware versions):
    - label_curve_endpoints(ax, [{"y_data":.., "text":..}], side='right') +
      fig.subplots_adjust(right=0.78) -- right-margin direct labels.
    - find_best_legend_loc(ax) -- if a legend IS used, pick the least-overlapping
      location; emit_paper_figure() then HARD-FAILS on any legend-data overlap
      (paper_figure_quality_rules('no_legend_overlap')).
""",
        "defaults": """\
[defaults]  -- Rougier rule 5, "Do Not Trust the Defaults".

  Out-of-the-box matplotlib / MATLAB is NOT publication-ready:
    - the tab10 cycle has a confusable red+orange (colourblind-unsafe);
    - pdf.fonttype defaults to Type-3 raster glyphs (fail IEEE / Elsevier preflight);
    - legends render a BOX; titles default on; jet was the historical image cmap.

  Two ways to fix it by STARTING from a journal style, not the default:
    - SciencePlots (github.com/garrettj403/SciencePlots) -- `pip install
      SciencePlots`; `import scienceplots; plt.style.use(['science','ieee'])`
      (or 'nature').  Sets the column width, serif/sans per journal, ticks-in,
      tight rcParams; also ships CJK font styles (e.g. cjk-jp) for Japanese labels.
    - the LAB path: radia_mcp.figure.paper_figure(profile=...) -- bakes in the exact
      IEEE / IEEJ / IGTE column geometry + 10pt@8cm + Okabe-Ito + Type-42 + frameless
      legend + the no-title / no-overlap / efficiency GATES (emit_paper_figure).
      Prefer this for lab papers; SciencePlots is the lighter general-purpose option.
""",
        "external_resources": """\
[external_resources]  -- where to LEARN 作図 (figure-making), curated 2026-06.

  BOOKS / PRINCIPLES
    - Rougier, "Scientific Visualization: Python + Matplotlib" (2021) -- FREE,
      github.com/rougier/scientific-visualization-book.  THE reference: figure
      anatomy, coordinate systems, typography, colour, design rules, layout.
    - Rougier, Droettboom & Bourne, "Ten Simple Rules for Better Figures" (2014),
      PLOS Comput. Biol. 10(9):e1003833 -- the principle checklist.
    - Tufte, "The Visual Display of Quantitative Information" (1983) -- data-ink.
    - Cleveland & McGill, "Graphical Perception" (1984), JASA 79(387):531 -- the
      channel-accuracy ranking.
    - Wong, "Points of view: Color blindness" (2011), Nature Methods 8:441 -- the
      Okabe-Ito palette.
    - Crameri, Shephard & Heron, "The misuse of colour in science communication"
      (2020), Nature Commun. 11:5444 -- perceptually-uniform colour maps.

  TOOLS / REPOS (matplotlib ecosystem)
    - SciencePlots   github.com/garrettj403/SciencePlots   -- journal style sheets.
    - cmcrameri      github.com/callumrollo/cmcrameri       -- Crameri colour maps.
    - adjustText     github.com/Phlya/adjustText            -- auto non-overlap labels.
    - matplotlib-label-lines github.com/cphyc/matplotlib-label-lines -- inline labels.
    - the matplotlib official "Choosing Colormaps" + "Annotations" user guides.

  These external sources inform the lab's own radia_mcp.figure (paper_figure, the
  quality gates, label_curve_endpoints, the Okabe-Ito default).
""",
        "sakuzu_vs_graph": """\
[sakuzu_vs_graph]  -- 作図 (figure DESIGN) vs グラフ (data PLOTTING): the split the
radia_mcp toolchain mirrors.

  作図 (sakuzu) = the DESIGN decisions made BEFORE / AROUND the plot: what message
    (rule 2), what to encode + which channel (perception), colour's role (color),
    what to strip (chartjunk), how to label (direct_labeling), medium adaptation
    (paper vs slide), the caption.
      -> THIS tool (figure_design_principles) for the canon, plus
         paper_figure_quality_rules for the lab MECHANICS / gates that ENFORCE the
         design (no-title, frameless, units-in-parens, axes-efficiency, no-overlap,
         Type-42, Okabe-Ito).

  グラフ (graph) = the data PLOTTING itself, the line / scatter / bar:
      -> radia_mcp.chart2d (22 chart types) drawn on a paper_figure() canvas, then
         emit_paper_figure() gates the result.

  ORDER: decide the 作図 (message, encoding, colour, labels) FIRST, then draw the
  グラフ on a paper_figure() canvas, then let emit_paper_figure() gate it against
  the design rules.  A good グラフ on a bad 作図 still fails review -- design first.
""",
    }
    q = (topic or "all").strip().lower()
    if q == "all":
        return "\n\n".join(principles.values()) + """\

Layering (see 'sakuzu_vs_graph'):
  figure_design_principles   -- the 作図 DESIGN canon (this tool)
  paper_figure_quality_rules -- the lab MECHANICS + save-time gates
  radia_mcp.chart2d / paper_figure -- the グラフ data-plotting + journal canvas
"""
    if q in principles:
        return principles[q]
    return (
        f"Unknown topic {topic!r}. Available: "
        f"{', '.join(['all'] + list(principles))}"
    )


@mcp.tool()
def figure_diagram_recipes(topic: str = "all") -> str:
    r"""Flowchart + conceptual/schematic DIAGRAM recipes -- the diagram-DRAWING skill
    (distinct from data-PLOTTING グラフ and from the general 作図 design canon in
    figure_design_principles).

    COMPLEMENTS figure_tikz_recipe (which owns the general TikZ schematic / geometry /
    flux-path template + PGFPlots + externalize + matlab2tikz).  THIS tool adds the
    FLOWCHART + CONCEPTUAL-diagram skill and the multi-tool ecosystem: the TikZ flowchart
    idiom (shapes.geometric node styles), GRAPHVIZ/DOT AUTO-layout (the big complement --
    figure_tikz_recipe is manual-coordinate TikZ only), schemdraw / Mermaid, tool
    selection, and diagram DESIGN conventions (ISO 5807 symbols, flow direction, crossing
    minimisation).

    Topics:
        'all'                - everything
        'tool_selection'     - TikZ vs Graphviz vs schemdraw vs Mermaid vs draw.io
        'tikz_flowchart'     - the TikZ shapes.geometric flowchart idiom (ready template)
        'graphviz'           - DOT auto-layout (digraph, rankdir, clusters, engines)
        'concept_diagram'    - conceptual / block-relationship diagrams (architecture style)
        'design'             - ISO 5807 symbols + flow direction + crossing/alignment rules
        'external_resources' - manuals, galleries, standards
    """
    recipes = {
        "tool_selection": r"""[tool_selection] -- pick the diagram tool by WHERE it goes + HOW it is laid out.

  Tool          Best for                              Layout     Source     Paper-native?
  ------------  ------------------------------------  ---------  ---------  -------------
  TikZ          paper figures: schematics, flux       MANUAL     .tex       YES (font + math
                paths, small flowcharts, concept                            match the body)
                blocks -- coordinate-precise.
  Graphviz/DOT  large flowcharts, DAGs, dependency     AUTO       .dot       via PDF/SVG
                / state graphs -- when manual node     (dot/...)             include (font
                placement is tedious.                                        differs)
  schemdraw     programmatic flowcharts / circuits     semi-auto  .py        via PDF/SVG
                from Python (loop over data).          (Python)
  Mermaid       README / docs / web flowcharts,        AUTO       md fence   NO (docs only)
                quick sketches in Markdown.
  draw.io       one-off GUI diagrams (NOT diffable --  GUI        .drawio    export PDF
                avoid for reproducible repo figures).

  LAB DEFAULT for a PAPER diagram: TikZ -- text-source/diffable, labels inherit the paper's
  Times/newtx font + math (see figure_tikz_recipe for the general schematic; 'tikz_flowchart'
  here for the flowchart idiom).  Reach for Graphviz when the graph is big or a hierarchy/DAG
  you do NOT want to place by hand, then include the rendered PDF.  Mermaid = repo READMEs
  only (never camera-ready).
""",
        "tikz_flowchart": r"""[tikz_flowchart] -- the TikZ shapes.geometric flowchart idiom (ISO-5807 node shapes,
auto-spaced with `positioning`).  For a non-flow schematic (geometry, flux, BC) use
figure_tikz_recipe('schematic') instead.

\documentclass[tikz,border=2mm]{standalone}
\usepackage{newtxtext,newtxmath}                 % match IEEE/IEEJ body font
\usetikzlibrary{shapes.geometric, shapes.misc, arrows.meta, positioning}  % shapes.misc = rounded rectangle
\begin{document}
\begin{tikzpicture}[
  node distance=8mm and 14mm, font=\footnotesize,
  start/.style   ={rounded rectangle, draw, fill=black!5, minimum height=7mm, inner xsep=3mm},
  process/.style ={rectangle, draw, fill=blue!4,  minimum height=7mm, text width=24mm, align=center},
  decision/.style={diamond, draw, fill=orange!12, aspect=2, inner sep=1pt, align=center},
  io/.style      ={trapezium, trapezium left angle=70, trapezium right angle=110,
                   draw, fill=black!4, minimum height=7mm},
  arrow/.style   ={-{Latex[length=2mm]}, semithick},
]
  \node[start]                    (a) {start};
  \node[io,       below=of a]     (b) {read .vol};
  \node[process,  below=of b]     (c) {assemble DtN};
  \node[decision, below=of c]     (d) {$p\ge n$?};
  \node[process,  below=of d]     (e) {refine order};
  \node[start,    right=24mm of d](f) {done};
  \draw[arrow] (a)--(b); \draw[arrow] (b)--(c); \draw[arrow] (c)--(d);
  \draw[arrow] (d)-- node[left]{no} (e);
  \draw[arrow] (e.west) -- ++(-7mm,0) |- (c.west);     % feedback loop, routed orthogonally
  \draw[arrow] (d)-- node[above]{yes} (f);
\end{tikzpicture}
\end{document}

KEYS: the SHAPE encodes the ISO-5807 meaning (see 'design'); `positioning` (`below=of`,
`right=of` + `node distance=A and B`) AUTO-spaces -- never hand-tune (x,y) for a flowchart;
route feedback edges orthogonally with `|-` / `-|`.  Compile standalone -> PDF, or paste the
tikzpicture into the paper and embed at the column width.
""",
        "graphviz": r"""[graphviz] -- DOT AUTO-layout: let the engine place the nodes.  Use when the graph is
large or a hierarchy/DAG you do NOT want to position by hand (the big complement to manual
TikZ; figure_tikz_recipe has no auto-layout).

  flow.dot:
    digraph G {
      rankdir=TB;                                  // TB top-down | LR left-right
      node [shape=box, style=rounded, fontname="Times", fontsize=10];
      edge [fontname="Times", fontsize=9];
      start [shape=stadium, label="start"];
      read  [shape=parallelogram, label="read .vol"];
      asm   [label="assemble DtN"];
      chk   [shape=diamond, label="p >= n ?"];
      ref   [label="refine order"];
      done  [shape=stadium, label="done"];
      start -> read -> asm -> chk;
      chk -> ref  [label="no"];
      ref -> asm  [constraint=false];              // feedback: do NOT affect ranking
      chk -> done [label="yes"];
      subgraph cluster_solve { label="solve loop"; style=dashed; asm; chk; ref; }
    }

  render:  dot -Tpdf flow.dot -o flow.pdf          # vector, for paper include
           dot -Tsvg flow.dot -o flow.svg          # web / docs

  LAYOUT ENGINES (pick by graph shape):
    dot    layered / hierarchical -> FLOWCHARTS, DAGs, call graphs   (the default choice)
    neato  spring model           -> small undirected relationship graphs
    fdp    force-directed         -> larger undirected / clustered graphs
    circo  circular               -> ring / cyclic topologies
    twopi  radial                 -> trees around a centre

  TIPS: `rankdir=LR` for wide-short page fits; `constraint=false` stops a feedback edge from
  distorting the ranking; `subgraph cluster_*` draws a labelled box round a group; set
  fontname="Times" to approach the paper body (still NOT an exact match -- for exact font
  use TikZ, see 'tool_selection').
""",
        "concept_diagram": r"""[concept_diagram] -- conceptual / block-relationship diagrams (architecture, data flow,
"X consumes Y"): boxes + LABELLED arrows + optional grouping; not a strict process flow.

  TikZ (paper, exact font; `fit`+`backgrounds` draw the group box):
    \usetikzlibrary{positioning, fit, backgrounds, arrows.meta}
    \begin{tikzpicture}[font=\footnotesize, >={Latex[length=2mm]},
      blk/.style={rectangle, draw, rounded corners, fill=black!4,
                  minimum height=8mm, text width=22mm, align=center}]
      \node[blk] (cubit) {Cubit\\hex mesh};
      \node[blk, right=14mm of cubit] (vol) {.vol};
      \node[blk, right=14mm of vol]   (ng)  {NGSolve\\FEM};
      \draw[->] (cubit) -- node[above]{export} (vol);
      \draw[->] (vol)   -- node[above]{Mesh()} (ng);
      \begin{scope}[on background layer]
        \node[draw=blue!40, dashed, rounded corners, fit=(vol)(ng),
              inner sep=3mm, label=below:{computation}] {};
      \end{scope}
    \end{tikzpicture}

  Graphviz alternative (auto-layout, good when there are many blocks):
    digraph { rankdir=LR; node[shape=box,style=rounded,fontname=Times];
      cubit->vol[label=export]; vol->ng[label="Mesh()"]; }

  The lab CLAUDE.md ASCII box-diagrams (the 4-Layer panel architecture, the accelerator-magnet
  pipeline) convert directly: one box per ASCII box, one arrow per `->`.  Keep ONE flow
  direction; group with a dashed `fit` box (TikZ) or a `cluster` (DOT).
""",
        "design": r"""[design] -- flowchart / diagram design rules (ISO 5807 symbols + layout craft; the
Tufte/Rougier canon in figure_design_principles applies here too).

  SHAPE = MEANING (ISO 5807:1985 flowchart symbols -- keep them consistent):
    terminator (start / end)    stadium / rounded rectangle
    process / action            rectangle
    decision / branch           diamond   (label EVERY outgoing edge: yes / no)
    input / output (data)       parallelogram
    predefined process (sub)    rectangle with double side bars
    connector (off-page / loop) small circle

  LAYOUT:
    - ONE dominant flow direction: top-to-bottom OR left-to-right, never both.
    - align nodes on a grid (TikZ `node distance` / DOT ranks); ragged placement reads
      as careless.
    - MINIMISE edge crossings; route feedback / loop edges orthogonally around the side
      (TikZ `-|` / `|-`; DOT `constraint=false`).
    - label decision branches AND meaningful edges; an unlabeled fork is ambiguous.
    - group related steps with a dashed box (TikZ `fit` / DOT `cluster`) + a group label.

  STILL THE LAB RULES (paper_figure_quality_rules / figure_design_principles):
    - NO in-figure title (-> the LaTeX caption).
    - light Okabe-Ito / greyscale fills, black outlines; colour carries MEANING
      (one hue per subsystem), not decoration.
    - 10 pt page text; spare -- erase ink that is not a node, an edge, or a label.
""",
        "external_resources": r"""[external_resources] -- where to learn diagram-making, curated 2026-06.

  TikZ
    - "TikZ & PGF" manual (pgf-tikz.github.io / CTAN) -- libraries shapes.geometric,
      arrows.meta, positioning, chains, fit, backgrounds.
    - TeXample.net (texample.net/tikz/examples) -- a large gallery of TikZ diagrams.
    - Overleaf "Creating Flowcharts" tutorial (overleaf.com/learn) -- the node-style idiom.
  Graphviz
    - graphviz.org -- the DOT language reference, attribute list, and gallery; the `dot`
      hierarchical engine for flowcharts.
  Python
    - schemdraw (schemdraw.readthedocs.io) -- flowcharts + circuits from Python.
  Markdown / web
    - Mermaid (mermaid.js.org) -- ```mermaid flowcharts in READMEs (docs only, not paper).
  Standards / principles
    - ISO 5807:1985 -- flowchart symbol semantics.
    - Rougier, "Scientific Visualization: Python + Matplotlib" (2021) -- the layout +
      figure-anatomy chapters (see figure_design_principles('external_resources')).
""",
    }
    q = (topic or "all").strip().lower()
    if q == "all":
        return "\n\n".join(recipes.values()) + r"""
See also: figure_tikz_recipe (general TikZ schematic + PGFPlots + externalize),
figure_design_principles (the 作図 design canon), paper_figure_quality_rules (the gates).
"""
    if q in recipes:
        return recipes[q]
    return (
        f"Unknown topic {topic!r}. Available: "
        f"{', '.join(['all'] + list(recipes))}"
    )


@mcp.tool()
def figure_audit_embeds(tex_path: str) -> str:
    """Lint every \\includegraphics in a LaTeX file for figure embeds that
    cannot guarantee on-page 10 pt @ 8 cm (the CEFC-2026 mistake class).

    Flags:
      * HEIGHT-constrained embeds (``\\includegraphics[height=Xcm]``) -- the
        on-page font size is then uncontrolled.
      * ``width=\\linewidth`` embeds -- the on-page font depends on the column
        width; only safe if the figure was authored AT that exact width.
      * figures not found next to the .tex.
      * figure PDFs that embed DejaVu instead of Times New Roman.

    The compliant pattern: author with ``radia_mcp.figure.lab_figure(
    embed_width_cm=W)``, save with ``save_lab_figure`` (fail-loud gates), and
    embed with the returned ``\\includegraphics[width=W cm]`` at 100%.

    Args:
        tex_path: path to the .tex file to audit.

    Returns a multi-line report (figure count, flagged count, per-figure risks).
    """
    import os
    if not os.path.isfile(tex_path):
        return f"File not found: {tex_path}"
    from ._lab_api import audit_tex_figures
    rep = audit_tex_figures(tex_path)
    lines = [f"figure-embed audit: {rep['tex']}",
             f"  {rep['n_figures']} figures, {rep['n_flagged']} flagged", ""]
    for r in rep["figures"]:
        if r["risks"]:
            lines.append(f"  [FLAG] {r['figure']}  (opts: {r['options'] or '-'})")
            for rk in r["risks"]:
                lines.append(f"         - {rk}")
        else:
            lines.append(f"  [ ok ] {r['figure']}")
    if rep["n_flagged"] == 0:
        lines.append("\nAll embeds use a fixed-cm width and a TNR figure -- clean.")
    else:
        lines.append("\nFix: author at the embed width with "
                     "lab_figure(embed_width_cm=W) + save_lab_figure, then "
                     "\\includegraphics[width=W cm] at 100%.")
    return "\n".join(lines)


# ============================================================
# Self-introspection (uniform with other radia_mcp servers)
# ============================================================

register_status_tool(
    mcp,
    server_name="mcp-server-figure",
    description="Sugahara Lab publication-figure style guide: "
                "IEEE / IEEJ font/size profiles, MATLAB + Matplotlib "
                "snippets, lab style rules (units in parentheses, no "
                "in-figure title, Times New Roman serif).",
    subpackage="radia_mcp.figure",
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
        print("figure MCP server self-test:")
        print(f"  Registered tools ({len(_REGISTERED)}):")
        for name in sorted(_REGISTERED):
            print(f"    - {name}")
        # Smoke-test both legacy tools across every profile so a typo
        # in the _PROFILES dict immediately surfaces in CI.
        profiles = sorted(_tools._PROFILES.keys())
        print(f"  Lab-style profiles ({len(profiles)}):")
        for prof in profiles:
            guide = _tools.figure_style_guide(prof)
            size = _tools.figure_size_for_target(prof)
            assert len(guide) > 200, (
                f"figure_style_guide({prof!r}) returned only "
                f"{len(guide)} chars (suspiciously short)"
            )
            assert "Unknown target" not in guide, (
                f"figure_style_guide({prof!r}) reported unknown target"
            )
            assert "Unknown target" not in size, (
                f"figure_size_for_target({prof!r}) reported unknown target"
            )
            print(f"    {prof:42s} guide={len(guide):4d} ch, "
                  f"size={len(size):4d} ch")
        # Sanity: 'all' returns everything
        full = _tools.figure_style_guide("all")
        assert all(p in full for p in profiles), \
            "figure_style_guide('all') missing some profiles"
        # Sanity: unknown target returns help text (not crash)
        unk = _tools.figure_style_guide("not-a-real-profile")
        assert "Unknown target" in unk and "Valid:" in unk
        print(f"  figure_style_guide('all')               -> "
              f"{len(full):5d} chars")
        print(f"  figure_style_guide('not-a-real-profile') -> "
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

        # Design-principles tool (the 作図 canon, distilled from external refs)
        principles = figure_design_principles("all")
        for sec in ("ten_rules", "perception", "color", "chartjunk",
                    "direct_labeling", "defaults", "external_resources",
                    "sakuzu_vs_graph"):
            assert f"[{sec}]" in principles, f"design topic {sec!r} missing"
        assert ("Rougier" in principles and "Cleveland" in principles
                and "Crameri" in principles), "design-principles citations missing"
        assert "Unknown topic" in figure_design_principles("nope"), \
            "design-principles unknown-topic help missing"
        print(f"  figure_design_principles('all')   -> "
              f"{len(principles):5d} chars")

        # Diagram-recipes tool (flowcharts + concept diagrams; TikZ + Graphviz)
        diagrams = figure_diagram_recipes("all")
        for sec in ("tool_selection", "tikz_flowchart", "graphviz",
                    "concept_diagram", "design", "external_resources"):
            assert f"[{sec}]" in diagrams, f"diagram topic {sec!r} missing"
        assert ("digraph" in diagrams and "tikzpicture" in diagrams
                and "ISO 5807" in diagrams), "diagram-recipes content missing"
        assert "Unknown topic" in figure_diagram_recipes("nope"), \
            "diagram-recipes unknown-topic help missing"
        print(f"  figure_diagram_recipes('all')     -> "
              f"{len(diagrams):5d} chars")

        # Try importing matplotlib + smoke-test the runtime helpers
        # (paper_figure / measure_figure_efficiency).  Skip if mpl is
        # absent (the MCP tools above already loaded without mpl).
        try:
            import matplotlib
            matplotlib.use("Agg")
            from radia_mcp.figure import (
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
