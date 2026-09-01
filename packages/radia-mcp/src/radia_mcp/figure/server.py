"""radia_mcp.figure.server -- figure tool implementations.

Sugahara Lab publication-figure MCP tools.  As of 2026-07-18 figure no
longer runs a standalone MCP server: it was a shared middle-layer server
for paper-writing + presentation, and since presentation merged into
paper-writing (2026-07-17), figure follows.  The mcp-server-figure entry
point is retired; every figure_* / paper_figure_* tool is now served by
mcp-server-paper-writing via radia_mcp.figure.register().

This module only DEFINES the six inline tool functions
(paper_figure_profiles / paper_figure_recipe / paper_figure_quality_rules
/ figure_design_principles / figure_diagram_recipes / figure_audit_embeds);
the figure_* helpers in radia_mcp.figure.tools are registered alongside
them by register().  Nothing runs at import time except the imports below.

The companion Python helpers (apply_lab_style, lab_figsize, lab_savefig,
...) are NOT MCP tools -- they have matplotlib-typed signatures that do
not round-trip through MCP-JSON; they live in radia_mcp.figure.tools for
direct import by analysis scripts.
"""

from . import tools as _tools  # noqa: F401  (kept per module contract)
from . import _paper_figure as _paper


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


def paper_figure_quality_rules(query: str = "all") -> str:
    """Why paper-quality figures need a margin-efficiency gate.

    Returns text on what 'axes-area / total-area fraction' is, why the
    paper-quality floor is around 0.78, what aspects to watch (units in
    parentheses not brackets, no in-figure title, TrueType font embed
    pdf.fonttype=42, 10 pt minimum visible text for IEEE/IEEJ figures), and how
    auto_tighten + emit_paper_figure compose to make the workflow
    refuse to ship a wasteful figure.

    Topics:
        'all'             - full text
        'efficiency'      - how axes_area_fraction is computed + thresholds
        'margins'         - per-margin breakdown reading guide
        'units'           - units in parentheses convention (IEEE/IEEJ)
        'font_family'     - Times New Roman, no fallback
        'font_embedding'  - Type 42 requirement
        'multipanel'      - 1x2 / 2x1 / 2x2 layout tactics
        'side_by_side'    - two figures in 8 cm -> each <= 4 cm,
                            every visible font >= 10 pt, no overlap
        'slide_169'       - author at 24 pt; actual slide display >= 20 pt
        'export_targets'  - format matrix: vector PDF (paper) / EMF
                            (Word/PowerPoint, MATLAB -dmeta) / PNG
                            400 dpi (draft/web) -- which format for
                            which venue, from real lab scripts
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
  For ordinary paper profiles, the save-time gate accepts 10.0--10.5 pt
  on the final page.  This is a target band, not merely a minimum: oversized
  tick labels and legends make the paper look typographically unbalanced.

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
  tick   = 10 pt (the 10 pt floor applies to every visible text item)

This rule pins every paper_figure profile's font_pt to 10.0 regardless
of column width.  If you find yourself wanting smaller text "to make
the axes fit", the answer is to use a wider column or simplify the
plot -- never to shrink font below 10 pt.

There is no small-legend exception. For two panels side-by-side in an
8 cm figure, legends and annotations also remain >=10 pt. If they do
not fit, simplify the content, use direct labels, or allocate more width.
""",
        "slide_169": """\
[slide_169]

SUGAHARA LAB 16:9 SLIDE FIGURE RULE:

  Author at 24 pt where practical; every font must be at least 20 pt
  after scaling to the actual pasted width.

This floor applies to axis labels, tick labels, legends, panel labels, and
free annotations.  Use `beamer_169_full`, `beamer_169_half`, or the
`presentation_slide` style profile.  These profiles set every default to
24 pt. `emit_paper_figure(..., embed_width_cm=<actual pasted width>)`
audits every text artist after scaling, so a 24 pt source that becomes
19 pt in PowerPoint is rejected before save.

If the 20 pt displayed floor does not fit, simplify the plot, reduce the
number of panels, or allocate more slide area.
The no-in-figure-title rule still applies: the title belongs in the slide
title band or Beamer frame title.
""",
        "font_family": """\
[font_family]

LAB RULE: FIGURES ARE GENERATED IN TIMES NEW ROMAN.

Times New Roman is not a preference here; it is the Sugahara Lab
standard for publication figures.  It matches the expected IEEE /
IEEJ / Elsevier paper typography and keeps labels, tick text, legends,
and math-adjacent text visually consistent across figures.

paper_figure(..., use_times_roman=True) is the default and now fails
loudly if matplotlib cannot resolve Times New Roman.  It does NOT
silently fall back to DejaVu Serif, Liberation Serif, or a generic
Times-like family.

Required rcParams (paper_figure() sets these):
    font.family = 'serif'
    font.serif  = ['Times New Roman']
    mathtext.fontset = 'stix'

emit_paper_figure() also checks the save path.  If a figure was built
outside paper_figure(), it requires rcParams to request Times New Roman
and scans the final PDF for common serif fallback fonts.

Use use_times_roman=False only for deliberate non-public scratch plots.
Do not use it for paper, slide, report, or NotebookLM source figures.
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
        "text_overflow": """\
[text_overflow]

LAB RULE: DIAGRAM TEXT MUST FIT INSIDE THE FIXED FIGURE CANVAS.

Concept diagrams often use `ax.text(...)`, `annotate(...)`, or
`fig.text(...)` instead of x/y axis labels.  These free labels can be
cut off at the canvas edge even when ordinary graph checks pass.

emit_paper_figure() defaults to `check_text_overflow=True` and scans
free diagram text before saving.  If a label extends past the figure
canvas, it raises with the side and overhang in points.

Fixes:
  - move the label inside the axes or figure,
  - shorten the label,
  - increase aspect / width,
  - reserve margin with `fig.subplots_adjust(...)`.

Use `audit_text_overflow(fig)` directly while building a diagram to
inspect offenders before calling `emit_paper_figure()`.
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
  For each axis-legend pair, test the rendered Line2D paths, markers,
  and collection offsets (including scatter points) against the legend
  bbox in display coordinates.  This catches both dense curves and a
  sparse two-point segment that crosses the legend.  Any intersection
  = overlap = fail the gate.

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
  - legend / annotation      : >= 10 pt. There is no small-text exception.
                              Simplify or move content if it does not fit.

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
                fontsize=10)
  emit_paper_figure(fig, 'out', 'ieee_double_column', on_fail='raise')

Or use the size/font recipe directly:
  figure_size_for_target('digest_double_column_side_by_side')
  # -> 8 cm wide, every visible font >= 10 pt
""",
        "export_targets": """\
[export_targets]

LAB EXPORT-FORMAT MATRIX (extracted 2026-06 from internal MATLAB
+ matplotlib scripts).  The correct
format depends on WHERE the figure is embedded -- the lab uses three
in practice, not just paper PDF:

  Venue / embed              Format          How
  -------------------------  --------------  --------------------------
  IEEE/IEEJ/IGTE paper       vector PDF      emit_paper_figure(...) ->
   (LaTeX includegraphics)    (Type-42)       .pdf  (this server default)
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

OBSERVED LAB HABIT (from internal scripts):
  - FEMM MATLAB (legacy): set(gca,'FontName','Times'); xlabel('{\\it X}
                   (m)'); print('-dmeta','f.emf');  % print is DEPRECATED
                   -> modern: exportgraphics(...,'ContentType','vector')
  - COMSOL/CoreformCubit matplotlib: figsize=(3,4), dpi=400, Times New
                   Roman 10 pt, inward ticks, savefig PNG.

RULE: NEVER embed a raster PNG in a CAMERA-READY paper -- re-render to
vector PDF.  PNG/EMF are for drafts + Office.  EMF (not PNG)
is the right Office format because it keeps the text vector inside Word.

VERIFICATION (2026-06): the lab's internal script archive AGREES with this
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


def figure_diagram_recipes(topic: str = "all") -> str:
    r"""Flowchart + conceptual/schematic DIAGRAM recipes -- the diagram-DRAWING skill
    (distinct from data-PLOTTING グラフ and from the general 作図 design canon in
    figure_design_principles).

    TikZ and Graphviz were abolished on 2026-09-01, so this tool now carries the
    sanctioned route only: schemdraw for circuits and flowcharts, matplotlib / MATLAB
    for plots, gmsh for fields, meshes and CAD -- plus the diagram DESIGN conventions
    (ISO 5807 symbols, flow direction, crossing minimisation) that outlive any tool.

    Topics:
        'all'                - everything
        'tool_selection'     - the two principles, and which tool for which figure
        'concept_diagram'    - conceptual / block-relationship diagrams (architecture style)
        'design'             - ISO 5807 symbols + flow direction + crossing/alignment rules
        'external_resources' - manuals, galleries, standards
    """
    recipes = {
        "tool_selection": r"""[tool_selection] -- TikZ and Graphviz are ABOLISHED (lab decision 2026-09-01).

  TWO PRINCIPLES.  A tool is sanctioned only if it satisfies both.

  1. DERIVED GEOMETRY.  The geometry comes from the data, from the model, or
     from the relative placement of elements.  NEVER hand-typed coordinates.
     Hand coordinates are how a wire ends up drawn through a box.  This matters
     more when an LLM is drawing: it emits coordinates fluently and cannot see
     the result, so it makes that mistake and does not notice.  schemdraw's
     .right() / .down() / .at() / anchors make the defect impossible instead of
     something a reviewer has to catch.

  2. THE FIGURE IS A FILE.  Every figure is a standalone PDF/PNG written by a
     committed script, independent of the document that includes it.  A figure
     that exists only as source inside a .tex cannot be put on a slide, cannot
     be audited (audit_pptx_figures and check_embedded_figure_text_size both
     read a file), cannot be reused by paper / slide / poster / notebook, and
     cannot even be looked at without compiling the whole document.

  Kind of figure                     Tool                    Geometry from
  ---------------------------------  ----------------------  --------------------
  data plot, scatter, sweep, Bode    matplotlib  or  MATLAB  the data
  circuit diagram                    schemdraw               relative placement
  flowchart / block diagram          schemdraw.flow          relative placement
  field / mesh / post-processing     gmsh                    the model
  geometry, CAD, concept solid       gmsh                    the model
                                     (author the solid with radia_mcp.build123d's
                                      modeling verbs / archetypes, export STEP,
                                      render in gmsh)

  NOT sanctioned: TikZ / PGF and matlab2tikz (hand coordinates, and the figure
  never becomes a file), Graphviz / DOT (same file objection; placement is
  automatic rather than derived from anything the author controls), draw.io
  (not diffable), Mermaid outside a README, and raw matplotlib primitives
  standing in for schemdraw.

  THE SCRIPT LIVES NEXT TO THE FIGURE.  A `make_figs.py` in the same directory,
  committed, that regenerates every figure there and prints the numbers the
  document quotes.  This is the same rule as the lab's data-persistence policy:
  a committed figure whose source is not committed beside it is not reproducible.

  Everything routes through radia_mcp.figure for sizing, saving and auditing:

      fig, axes = paper_figure("beamer_169_full", aspect=0.40)   # or a paper profile
      ax = axes.ravel()[0]
      with schemdraw.Drawing(canvas=ax, show=False) as d:        # circuits/flowcharts
          d.config(fontsize=24, lw=2.2, color=INK)
          d += elm.Inductor().right().label("$L_1$").color(BULK)
      lab_savefig(fig, out, medium="presentation", embed_width_cm=15.0)

  schemdraw validates colour names itself and REJECTS matplotlib's "tab:blue"
  form -- pass hex ("#1f77b4") or a plain CSS name.  Define the palette once and
  share it with the data plots.  One source then gives both the coloured slide
  version and the mono paper version: the palette is the only thing that changes.

  schemdraw gotchas worth knowing before you start:
    - flow.Box hands the cursor to its EAST anchor.  Going downward, say so:
      box.at(arrow.end).anchor("N").
    - two-terminal label locations are 'top'/'bot'/'lft'/'rgt'.  Passing 'right'
      silently centres the caption on whatever is below.  For an arrow caption,
      derive the position from arrow.center and place it with ax.text.
""",
        "concept_diagram": r"""[concept_diagram] -- conceptual / block-relationship diagrams
(architecture, data flow, "X consumes Y"): boxes + LABELLED arrows + optional grouping;
not a strict process flow.  schemdraw.flow, drawn onto a radia_mcp.figure axes.

    import matplotlib.pyplot as plt
    import schemdraw
    from schemdraw import flow
    from radia_mcp.figure import lab_savefig, paper_figure

    fig, axes = paper_figure("ieee_single_column", aspect=0.42)
    ax = axes.ravel()[0]

    with schemdraw.Drawing(canvas=ax, show=False) as d:
        d.config(fontsize=10, lw=1.0)
        cubit = flow.Box(w=2.6, h=1.1).label("Cubit\nhex mesh")
        d += cubit
        d += flow.Arrow().right().length(1.4).label("export", ofst=0.75)
        vol = flow.Box(w=1.8, h=1.1).label(".vol")
        d += vol
        d += flow.Arrow().right().length(1.4).label("Mesh()", ofst=0.75)
        ng = flow.Box(w=2.6, h=1.1).label("NGSolve\nFEM")
        d += ng

    # the group box comes from the members' own anchors, never from typed numbers
    left, right = vol.W[0], ng.E[0]
    top, bot, pad = ng.N[1], ng.S[1], 0.28
    ax.add_patch(plt.Rectangle((left - pad, bot - pad),
                               (right - left) + 2 * pad, (top - bot) + 2 * pad,
                               fill=False, ls="--", lw=1.0, edgecolor="#1f77b4"))
    ax.text((left + right) / 2, bot - pad - 0.28, "computation",
            ha="center", va="top", color="#1f77b4", fontsize=10)

    ax.set_aspect("equal"); ax.axis("off"); fig.tight_layout(pad=0.15)
    lab_savefig(fig, "fig/pipeline", medium="paper", embed_width_cm=8.89)

  This template was run before being written down.  Two things it earns you:
    - ofst=0.75 on the arrow captions.  A horizontal arrow's caption defaults to
      just above the shaft, which is INSIDE tall neighbouring boxes.  Clear h/2.
    - the group box from vol.W / ng.E / ng.N / ng.S.  Typing its corners is the
      failure this whole policy exists to prevent.

  The lab CLAUDE.md ASCII box-diagrams (the 4-Layer panel architecture, the
  accelerator-magnet pipeline) convert directly: one flow.Box per ASCII box, one
  flow.Arrow per `->`.  Keep ONE flow direction.
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
    - align nodes on a grid (schemdraw's unit spacing and anchors do this for you);
      ragged placement reads as careless.
    - MINIMISE edge crossings; route feedback / loop edges orthogonally around the
      side (schemdraw: Line().right().down().left(), or .to() an anchor).
    - label decision branches AND meaningful edges; an unlabeled fork is ambiguous.
    - group related steps with a dashed box + a group label; derive its corners from
      the members' anchors (see 'concept_diagram'), never type them.

  STILL THE LAB RULES (paper_figure_quality_rules / figure_design_principles):
    - NO in-figure title (-> the LaTeX caption).
    - light Okabe-Ito / greyscale fills, black outlines; colour carries MEANING
      (one hue per subsystem), not decoration.
    - 10 pt page text; spare -- erase ink that is not a node, an edge, or a label.
""",
        "external_resources": """\
[external_resources] -- where to learn diagram-making, curated 2026-09.

  Python (the sanctioned route)
    - schemdraw (schemdraw.readthedocs.io) -- circuits and flowcharts from Python;
      relative placement, so overlaps cannot arise from a mistyped coordinate.
    - matplotlib / MATLAB for plots; gmsh for fields, meshes and CAD.
  Markdown / web
    - Mermaid (mermaid.js.org) -- ```mermaid flowcharts in READMEs (docs only,
      never a paper or slide figure: it is not a file).
  Standards / principles
    - ISO 5807:1985 -- flowchart symbol semantics.
    - Rougier, "Scientific Visualization: Python + Matplotlib" (2021) -- the layout +
      figure-anatomy chapters (see figure_design_principles('external_resources')).

  TikZ and Graphviz were abolished on 2026-09-01; their references are deliberately
  not listed here.  See 'tool_selection' for the two principles behind that.
""",
    }
    q = (topic or "all").strip().lower()
    if q == "all":
        return "\n\n".join(recipes.values()) + r"""
See also:
figure_design_principles (the 作図 design canon), paper_figure_quality_rules (the gates).
"""
    if q in recipes:
        return recipes[q]
    return (
        f"Unknown topic {topic!r}. Available: "
        f"{', '.join(['all'] + list(recipes))}"
    )


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
