"""Catalog of 22 2D chart specs for radia_mcp.chart2d.

Each ChartSpec describes:
  - the matplotlib call name to use (and any required imports)
  - which kwargs the user provides (data + per-chart-type style knobs)
  - a recipe template (Python code with placeholders) for the
    text-return mode
  - a renderer that converts the user kwargs into a (fig, axes)
    via radia_mcp.graph.paper_figure so the output inherits the
    lab paper-quality style (CVD-safe Okabe-Ito palette, 10 pt @ 8 cm,
    Type-42 font embed, units-in-parentheses).

The catalog is consumed by:
  - server.py to register one @mcp.tool() per chart type
  - _render.py to dispatch type -> renderer
  - _recipe.py to dispatch type -> recipe template

22 chart types, grouped:

  TIME-SERIES / LINE (8):
    line, loglog, semilogx, semilogy, step, errorbar, fill_between, bode
  DISTRIBUTION / STATISTICS (5):
    histogram, bar, box, violin, ecdf
  SCIENTIFIC 2D FIELDS (7):
    contour, contourf, pcolormesh, quiver, streamplot, imshow, polar
  POINT (2):
    scatter, phase

phase = scatter on the complex plane (Re vs Im) -- useful for impedance
loci, root-locus, Nyquist; not strictly a separate matplotlib call but
deserves a dedicated tool so the AI can ask for "Nyquist" / "Bode"
naturally without re-deriving the recipe each time.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional


@dataclass(frozen=True)
class ChartSpec:
    """Describe one chart type.

    Attributes:
        name: short canonical name (matches the MCP tool name suffix).
        mpl_call: matplotlib Axes method invoked by the renderer.
        group: 'line' | 'dist' | 'field' | 'point' -- for catalog grouping.
        description: 1-line human description, surfaced via chart2d_catalog.
        x_required: does the chart need an x array? (some accept y only)
        x_units_hint: example unit string for the x axis (for the
            UI helper text, e.g. '(Hz)' for Bode plots).
        y_units_hint: same for y.
        aliases: tuple of synonyms (e.g. 'nyquist' for 'phase').
        special_kwargs: list of (name, default, doc) tuples beyond
            the universal x/y/labels/xlabel/ylabel/...
        notes: longer guidance shown when the AI asks for the recipe.
    """

    name: str
    mpl_call: str
    group: str
    description: str
    x_required: bool = True
    x_units_hint: str = ""
    y_units_hint: str = ""
    aliases: tuple = ()
    special_kwargs: tuple = ()
    notes: str = ""


# ===========================================================
# 8 LINE / TIME-SERIES
# ===========================================================

LINE = ChartSpec(
    name="line",
    mpl_call="plot",
    group="line",
    description="Standard line plot; default for time-series + 1D sweeps.",
)
LOGLOG = ChartSpec(
    name="loglog",
    mpl_call="loglog",
    group="line",
    description="Log-log plot; both axes log10.",
    notes="Use for power-law data: y = k * x^alpha appears as a straight"
          " line of slope alpha.  Pair with add_slope_guide() from"
          " radia_mcp.graph to annotate slope reference triangles.",
)
SEMILOGX = ChartSpec(
    name="semilogx",
    mpl_call="semilogx",
    group="line",
    description="Linear y, log10 x (Bode-magnitude style).",
    x_units_hint="(Hz)",
    notes="Default for any frequency-domain magnitude or phase plot"
          " when frequency spans > 1 decade.  Bode-gain top half.",
)
SEMILOGY = ChartSpec(
    name="semilogy",
    mpl_call="semilogy",
    group="line",
    description="Linear x, log10 y (exponential decay / growth).",
    notes="Use for decaying / growing quantities (relaxation curves,"
          " exponential noise tails).  Straight line on this scale"
          " => purely exponential.",
)
STEP = ChartSpec(
    name="step",
    mpl_call="step",
    group="line",
    description="Stair-step (constant between samples).",
    special_kwargs=(("where", "post",
                     "'pre'|'post'|'mid' -- which side of each x"
                     " the step lands"),),
    notes="For switching waveforms (PWM), discrete-time signals, or"
          " histograms drawn as steps rather than bars.",
)
ERRORBAR = ChartSpec(
    name="errorbar",
    mpl_call="errorbar",
    group="line",
    description="Line + symmetric or asymmetric error bars.",
    special_kwargs=(
        ("yerr", None, "y-error magnitudes (scalar / array / asymmetric)"),
        ("xerr", None, "x-error magnitudes"),
        ("capsize", 3.0, "cap length at error-bar tip (pt)"),
        ("fmt", "o", "marker / line spec"),
    ),
    notes="Use for measured data with uncertainty.  Pin `fmt='o'`"
          " for points-only (no connecting line) or `fmt='-o'` for"
          " markers + connecting line.",
)
FILL_BETWEEN = ChartSpec(
    name="fill_between",
    mpl_call="fill_between",
    group="line",
    description="Filled band between two y arrays (confidence/hysteresis).",
    special_kwargs=(
        ("y2", None, "lower bound of the band (default 0)"),
        ("alpha", 0.3, "fill opacity"),
    ),
    notes="Use for confidence intervals (y +/- 1*sigma -> y1 = mean-sigma,"
          " y2 = mean+sigma).  Hysteresis loops: pass two separate"
          " fill_between calls for the ascending vs descending branch.",
)
BODE = ChartSpec(
    name="bode",
    mpl_call="bode_pair",   # SYNTHETIC -- handled specially in renderer
    group="line",
    description="Bode pair: magnitude (semilogx) above, phase (semilogx)"
                " below in a 2x1 stack.",
    x_units_hint="(Hz)",
    aliases=("bode_plot",),
    special_kwargs=(
        ("freq", None, "frequency array (Hz), shared x for both panels"),
        ("magnitude", None, "|H(jw)| array (linear or dB)"),
        ("phase", None, "arg(H(jw)) array (degrees)"),
        ("mag_db", True, "True: y label '|H| (dB)'; False: '|H|'"),
        ("phase_unit", "deg",
         "'deg' or 'rad' -- axis label adapts automatically"),
    ),
    notes="Builds a 2x1 paper_figure layout, shares the frequency axis"
          " between the gain and phase subplots, drops the x tick"
          " labels on the top panel (cleaner).",
)


# ===========================================================
# 5 DISTRIBUTION / STATISTICS
# ===========================================================

HISTOGRAM = ChartSpec(
    name="histogram",
    mpl_call="hist",
    group="dist",
    description="1D histogram.",
    x_required=False,           # accepts `data` directly
    special_kwargs=(
        ("data", None, "1D data array (positional or 'data=')"),
        ("bins", 30, "int bin count or 1D edges array"),
        ("density", False,
         "True: density (integral = 1); False: counts"),
    ),
)
BAR = ChartSpec(
    name="bar",
    mpl_call="bar",
    group="dist",
    description="Vertical bar chart; categorical x.",
    special_kwargs=(
        ("categories", None, "x-axis category labels (str list)"),
        ("heights", None, "bar heights"),
        ("width", 0.8, "bar width as fraction of category slot"),
    ),
)
BOX = ChartSpec(
    name="box",
    mpl_call="boxplot",
    group="dist",
    description="Box-and-whisker plot.",
    x_required=False,
    special_kwargs=(
        ("data", None, "list of 1D arrays (one per group) or 2D ndarray"),
        ("labels", None, "category labels for each box"),
        ("showfliers", True, "True: show outlier points beyond whiskers"),
    ),
)
VIOLIN = ChartSpec(
    name="violin",
    mpl_call="violinplot",
    group="dist",
    description="Violin plot (kernel density estimate as fill).",
    x_required=False,
    special_kwargs=(
        ("data", None, "list of 1D arrays or 2D ndarray"),
        ("labels", None, "category labels (set as xticks)"),
        ("showmedians", True, "True: draw median bar"),
    ),
)
ECDF = ChartSpec(
    name="ecdf",
    mpl_call="ecdf",
    group="dist",
    description="Empirical cumulative distribution.",
    x_required=False,
    special_kwargs=(("data", None, "1D data array"),),
    notes="matplotlib 3.8+ has ax.ecdf() natively.  Useful when"
          " comparing distributions where the shape of the tail"
          " matters more than the mode.",
)


# ===========================================================
# 7 SCIENTIFIC 2D FIELDS
# ===========================================================

CONTOUR = ChartSpec(
    name="contour",
    mpl_call="contour",
    group="field",
    description="Contour lines of a 2D scalar field Z(X, Y).",
    special_kwargs=(
        ("X", None, "2D meshgrid array"),
        ("Y", None, "2D meshgrid array"),
        ("Z", None, "2D scalar field"),
        ("levels", None, "int (count) or 1D contour values"),
    ),
)
CONTOURF = ChartSpec(
    name="contourf",
    mpl_call="contourf",
    group="field",
    description="Filled contour (heatmap-with-levels).",
    special_kwargs=(
        ("X", None, "2D meshgrid"),
        ("Y", None, "2D meshgrid"),
        ("Z", None, "2D scalar field"),
        ("levels", None, "level edges"),
        ("cmap", "viridis",
         "perceptually uniform colormap (viridis, magma, plasma, cividis)"),
    ),
    notes="Default cmap viridis is perceptually uniform AND CVD-safe."
          " Avoid `jet` -- non-uniform luminance, confuses readers"
          " (Borland & Taylor 2007).",
)
PCOLORMESH = ChartSpec(
    name="pcolormesh",
    mpl_call="pcolormesh",
    group="field",
    description="2D heatmap on irregular grid (per-cell color).",
    special_kwargs=(
        ("X", None, "1D or 2D x edges/centers"),
        ("Y", None, "1D or 2D y edges/centers"),
        ("Z", None, "2D scalar field"),
        ("cmap", "viridis", "perceptually uniform colormap"),
        ("shading", "auto", "'auto'|'flat'|'gouraud'"),
    ),
)
QUIVER = ChartSpec(
    name="quiver",
    mpl_call="quiver",
    group="field",
    description="2D vector field as arrows.",
    special_kwargs=(
        ("X", None, "x positions (1D or 2D)"),
        ("Y", None, "y positions"),
        ("U", None, "x components of vectors"),
        ("V", None, "y components of vectors"),
        ("scale", None, "arrow scale (None = auto)"),
        ("scale_units", "xy", "'xy'|'width'|'height'|None"),
    ),
)
STREAMPLOT = ChartSpec(
    name="streamplot",
    mpl_call="streamplot",
    group="field",
    description="Field lines of a 2D vector field.",
    special_kwargs=(
        ("X", None, "1D x grid (REGULAR)"),
        ("Y", None, "1D y grid (REGULAR)"),
        ("U", None, "x component (2D)"),
        ("V", None, "y component (2D)"),
        ("density", 1.0, "stream density"),
        ("color", None, "single color or 2D scalar for color-mapping"),
    ),
    notes="streamplot requires a REGULAR (uniform) X/Y grid -- if"
          " your data is irregular, interpolate first via"
          " scipy.interpolate.griddata.",
)
IMSHOW = ChartSpec(
    name="imshow",
    mpl_call="imshow",
    group="field",
    description="Display 2D array as image / heatmap.",
    special_kwargs=(
        ("data", None, "2D array (or 3D for RGB)"),
        ("cmap", "viridis", "perceptually uniform colormap"),
        ("origin", "lower",
         "'lower' for math axes; 'upper' for image-style (row 0 on top)"),
        ("extent", None, "(xmin, xmax, ymin, ymax) for axis range"),
        ("aspect", "auto", "'auto'|'equal'|<float>"),
    ),
)
POLAR = ChartSpec(
    name="polar",
    mpl_call="plot",        # uses projection='polar' axis
    group="field",
    description="Polar / radial plot (r vs theta).",
    special_kwargs=(
        ("theta", None, "angle array (rad)"),
        ("r", None, "radial values"),
    ),
    notes="Pair with theta_zero_location='N' (north) and"
          " theta_direction=-1 (clockwise) for compass-like plots.",
)


# ===========================================================
# 2 POINT
# ===========================================================

SCATTER = ChartSpec(
    name="scatter",
    mpl_call="scatter",
    group="point",
    description="Scatter plot of (x, y) points.",
    special_kwargs=(
        ("s", None, "marker size (scalar or per-point array)"),
        ("c", None,
         "color (single value or per-point array for color-mapping)"),
        ("alpha", 0.7, "marker opacity"),
        ("cmap", "viridis", "colormap when c is an array"),
    ),
)
PHASE = ChartSpec(
    name="phase",
    mpl_call="plot",        # complex Re vs Im
    group="point",
    description="Complex-plane (Re vs Im) plot.  Nyquist / impedance"
                " locus / root locus.",
    aliases=("nyquist", "complex_plane", "argand", "smith_like"),
    special_kwargs=(
        ("z", None,
         "complex-valued array OR (real, imag) tuple of arrays"),
    ),
    notes="Sets equal aspect ratio so circles look circular.  For a"
          " Smith-chart-grade plot use a dedicated package (skrf,"
          " scikit-rf); this tool is for general impedance locus."
          " Pair with chart2d_polar if you want |Z| vs arg(Z).",
)


# ===========================================================
# Master catalog
# ===========================================================

CHARTS: dict[str, ChartSpec] = {
    s.name: s for s in [
        LINE, LOGLOG, SEMILOGX, SEMILOGY, STEP, ERRORBAR,
        FILL_BETWEEN, BODE,
        HISTOGRAM, BAR, BOX, VIOLIN, ECDF,
        CONTOUR, CONTOURF, PCOLORMESH, QUIVER, STREAMPLOT,
        IMSHOW, POLAR,
        SCATTER, PHASE,
    ]
}

# Alias resolution: nyquist -> phase, etc.
ALIASES: dict[str, str] = {}
for _s in CHARTS.values():
    for _a in _s.aliases:
        ALIASES[_a] = _s.name


def get_spec(name: str) -> ChartSpec:
    """Resolve a chart type by name OR alias."""
    key = name.strip().lower()
    if key in CHARTS:
        return CHARTS[key]
    if key in ALIASES:
        return CHARTS[ALIASES[key]]
    raise ValueError(
        f"Unknown chart type {name!r}.  Available: "
        f"{', '.join(sorted(CHARTS))}.  Aliases: "
        f"{', '.join(sorted(ALIASES))}"
    )


def chart_groups() -> dict[str, list[str]]:
    """Return {group_name: [chart_name, ...]} for the catalog tool."""
    g = {}
    for name, spec in CHARTS.items():
        g.setdefault(spec.group, []).append(name)
    return g
