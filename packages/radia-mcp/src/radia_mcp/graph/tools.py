"""Tool functions for radia_mcp.graph (publication-figure style guide).

Promoted on 2026-05-26 from
  s:/mcp-server/src/mcp_server_document/graph/tools.py
into radia-mcp.  The original docstring `mcp_server_document.graph.tools`
import paths in examples have been updated to `radia_mcp.graph` below.

Lab style standards (Sugahara Lab, Kindai University):

  General
    - **No title inside the figure**.  Place the description in the
      LaTeX caption.
    - Times New Roman (serif, matches IEEE/Elsevier conventions).
    - Black axes, thin gridlines, white background.
    - Use parula / viridis colormap for sequential families;
      distinct named colors for categorical comparisons.

  Sizing rule (the key one)
    - Choose the **physical width at which the figure will be embedded**
      in the document, then back-out the font sizes so that text
      appears the right point size **at the embedded width**.
    - Double-column digest, side-by-side panels, embedded width 8 cm:
        font 10 pt, legend 8 pt minimum.
    - Double-column digest, full-width spanning, embedded width 16 cm:
        font 20 pt (so it remains readable when the page is printed at
        100 % or shrunk to single-column reprint).
    - Single-column journal paper, embedded width 8.4 cm: font 8-9 pt.

  Output
    - Save at >= 300 DPI (use 600 DPI for line-art with thin strokes).
    - Prefer PNG with alpha for raster, PDF for vector.
"""
from __future__ import annotations


# ----- Style profile registry -----------------------------------------------

_PROFILES = {
    "digest_double_column_side_by_side": {
        "embed_width_cm": 8.0,
        "font_pt": 10,
        "legend_pt_min": 8,
        "legend_pt_max": 9,
        "linewidth_pt": 1.0,
        "axes_linewidth_pt": 0.8,
        "marker_size_pt": 6,
        "notes": (
            "TWO panels placed SIDE-BY-SIDE (1 row x 2 cols, horizontal) "
            "inside a single 8 cm-wide figure. RULE: each sub-panel must "
            "be <= 4 cm wide (8 cm / 2 = 4 cm; account for the inter-panel "
            "gap so each is slightly under 4 cm, ~3.5-4 cm). Body/axis "
            "font stays 10 pt (the lab absolute font rule -- do NOT shrink "
            "it for the narrow panel). Legend is REDUCED to 8-9 pt (not "
            "10 pt) because a 10 pt legend eats too much of a 4 cm panel; "
            "8-9 pt keeps it readable without crowding. The legend MUST "
            "NOT overlap the graph (curves/markers) -- use "
            "find_best_legend_loc() or label_curve_endpoints(), and "
            "emit_paper_figure() will reject an overlapping legend. Keep "
            "ticks sparse (~4-5 per axis) so the 4 cm panel doesn't crowd."
        ),
    },
    "digest_single_column": {
        "embed_width_cm": 8.0,
        "font_pt": 9,
        "legend_pt_min": 8,
        "linewidth_pt": 1.0,
        "axes_linewidth_pt": 0.8,
        "marker_size_pt": 5,
        "inline_label_pt": 7.5,    # direct-line / endpoint labels
        "annotation_pt": 7.5,      # f_N, slope-guide labels, arrows
        "notes": (
            "Single panel filling one column of a 2-column digest, "
            "embedded at 8 cm via \\begin{figure}.  Use when the figure "
            "has ONE axes (not subplots) -- the panel gets the full "
            "column width so axis-label / tick fonts can match IEEE "
            "single-column convention (9 pt).  Inline annotations "
            "(curve labels at line endpoints, slope-guide labels, "
            "f_N text) target ~7.5 pt -- larger than the 6 pt typical "
            "for the side-by-side 2-panel layout, since each axes is "
            "now twice as wide and supports bigger annotations.  "
            "Aspect ratio 0.45-0.55 typical (vertical height fits "
            "1-page digest with text wrapping around)."
        ),
    },
    "digest_double_column_full_width": {
        "embed_width_cm": 16.0,
        "font_pt": 20,
        "legend_pt_min": 16,
        "linewidth_pt": 1.6,
        "axes_linewidth_pt": 1.2,
        "marker_size_pt": 10,
        "notes": (
            "Single graph spanning both columns of the digest, embedded "
            "at 16 cm. Use figure* in LaTeX. Fonts inflated so the "
            "figure is still readable when reduced to half page in a "
            "single-column reprint."
        ),
    },
    "matlab_oversized_for_8cm_embed": {
        "embed_width_cm": 8.5,
        "font_pt": 30,
        "legend_pt_min": 24,
        "linewidth_pt": 2.0,
        "axes_linewidth_pt": 1.5,
        "marker_size_pt": 12,
        "notes": (
            "MATLAB-specific work-around: render at 24 cm x 9 cm with "
            "30 pt fonts (MATLAB layout breaks at small physical sizes), "
            "then LaTeX scales down by 8.5/24 to embed at single-column "
            "8.5 cm width. Visible font ~ 10.6 pt, legend ~ 8.5 pt."
        ),
    },
    "paper_single_column": {
        "embed_width_cm": 8.89,   # IEEE 3.5 in = 88.9 mm exactly
        "font_pt": 9,
        "legend_pt_min": 8,
        "linewidth_pt": 1.0,
        "axes_linewidth_pt": 0.8,
        "marker_size_pt": 5,
        "notes": (
            "Single column of an IEEE / IEEJ Transactions page: "
            "3.5 in = 88.9 mm = 21 picas (canonical IEEE value). "
            "Use with \\begin{figure} in IEEEtran."
        ),
    },
    "paper_double_column": {
        "embed_width_cm": 18.2,   # IEEE 7.16 in = 182 mm exactly
        "font_pt": 16,
        "legend_pt_min": 14,
        "linewidth_pt": 1.4,
        "axes_linewidth_pt": 1.0,
        "marker_size_pt": 8,
        "notes": (
            "Full-page width spanning both columns of IEEE / IEEJ "
            "Transactions: 7.16 in = 182 mm = 43 picas (canonical "
            "IEEE value). Use with \\begin{figure*} in IEEEtran."
        ),
    },
    "presentation_slide": {
        "embed_width_cm": 25.0,
        "font_pt": 22,
        "legend_pt_min": 18,
        "linewidth_pt": 2.0,
        "axes_linewidth_pt": 1.5,
        "marker_size_pt": 12,
        "notes": (
            "Slide deck (PowerPoint / Beamer). Big fonts so the back row "
            "can read; thicker line widths so contrast survives projector."
        ),
    },
}


_GENERAL_RULES = """\
Lab style — graph rules (Sugahara Lab):

1.  No title inside the figure. Description goes in the LaTeX caption only.
2.  Times New Roman font for both axis labels and legend.
3.  Italic math symbols, upright units: '{\\itB} (T)', '{\\itH} (kA/m)'.
4.  Sparse ticks (~5 per axis). Subdued grid (0.7 alpha gray).
5.  Legend in upper-left or lower-right (whichever is empty), 'Box','off'.
6.  Choose font size for the EMBEDDED width, not the rendered pixel size.
7.  Save at >= 300 DPI PNG; use vector PDF when math symbols are heavy.
8.  Box ON (all four spines), ticks INWARD (xtick.direction='in').
9.  *** UNITS IN PARENTHESES ***, NOT square brackets:
    correct: 'f (Hz)', 'B (T)', 'Temperature (K)'
    wrong:   'f [Hz]', 'B [T]', 'Temperature [K]'
    Per IEEE Editorial Style Manual ('Magnetization (kA/m)',
    'Temperature (K)') and IEEJ practice.  ISO 80000 'f / Hz' is a
    distinct convention NOT used in IEEE/IEEJ engineering papers.
10. *** TWO FIGURES IN 8 cm -> SIDE BY SIDE, each <= 4 cm ***:
    When you must place 2 graphs within an 8 cm width, lay them out
    HORIZONTALLY (1 row x 2 columns), each sub-panel <= 4 cm wide
    (8 cm / 2; leave a small inter-panel gap so each is ~3.5-4 cm).
    - Body / axis-label font: 10 pt (the absolute lab rule -- do NOT
      shrink it just because the panel is narrow).
    - Legend font: 8-9 pt (REDUCED from 10 pt; a 10 pt legend crowds a
      4 cm panel).  This is the ONE sanctioned exception to the
      "legend matches body" rule.
    - Legend MUST NOT overlap the graph -- place it in an empty corner
      ('Box','off') or use direct endpoint labels; verify with
      find_best_legend_loc() / emit_paper_figure()'s overlap gate.
    Use profile 'digest_double_column_side_by_side'.
"""


# ----- Tool functions -------------------------------------------------------

def graph_style_guide(target: str = "all") -> str:
    """Return the lab-standard graph style guide.

    Args:
        target: One of:
            "all"                                    full guide
            "digest_double_column_side_by_side"      side-by-side at 8 cm
            "digest_double_column_full_width"        spanning at 16 cm
            "paper_single_column"                    journal single column
            "paper_double_column"                    journal double column
            "presentation_slide"                     slide deck

    Returns:
        Style rules + numerical settings for the chosen target.
    """
    target = target.strip().lower()

    if target == "all":
        out = [_GENERAL_RULES, "", "Available profiles:", ""]
        for name, prof in _PROFILES.items():
            out.append(f"  [{name}]")
            for k, v in prof.items():
                out.append(f"    {k:20s} = {v}")
            out.append("")
        return "\n".join(out)

    if target not in _PROFILES:
        return (f"Unknown target '{target}'. "
                f"Valid: {', '.join(_PROFILES.keys())}, all")

    prof = _PROFILES[target]
    lines = [_GENERAL_RULES, "", f"Profile: {target}", ""]
    for k, v in prof.items():
        lines.append(f"  {k:20s} = {v}")
    return "\n".join(lines)


def graph_size_for_target(target: str = "digest_double_column_side_by_side",
                          embed_width_cm: float | None = None) -> str:
    """Recommend output figure size + font settings for a target embedding.

    The rule: when the figure is finally embedded at width
    ``embed_width_cm`` in the document, the rendered font must look
    right (e.g. 10 pt at 8 cm).  This tool tells you the size to set
    in the plotting library and the font sizes to use, given the
    embed width.

    Args:
        target: Style profile name (see graph_style_guide).
        embed_width_cm: Override the profile's embed width.

    Returns:
        Recommended figure dimensions + font sizes + Matplotlib /
        MATLAB snippet to apply them.
    """
    target = target.strip().lower()
    if target not in _PROFILES:
        return (f"Unknown target '{target}'. "
                f"Valid: {', '.join(_PROFILES.keys())}")

    prof = _PROFILES[target]
    w_cm = float(embed_width_cm) if embed_width_cm else prof["embed_width_cm"]
    h_cm = w_cm * 0.4   # 5:2 aspect for 2-panel side-by-side; tweak as needed

    font = prof["font_pt"]
    legend = prof["legend_pt_min"]
    legend_max = prof.get("legend_pt_max")
    legend_str = (f"{legend}-{legend_max}" if legend_max else f">= {legend}")
    lw = prof["linewidth_pt"]
    aw = prof["axes_linewidth_pt"]
    ms = prof["marker_size_pt"]

    # Render-units conversion (point = 1/72 inch)
    inches = w_cm / 2.54
    pixels_300dpi = int(round(inches * 300))

    matlab_snippet = (
        "fig = figure('Units','centimeters','Position',[1 1 "
        f"{w_cm:.1f} {h_cm:.1f}], "
        "'PaperUnits','centimeters','PaperSize',["
        f"{w_cm:.1f} {h_cm:.1f}], 'Color','w');\n"
        "set(fig,'DefaultAxesFontName','Times New Roman');\n"
        f"set(fig,'DefaultAxesFontSize',{font});\n"
        "set(fig,'DefaultTextFontName','Times New Roman');\n"
        f"set(fig,'DefaultLineLineWidth',{lw});\n"
        f"set(fig,'DefaultAxesLineWidth',{aw});\n"
        "% legend(...,'FontSize'," f"{legend}" ",'Box','off');\n"
        "% exportgraphics(fig, 'out.png', 'Resolution', 300);"
    )

    matplotlib_snippet = (
        "import matplotlib.pyplot as plt\n"
        f"fig, ax = plt.subplots(figsize=({inches:.2f}, {h_cm/2.54:.2f}), dpi=300)\n"
        "plt.rcParams.update({\n"
        "    'font.family': 'Times New Roman',\n"
        f"    'font.size': {font},\n"
        f"    'axes.linewidth': {aw},\n"
        f"    'lines.linewidth': {lw},\n"
        f"    'lines.markersize': {ms},\n"
        f"    'legend.fontsize': {legend},\n"
        "})\n"
        "# ax.set_title('')   # NO title in figure\n"
        "# fig.savefig('out.png', dpi=300, bbox_inches='tight')"
    )

    return (
        f"Target: {target}\n"
        f"  embed width      : {w_cm:.1f} cm\n"
        f"  recommended size : {w_cm:.1f} cm x {h_cm:.1f} cm "
        f"({pixels_300dpi} px wide @ 300 DPI)\n"
        f"  font / legend    : {font} pt / {legend_str} pt\n"
        f"  line / axis      : {lw} / {aw} pt\n"
        f"  marker           : {ms} pt\n"
        f"  notes            : {prof['notes']}\n"
        "\n"
        "MATLAB snippet:\n" + matlab_snippet + "\n"
        "\n"
        "Matplotlib snippet:\n" + matplotlib_snippet
    )


# ----- Python callable helpers (NEW) ---------------------------------------

def lab_figsize(target: str = "digest_double_column_side_by_side",
                embed_width_cm: float | None = None,
                aspect: float = 0.4) -> tuple[float, float]:
    """Return matplotlib figsize (width_inches, height_inches) for a target.

    Use directly with ``plt.subplots(figsize=lab_figsize(...))``.

    Args:
        target: Profile name (see _PROFILES keys).
        embed_width_cm: Override the profile's embed width.
        aspect: height / width ratio. 0.4 is the default for side-by-side
            2-panel digest figures (5:2). For single-axes use ~0.7-0.75.

    Returns:
        (width_in, height_in) in inches, ready for matplotlib.
    """
    target = target.strip().lower()
    if target not in _PROFILES:
        raise ValueError(
            f"Unknown target '{target}'. "
            f"Valid: {', '.join(_PROFILES.keys())}"
        )
    prof = _PROFILES[target]
    w_cm = float(embed_width_cm) if embed_width_cm else float(prof["embed_width_cm"])
    h_cm = w_cm * float(aspect)
    return (w_cm / 2.54, h_cm / 2.54)


def apply_lab_style(target: str = "digest_double_column_side_by_side",
                    embed_width_cm: float | None = None,
                    aspect: float = 0.4,
                    use_times_roman: bool = True) -> tuple[float, float]:
    """Set matplotlib rcParams to Sugahara-Lab style and return figsize.

    Side effect: mutates ``plt.rcParams`` to match the lab style guide:
    Times New Roman (or sans fallback), tight axes, proper line/marker
    sizes for the chosen embed width.  Idempotent — safe to call multiple
    times within a script.

    Returns:
        (width_in, height_in) for use in ``plt.subplots(figsize=...)``.

    Notes:
        - If Times New Roman is not installed on the system, falls back
          to DejaVu Serif (still serif, IEEE-compatible).  Set
          ``use_times_roman=False`` to skip the font setting altogether
          (use whatever rcParams default).
        - Caller should still:
          * call ``ax.set_title('')`` to remove any default title (lab
            style forbids in-figure titles; description goes in LaTeX
            caption only).
          * call ``fig.savefig(..., dpi=300, bbox_inches='tight')`` for
            production output.

    Example:
        >>> import matplotlib.pyplot as plt
        >>> from radia_mcp.graph import apply_lab_style
        >>> figsize = apply_lab_style("digest_double_column_side_by_side")
        >>> fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        >>> # ... plot ...
        >>> fig.savefig("out.pdf", dpi=300, bbox_inches="tight")
    """
    import matplotlib.pyplot as plt
    import matplotlib

    target = target.strip().lower()
    if target not in _PROFILES:
        raise ValueError(
            f"Unknown target '{target}'. "
            f"Valid: {', '.join(_PROFILES.keys())}"
        )
    prof = _PROFILES[target]

    rc = {
        "font.size": prof["font_pt"],
        "axes.labelsize": prof["font_pt"],
        "axes.titlesize": prof["font_pt"],
        "xtick.labelsize": prof["font_pt"] - 1,
        "ytick.labelsize": prof["font_pt"] - 1,
        "legend.fontsize": prof["legend_pt_min"],
        "axes.linewidth": prof["axes_linewidth_pt"],
        "lines.linewidth": prof["linewidth_pt"],
        "lines.markersize": prof["marker_size_pt"],
        "axes.grid": True,
        "grid.linestyle": ":",
        "grid.alpha": 0.4,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        # --- lab style: box ON, ticks INWARD --------------------------------
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
        # --- legend: no frame (lab style: 'Box','off') ---------------------
        "legend.frameon": False,
        # --- font embedding: TrueType (Type 42) instead of Type 3 ----------
        # Default matplotlib backend produces Type 3 raster glyphs which are
        # CRITICAL reject reason at IEEE / Elsevier (blur at scale).  Force
        # TrueType (Type 42) for both pdf and ps backends so embedded fonts
        # remain vector and editor-friendly.
        # Ref: matplotlib FAQ "Producing PDF/EPS with embedded fonts".
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
    if use_times_roman:
        # Try Times New Roman; fall back gracefully if unavailable.
        # ALSO unify mathtext fontset so $...$ math renders in matching serif.
        available = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
        if "Times New Roman" in available:
            rc["font.family"] = "serif"
            rc["font.serif"] = ["Times New Roman"] + plt.rcParams.get(
                "font.serif", [])
            rc["mathtext.fontset"] = "stix"   # matches Times metrics
        elif "Times" in available:
            rc["font.family"] = "serif"
            rc["font.serif"] = ["Times"] + plt.rcParams.get("font.serif", [])
            rc["mathtext.fontset"] = "stix"
        else:
            rc["font.family"] = "serif"
            rc["font.serif"] = ["DejaVu Serif"]
            rc["mathtext.fontset"] = "dejavuserif"
        # Italic math variables (lab style: '{\\itB} (T)', '{\\itH} (kA/m)')
        rc["mathtext.default"] = "it"

    plt.rcParams.update(rc)
    return lab_figsize(target, embed_width_cm=embed_width_cm, aspect=aspect)


def lab_savefig(fig, path: str, dpi: int = 300, **kwargs) -> None:
    """Save a figure using lab-style defaults: dpi=300, bbox_inches='tight'.

    Convenience wrapper around ``fig.savefig(...)``.  Saves both .pdf and
    .png if path ends in either extension; saves only the requested
    extension otherwise.

    Args:
        fig: matplotlib Figure.
        path: Output filename (with or without extension).
        dpi: DPI for raster output (default 300).
        **kwargs: Forwarded to fig.savefig.

    Example:
        >>> lab_savefig(fig, "circle_warburg")     # saves both .pdf and .png
        >>> lab_savefig(fig, "out.png", dpi=600)   # saves only the .png
    """
    from pathlib import Path
    p = Path(path)
    common = {"dpi": dpi, "bbox_inches": "tight", **kwargs}
    if p.suffix == "":
        fig.savefig(p.with_suffix(".pdf"), **common)
        fig.savefig(p.with_suffix(".png"), **common)
    else:
        fig.savefig(p, **common)


_TIGHTEN_PRESETS = {
    # IGTE'26 digest (2 panels side-by-side, 8 cm embed, no inset legend).
    # Verified to give ~20% more axes width than matplotlib's default
    # tight_layout alone on the 'digest_double_column_side_by_side' profile.
    "digest_2panel_tight": dict(
        pad=0.15, w_pad=0.05, h_pad=0.0,
        left=0.075, right=0.995, top=0.99, bottom=0.20,
    ),
    # IEEE single column (one panel, 8.89 cm wide).
    "ieee_single_column": dict(
        pad=0.20, w_pad=0.0, h_pad=0.0,
        left=0.18, right=0.97, top=0.97, bottom=0.20,
    ),
    # IEEE double column (one panel, 18.2 cm wide, figure*).
    "ieee_double_column": dict(
        pad=0.20, w_pad=0.05, h_pad=0.0,
        left=0.09, right=0.98, top=0.97, bottom=0.16,
    ),
    # Permissive default: matplotlib tight_layout, then a small margin
    # tighten.  Use when you have not measured the figure yet.
    "default": dict(
        pad=0.2, w_pad=0.2, h_pad=0.2,
        left=0.10, right=0.97, top=0.97, bottom=0.18,
    ),
}


def tighten_margins(fig, preset: str = "default", *,
                    pad: float | None = None,
                    w_pad: float | None = None,
                    h_pad: float | None = None,
                    left: float | None = None,
                    right: float | None = None,
                    top: float | None = None,
                    bottom: float | None = None) -> dict:
    """Reduce figure whitespace to lab-tested tight margins.

    Combines :func:`tight_layout` (auto-fit labels/ticks) with explicit
    :func:`subplots_adjust` (known-good tight outer margins).  The
    explicit subplots_adjust step is necessary because matplotlib's
    tight_layout leaves a conservative outer margin; on a 2-panel
    IEEE-2-column digest this margin alone costs ~20% of axes width.

    Picking the right preset is the easy way to get tight margins
    without measuring the figure yourself:

    - ``"digest_2panel_tight"``: IGTE'26 digest 2 panels side-by-side,
      8 cm embed, no inset legend.  Most aggressive — only safe when
      curves are labelled directly at their L/R endpoints (no legend
      reserved inside the axes).
    - ``"ieee_single_column"``: IEEE single column (8.89 cm), one panel.
    - ``"ieee_double_column"``: IEEE double column (18.2 cm) via figure*.
    - ``"default"``: conservative, leaves ~3% outer margin all around.

    Individual values can be overridden via keyword arguments — pass
    ``None`` (default) to use the preset's value, or a numeric value
    to override.

    Args:
        fig: matplotlib Figure.
        preset: One of ``_TIGHTEN_PRESETS`` keys.  Default ``"default"``.
        pad / w_pad / h_pad: tight_layout padding overrides.
        left / right / top / bottom: subplots_adjust overrides
            (figure-fraction, 0..1).

    Returns:
        dict of the actually-applied {pad, w_pad, h_pad, left, right,
        top, bottom} values — useful for logging or for caching the
        layout to reuse across figures.

    Example:
        >>> tighten_margins(fig, "digest_2panel_tight")
        >>> # Override only the left margin:
        >>> tighten_margins(fig, "digest_2panel_tight", left=0.09)
    """
    if preset not in _TIGHTEN_PRESETS:
        raise ValueError(
            f"unknown preset {preset!r}.  Available: "
            f"{sorted(_TIGHTEN_PRESETS)}"
        )
    p = dict(_TIGHTEN_PRESETS[preset])
    # Apply explicit overrides (only non-None)
    for k, v in dict(pad=pad, w_pad=w_pad, h_pad=h_pad,
                     left=left, right=right, top=top, bottom=bottom).items():
        if v is not None:
            p[k] = v

    fig.tight_layout(pad=p["pad"], w_pad=p["w_pad"], h_pad=p["h_pad"])
    fig.subplots_adjust(left=p["left"], right=p["right"],
                        top=p["top"], bottom=p["bottom"])
    return p


def label_curve_endpoints(ax, labels, *,
                          side: str = "right",
                          offset_axes_frac: float = 0.015,
                          fontsize: float = 7.5,
                          va: str = "center",
                          clip_on: bool = False):
    """Place curve labels JUST OUTSIDE the axes — no overlap with data.

    This is a "right-margin direct-labeling" pattern: rather than place
    line labels at the curve endpoints INSIDE the axes (where the
    label's bbox unavoidably covers the curve underneath), each label
    is anchored at ``axes_x = 1 + offset_axes_frac`` (just past the
    right edge) with ``ha="left"`` and ``clip_on=False``.  The label
    text extends into the figure's right margin -- a strict
    no-overlap-with-curves zone.

    The caller must reserve figure-margin space for the labels by
    calling ``fig.subplots_adjust(right=...)`` with a value smaller
    than ``1.0`` (typical 0.75-0.85 depending on longest label).

    Args:
        ax: matplotlib Axes.
        labels: iterable of dicts, each with keys:
            ``y_data`` (float, required): y coordinate in DATA units.
            ``text`` (str, required): label string.
            ``color`` (str): text color (default "k").
            ``fontweight`` (str): default "normal".
            ``fontsize_override`` (float): override default ``fontsize``
                for this one label.
            ``alpha`` (float): default 1.0.
        side: ``"right"`` (default) places labels past the right edge;
            ``"left"`` places past the left edge with ``ha="right"``.
        offset_axes_frac: distance past the axes edge to anchor the
            label (axes-fraction units; default 0.015 = 1.5%).
        fontsize: default fontsize in pt (default 7.5, matching the
            ``digest_single_column`` profile ``inline_label_pt`` value).
        va: vertical alignment ("center", "top", "bottom"; default
            "center").
        clip_on: whether the text gets clipped at the axes boundary
            (default False, so the label can extend into the margin).

    Example:
        >>> # Single-panel digest figure with 4 curve labels
        >>> label_curve_endpoints(ax, [
        ...     {"y_data": 141, "text": "CLN (R-term)", "color": "C1"},
        ...     {"y_data":  27, "text": r"$K/\\sqrt{\\omega}$",
        ...      "color": "k", "alpha": 0.7, "fontsize_override": 7.0},
        ...     {"y_data": 8.5, "text": r"Schur-F = $Y$",
        ...      "color": "C3", "fontweight": "bold"},
        ...     {"y_data": 0.55, "text": "CLN (L-term)", "color": "C0"},
        ... ])
        >>> # Leave room for labels in the right margin:
        >>> fig.subplots_adjust(right=0.78)
    """
    if side == "right":
        x_anchor = 1.0 + offset_axes_frac
        ha = "left"
    elif side == "left":
        x_anchor = -offset_axes_frac
        ha = "right"
    else:
        raise ValueError(f"side must be 'right' or 'left', got {side!r}")

    yax = ax.get_yaxis_transform()
    for spec in labels:
        y = spec["y_data"]
        text = spec["text"]
        color = spec.get("color", "k")
        weight = spec.get("fontweight", "normal")
        fs = spec.get("fontsize_override", fontsize)
        alpha = spec.get("alpha", 1.0)
        ax.text(x_anchor, y, text,
                transform=yax,
                ha=ha, va=va,
                color=color, fontsize=fs, fontweight=weight,
                alpha=alpha, clip_on=clip_on)


def add_slope_guide(ax, x_center, y_center, slope: float, *,
                    dx_decades: float = 0.3,
                    color: str = "gray",
                    lw: float = 0.6,
                    label: str | None = None,
                    label_offset: tuple[float, float] = (0.0, -0.15),
                    label_fontsize: float = 5.5):
    """Draw a slope-indicator triangle on a log-log axis.

    On a log-log plot, a curve segment with power-law $y\\propto x^\\alpha$
    appears as a straight line of slope $\\alpha$.  This helper draws a
    small triangle centred on ``(x_center, y_center)`` with horizontal
    base spanning ``dx_decades`` of $\\log_{10}(x)$ and vertical leg of
    height ``slope * dx_decades`` decades of $\\log_{10}(y)$, plus an
    optional text label.

    Args:
        ax: matplotlib log-log Axes (call ``set_xscale("log")``,
            ``set_yscale("log")`` first).
        x_center, y_center: centre of the guide in DATA coordinates.
        slope: power-law exponent the guide depicts (e.g. -0.5 for
            $y\\!\\propto\\!1/\\sqrt{x}$).
        dx_decades: half-width of the guide in decades of x (default 0.3).
        color: stroke colour for the triangle.
        lw: line width.
        label: optional text label (e.g. r"$-1/2$").  None to omit.
        label_offset: (dx_decades, dy_decades) offset of label from
            (x_center, y_center) in log10 units.
        label_fontsize: fontsize of the label (default 5.5 pt).

    Example:
        >>> ax.loglog(f, y_exact)
        >>> add_slope_guide(ax, x_center=1e6, y_center=27, slope=-0.5,
        ...                 label=r"$-1/2$")
    """
    import numpy as np
    # Triangle vertices in log10 space, then back-transform
    lx0 = np.log10(x_center) - dx_decades
    lx1 = np.log10(x_center) + dx_decades
    ly_left  = np.log10(y_center) - slope * dx_decades
    ly_right = np.log10(y_center) + slope * dx_decades
    # Two-segment "L" shape: horizontal base + sloped hypotenuse
    xs = [10**lx0, 10**lx1, 10**lx0, 10**lx0]
    ys = [10**ly_left, 10**ly_right, 10**ly_right, 10**ly_left]
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="butt",
            zorder=2)
    if label is not None:
        lx_lab = np.log10(x_center) + label_offset[0]
        ly_lab = np.log10(y_center) + label_offset[1]
        ax.text(10**lx_lab, 10**ly_lab, label,
                color=color, fontsize=label_fontsize,
                ha="center", va="center")


def check_legend_overlap(ax, n_samples_per_line: int = 200,
                         margin_px: float = 0.0):
    """Detect whether any plotted curve passes under the legend bbox.

    Use this BEFORE saving a publication figure: a legend that visually
    overlaps a curve is the single most common "obvious flaw" reviewers
    spot.  Programmatic detection is reliable; eyeball-check is not.

    Args:
        ax: A matplotlib Axes containing a legend (created via
            ``ax.legend(...)`` or ``fig.legend(...)``).
        n_samples_per_line: For each Line2D in the axes, sample this many
            equally spaced points along the line and test each against
            the legend bbox.  200 catches most overlaps cheaply; bump to
            500 for very wiggly curves.
        margin_px: Extra padding around the legend bbox in pixels; treat
            curve points inside legend_bbox EXPANDED by margin_px as
            overlapping.  Use 2-4 px for "almost touching" detection.

    Returns:
        A list of dicts, one per overlapping line:
            [{"label": <line label>,
              "n_overlap_points": <count>,
              "fraction": <fraction of sampled points inside bbox>}, ...]
        Empty list means no overlap detected.

    Example:
        >>> ax.plot(x, y1, label="data 1")
        >>> ax.plot(x, y2, label="data 2")
        >>> leg = ax.legend(loc="upper right")
        >>> fig.canvas.draw()      # required before bbox is valid
        >>> from radia_mcp.graph import check_legend_overlap
        >>> overlaps = check_legend_overlap(ax)
        >>> if overlaps:
        ...     for o in overlaps:
        ...         print(f"!! legend overlaps {o['label']} "
        ...               f"({o['fraction']*100:.0f}% of curve)")
        ...     # try ax.legend(loc='lower left') and re-check
    """
    import numpy as np
    from matplotlib.transforms import Bbox

    fig = ax.get_figure()
    # The legend may live on the Axes or on the Figure (fig.legend).
    legend = ax.get_legend()
    if legend is None:
        # Check figure-level legend.
        legends = [c for c in fig.get_children()
                   if c.__class__.__name__ == "Legend"]
        if not legends:
            raise ValueError(
                "Axes has no legend; call ax.legend(...) or "
                "fig.legend(...) before check_legend_overlap.")
        legend = legends[0]

    # Force a draw so legend bbox is finalized.
    fig.canvas.draw()
    leg_bbox = legend.get_window_extent()  # pixel coordinates
    if margin_px > 0:
        leg_bbox = leg_bbox.expanded(
            (leg_bbox.width + 2 * margin_px) / leg_bbox.width,
            (leg_bbox.height + 2 * margin_px) / leg_bbox.height,
        )

    overlaps = []
    for line in ax.lines:
        if not line.get_visible():
            continue
        label = line.get_label()
        if label.startswith("_"):
            # Lines whose label starts with "_" are excluded from legends
            # by matplotlib convention; skip them too.
            continue
        xdata = np.asarray(line.get_xdata(), dtype=float)
        ydata = np.asarray(line.get_ydata(), dtype=float)
        if xdata.size < 2:
            continue
        # Subsample to n_samples_per_line evenly.
        if xdata.size > n_samples_per_line:
            idx = np.linspace(0, xdata.size - 1, n_samples_per_line).astype(int)
            xdata = xdata[idx]
            ydata = ydata[idx]
        # Transform data coords to display pixels.
        pts_data = np.column_stack([xdata, ydata])
        pts_disp = ax.transData.transform(pts_data)
        inside = (
            (pts_disp[:, 0] >= leg_bbox.x0) & (pts_disp[:, 0] <= leg_bbox.x1) &
            (pts_disp[:, 1] >= leg_bbox.y0) & (pts_disp[:, 1] <= leg_bbox.y1)
        )
        n_overlap = int(inside.sum())
        if n_overlap > 0:
            overlaps.append({
                "label": label,
                "n_overlap_points": n_overlap,
                "fraction": float(n_overlap) / pts_disp.shape[0],
            })
    return overlaps


def find_best_legend_loc(ax, candidates=("upper right", "upper left",
                                          "lower right", "lower left",
                                          "center left", "center right")):
    """Try each candidate legend location and return the one with the
    fewest curve overlaps.

    Args:
        ax: Axes with data already plotted.  Calls ax.legend() internally.
        candidates: List of matplotlib `loc` strings to try.

    Returns:
        (best_loc, overlap_count_dict) where best_loc is the `loc` string
        with the minimum total overlap fraction, and overlap_count_dict
        is the full per-location overlap summary for diagnostics.

    Example:
        >>> # plot data, then:
        >>> best, summary = find_best_legend_loc(ax)
        >>> ax.legend(loc=best, frameon=False)
        >>> print(f"chose {best!r}, full summary: {summary}")

    Notes:
        Run AFTER all plotting calls and BEFORE the final savefig.
        The function mutates ax.legend() temporarily but leaves the
        final call to the caller.
    """
    summary = {}
    for loc in candidates:
        ax.legend(loc=loc, frameon=False)
        overlaps = check_legend_overlap(ax)
        total = sum(o["fraction"] for o in overlaps)
        summary[loc] = {"total_fraction": total, "overlaps": overlaps}
    # Remove the temporary legend so the caller can re-create it cleanly.
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()
    best = min(summary, key=lambda k: summary[k]["total_fraction"])
    return best, summary


# ============================================================
# ROM-paper helpers (Cauer / Multi-K / Schur-F asymptote plots)
# ============================================================

def plot_asymptote_ratio_sweep(ax, f_Hz, Y_curves, K_SIBC, *,
                                target_line: bool = True,
                                xlabel: str = r"$f$ (Hz)",
                                ylabel: str | None = None):
    """Plot r(f) = |Y_R(jw)| * sqrt(omega) / K_SIBC log-log.

    Canonical "C-axis" reviewer figure for eddy-current ROM papers:
    rational ROMs (CLN-N, Multi-K) decay as r ~ omega^{-(N-1/2)} above
    their wall band, while a Schur-F-augmented ROM has r flat at 1
    structurally.

    Args:
        ax: matplotlib Axes.
        f_Hz: 1D array of frequencies in Hz.
        Y_curves: dict mapping curve label -> Y(jw) complex 1D array.
                  Recommend including 'Exact' as the reference curve
                  (drawn black, solid, on top).
        K_SIBC: SIBC asymptote coefficient (``S*sqrt(sigma/mu)``
                in Foster-series convention).
        target_line: draw a horizontal r=1 reference dotted line.
        xlabel, ylabel: axis labels (defaults set automatically).

    Curve styling: 'Exact' gets black solid; first 'Schur'-containing
    label gets red solid; remaining labels get cycled dashed/dotted
    styles in green/magenta/blue.  Override by post-hoc ax.lines
    manipulation if needed.

    Returns:
        list of Line2D objects, in the order they were drawn.
    """
    import numpy as np
    f_Hz = np.asarray(f_Hz, dtype=float)
    omega = 2 * 3.141592653589793 * f_Hz

    # Drawing order: exact last (on top), Schur after exact, others before.
    exact_keys = [k for k in Y_curves if "exact" in k.lower()]
    schur_keys = [k for k in Y_curves if "schur" in k.lower()
                   and k not in exact_keys]
    other_keys = [k for k in Y_curves
                   if k not in exact_keys and k not in schur_keys]
    style_cycle = [
        ("b--", 1.2), ("m:", 1.2), ("g-.", 1.2), ("c--", 1.2),
        ("y:", 1.2),
    ]

    lines = []
    for i, k in enumerate(other_keys):
        Y = np.asarray(Y_curves[k])
        r = np.abs(Y) * np.sqrt(omega) / K_SIBC
        style, lw = style_cycle[i % len(style_cycle)]
        lines.append(ax.loglog(f_Hz, r, style, label=k, linewidth=lw,
                                zorder=1 + i)[0])
    for k in schur_keys:
        Y = np.asarray(Y_curves[k])
        r = np.abs(Y) * np.sqrt(omega) / K_SIBC
        lines.append(ax.loglog(f_Hz, r, "r-", label=k, linewidth=2.0,
                                zorder=8)[0])
    for k in exact_keys:
        Y = np.asarray(Y_curves[k])
        r = np.abs(Y) * np.sqrt(omega) / K_SIBC
        lines.append(ax.loglog(f_Hz, r, "k-", label=k, linewidth=1.6,
                                zorder=9)[0])

    if target_line:
        ax.axhline(1.0, color="gray", linewidth=0.5, linestyle=":",
                    zorder=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel if ylabel is not None
                   else r"$r(f) = |Y_R(j\omega)|\sqrt{\omega}/K_{\rm SIBC}$")
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    return lines


def plot_basis_size_convergence(ax, basis_sizes, errors_by_method, *,
                                  xlabel: str = "Basis size $N$",
                                  ylabel: str = "Wall-band error",
                                  show_ceiling: bool = True):
    """Plot wall-band error vs basis size (log-log) for multiple ROMs.

    Standard reviewer figure showing that rational ROMs hit a
    structural error ceiling that Schur-F bypasses.

    Args:
        ax: matplotlib Axes.
        basis_sizes: 1D array of basis sizes.
        errors_by_method: dict mapping method name -> 1D error array
                          aligned with basis_sizes (NaN for unsampled).
        xlabel, ylabel: axis labels.
        show_ceiling: if True, dotted horizontal line at the minimum
                      error across all rational methods (the
                      structural ceiling).

    Returns:
        list of Line2D objects.
    """
    import numpy as np
    basis_sizes = np.asarray(basis_sizes, dtype=float)
    lines = []
    style_cycle = ["o-", "s--", "^:", "d-.", "v-"]
    for i, (label, err) in enumerate(errors_by_method.items()):
        err = np.asarray(err, dtype=float)
        style = style_cycle[i % len(style_cycle)]
        if "schur" in label.lower():
            lines.append(ax.loglog(basis_sizes, err, "r" + style[1:],
                                    label=label, linewidth=2.0,
                                    markersize=6, zorder=5)[0])
        else:
            lines.append(ax.loglog(basis_sizes, err, style,
                                    label=label, linewidth=1.2,
                                    markersize=5, zorder=2 + i)[0])
    if show_ceiling:
        # Compute ceiling: min error across rational methods at largest basis
        rational_errs = [
            np.asarray(e, dtype=float)
            for k, e in errors_by_method.items() if "schur" not in k.lower()
        ]
        if rational_errs:
            ceiling = min(np.nanmin(e) for e in rational_errs)
            ax.axhline(ceiling, color="gray", linewidth=0.5,
                        linestyle=":", zorder=0,
                        label=f"rational ceiling $\approx {ceiling:.1e}$")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    return lines
