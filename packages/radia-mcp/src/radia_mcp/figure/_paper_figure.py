"""radia_mcp.figure.paper_figure — paper-grade figure scaffolds + gates.

What this module adds on top of `tools.py`:

  PROFILES (PaperProfile dataclass)
    IEEE_SINGLE_COLUMN / IEEE_DOUBLE_COLUMN      — IEEE Transactions
    IEEJ_SINGLE_COLUMN / IEEJ_DOUBLE_COLUMN      — IEEJ-D / IEEJ-B
    IGTE_DIGEST_DOUBLE                           — IGTE / Compumag digest
    plus the looser "*_loose" variant of each for first-pass drafts

  paper_figure(profile, nrows=1, ncols=1, aspect=None, panel_labels=False)
      -> (fig, axes_2d)
    One-shot scaffold: applies rcParams, computes figsize at the exact
    journal width, calls plt.subplots with WSPACE / HSPACE tuned per
    layout (1x1, 1x2, 2x1, 2x2, ...) and per profile.  Returns axes as
    a 2D ndarray for uniform indexing regardless of nrows/ncols.

  measure_figure_efficiency(fig)
      -> dict
    Reports axes_area / total_area, per-margin (top/bottom/left/right)
    figure-fraction breakdown, estimated wspace/hspace.  Used by
    emit_paper_figure as the gate metric and by humans to spot waste.

  auto_tighten(fig, target_axes_fraction=0.80, max_iter=8)
      -> dict
    Iteratively shrinks subplots_adjust margins until target reached
    or labels would clip the figure.  Returns the final per-margin
    values + a trace of each iteration's efficiency.  Safe to call
    before savefig.

  add_panel_labels(axes, labels=None, position='top-left',
                   offset_axes=(0.02, 0.97), fontsize=None, **kwargs)
    Place (a), (b), (c)... at consistent in-axes positions across a
    multi-panel figure.  Uses bold sans-serif by IEEE convention.

  emit_paper_figure(fig, path, profile, min_axes_fraction=0.78,
                    on_fail='raise' | 'warn' | 'auto_tighten',
                    dpi=600, ...)
      -> dict
    The validation gate.  Saves PDF + PNG at the profile's exact width
    in mm; measures efficiency; if below floor, behaves per `on_fail`:
      'raise'         (default) -- ValueError + per-margin suggestion
      'warn'          -- warnings.warn() but still saves
      'auto_tighten'  -- run auto_tighten once, re-measure, re-save

  Why this is a separate module from tools.py:
    tools.py covers per-axis helpers (label_curve_endpoints,
    add_slope_guide, ...) and the original digest-oriented presets.
    paper_figure.py is the high-level "I am writing an IEEE / IEEJ /
    Compumag paper, give me a tight figure" entry point that ties the
    style, the layout, and the measurement gate together so the user
    can't accidentally ship a 60%-axes-area figure to reviewers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional

# Lazy matplotlib import — keep the MCP server import-cost low.
# All public functions in this module need matplotlib at runtime; we
# defer the import to the first call so `mcp-server-figure --selftest`
# can run without matplotlib if it ever needs to.
_MPL = None


def _get_mpl():
    """Lazy-import matplotlib + return the (plt, matplotlib) tuple."""
    global _MPL
    if _MPL is None:
        import matplotlib
        import matplotlib.pyplot as plt
        _MPL = (plt, matplotlib)
    return _MPL


def _assert_times_new_roman_available() -> str:
    """Return the resolved Times New Roman font path, or raise.

    Sugahara Lab standard: publication figures are generated in Times
    New Roman.  Silent fallback to DejaVu / Liberation / generic Times
    is not acceptable because it changes the delivered typography and
    can desync figure text from the paper body.
    """
    _, matplotlib = _get_mpl()
    import os
    try:
        path = matplotlib.font_manager.findfont(
            "Times New Roman", fallback_to_default=False
        )
    except Exception as e:
        raise RuntimeError(
            "Times New Roman is required by the Sugahara Lab figure "
            "standard, but matplotlib cannot resolve it. Install Times "
            "New Roman and rebuild the matplotlib font cache "
            "(matplotlib.font_manager._load_fontmanager("
            "try_read_cache=False))."
        ) from e
    if not path or "times" not in os.path.basename(str(path)).lower():
        raise RuntimeError(
            "Times New Roman is required by the Sugahara Lab figure "
            f"standard, but matplotlib resolved {path!r}."
        )
    return str(path)


def _rcparams_request_times_new_roman() -> bool:
    """True when rcParams are configured to render serif text as TNR."""
    plt, _ = _get_mpl()
    fam = plt.rcParams.get("font.family", [])
    fam = [fam] if isinstance(fam, str) else list(fam)
    serifs = plt.rcParams.get("font.serif", [])
    serifs = [serifs] if isinstance(serifs, str) else list(serifs)
    names = [str(x).lower() for x in fam + serifs]
    return "times new roman" in names


# ============================================================
# Color palette: Okabe-Ito 8-color colorblind-safe (Wong 2011 Nature
# Methods).  Lab default for every paper profile.  Distinguishable
# under all common color-vision deficiencies (deuteranopia, protanopia,
# tritanopia) AND when printed in greyscale.
#
# Reference: Wong B (2011) "Color blindness", Nature Methods 8, 441;
#            Okabe & Ito, https://jfly.uni-koeln.de/color/
# Equivalent to SciencePlots `bright` color cycle.
# ============================================================
OKABE_ITO = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky-blue
    "#009E73",  # bluish-green
    "#F0E442",  # yellow         (use sparingly — low contrast on white)
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish-purple
]


def _make_axes_prop_cycle(colors):
    """Build matplotlib's cycler for axes.prop_cycle from a color list."""
    import matplotlib
    from cycler import cycler
    return cycler(color=colors)


# ============================================================
# PaperProfile dataclass + canonical profiles
# ============================================================


@dataclass(frozen=True)
class PaperProfile:
    """Production profile for a target journal column geometry.

    All widths in **millimetres** (exact journal spec, not approximate
    cm).  Font sizes in points.  Aspects below 0.5 = wide, > 0.7 =
    tall; default 0.62 is the golden-ratio-ish IEEE-friendly aspect
    where one axes fills the column nicely.

    The *_tight subplots_adjust values are empirically derived for
    a SINGLE axes filling the column at the profile's chosen font_pt.
    For nrows*ncols > 1, paper_figure() further shrinks the margins
    (the more panels share a row/column, the tighter we can be on the
    interior margins because labels are shared/dropped).
    """

    name: str
    full_name: str
    width_mm: float
    column: Literal["single", "double", "page"]   # for human reference
    font_pt: float
    legend_pt: float
    tick_pt: float
    linewidth_pt: float
    axes_linewidth_pt: float
    marker_size_pt: float
    # Default subplots_adjust for a SINGLE axes.  paper_figure() applies
    # delta-corrections for multi-panel layouts on top of this.
    margin_left: float
    margin_right: float
    margin_top: float
    margin_bottom: float
    # Inter-panel spacing fractions for multi-panel layouts.
    wspace: float = 0.22
    hspace: float = 0.32
    # Default aspect (height / width) for a single axes filling the column.
    default_aspect: float = 0.62
    # Color cycle: Okabe-Ito by default (colorblind-safe + greyscale-safe).
    # Override with a custom list if a project mandates a brand palette.
    color_cycle: tuple = field(default_factory=lambda: tuple(OKABE_ITO))
    # Citation hint.
    spec_url: str = ""

    @property
    def width_in(self) -> float:
        return self.width_mm / 25.4

    @classmethod
    def from_base(cls, *, name, full_name, width_mm, column,
                  base_pt=10.0, small_offset=0.0,
                  linewidth_pt=1.0, axes_linewidth_pt=0.7,
                  marker_size_pt=4.5,
                  margin_left=0.165, margin_right=0.98,
                  margin_top=0.97, margin_bottom=0.20,
                  wspace=0.20, hspace=0.32, default_aspect=0.72,
                  color_cycle=None, spec_url=""):
        """tueplots-style derivation: declare ONE base font; derive
        legend_pt = base_pt and tick_pt = base_pt - small_offset.

        This pins the font-derivation rule so hand-edits to a profile cannot
        accidentally desync the legend or tick sizes. The lab production
        default keeps ``small_offset=0`` because every visible paper label
        must remain at least 10 pt after embedding.

        Args:
            base_pt: body text size (lab rule: 10 pt @ 8 cm).
            small_offset: tick label size = base_pt - small_offset. Keep zero
                for production paper and presentation figures.
            ... (other args mirror the dataclass fields)
        """
        return cls(
            name=name, full_name=full_name, width_mm=width_mm,
            column=column,
            font_pt=base_pt,
            legend_pt=base_pt,                # ALWAYS match body
            tick_pt=base_pt - small_offset,
            linewidth_pt=linewidth_pt,
            axes_linewidth_pt=axes_linewidth_pt,
            marker_size_pt=marker_size_pt,
            margin_left=margin_left, margin_right=margin_right,
            margin_top=margin_top, margin_bottom=margin_bottom,
            wspace=wspace, hspace=hspace,
            default_aspect=default_aspect,
            color_cycle=tuple(color_cycle) if color_cycle else tuple(OKABE_ITO),
            spec_url=spec_url,
        )


# --- IEEE Transactions ---------------------------------------------------

# ABSOLUTE FONT RULE (Sugahara Lab, 2026-05-26):
#
#   The figure-text font is **10 pt** regardless of column width.
#   This is the LAB STANDARD: when the figure is finally embedded at
#   8 cm in the document, the rendered text must be 10 pt so that it
#   matches body-text readability after the figure is printed at
#   100% scale.  Wider columns (e.g. 18 cm double-column) keep the
#   SAME 10 pt -- the axes box grows but the text stays at the
#   absolute readable size, not relative to figure size.
#
# This sets a hard floor on margins: 10 pt = ~3.5 mm tall.  At an
# 8 cm column, an x-axis label + tick labels eat ~13-18 mm of figure
# height (bottom margin ~ 0.16-0.22).  Trying to compress that
# further would shrink the text below 10 pt and break the rule.
#
# DERIVED CONSTANTS at 8 cm column (88.9 mm):
#   font height in fig fraction = (10 pt / 88.9 mm) * (1 pt / 0.3528 mm)
#                                = 3.53 mm / 88.9 mm = 0.0397
#   typical ytick label width (2-3 chars) = ~5 pt wide × 3 = 15 pt
#                                = 5.3 mm / 88.9 mm = 0.060
#   ylabel rotated 90°: needs ~ 0.045 horizontal space (10 pt + padding)
#   tick mark + label + ylabel + padding ~ 0.15 of figure width
#   xlabel + tick labels + tickmark ~ 0.18 of figure height (at default
#                                                            aspect 0.72)


IEEE_SINGLE_COLUMN = PaperProfile(
    name="ieee_single_column",
    full_name="IEEE Transactions, single column",
    width_mm=88.9,                # 3.5 in = 88.9 mm exactly
    column="single",
    font_pt=10.0,                 # LAB STANDARD: 10 pt @ 8 cm column.
                                  #   IEEE Author Center allows 8-10 pt
                                  #   figure text; we pin 10 pt for
                                  #   matched body-text readability.
    legend_pt=10.0,               # Same as body — legends must be
                                  #   readable at print scale, not a
                                  #   shrunken afterthought.
    tick_pt=10.0,
    linewidth_pt=1.0,
    axes_linewidth_pt=0.7,
    marker_size_pt=4.5,
    margin_left=0.165,            # 10 pt ylabel + 9 pt tick labels +
                                  #   tickmark + padding ~ 0.15-0.17
    margin_right=0.98,
    margin_top=0.97,
    margin_bottom=0.200,          # 10 pt xlabel + 9 pt tick labels +
                                  #   tickmark + padding ~ 0.18-0.22 of
                                  #   figure height at default aspect
    wspace=0.20, hspace=0.32,
    default_aspect=0.72,          # taller-than-wide single column is canonical
    spec_url="https://journals.ieeeauthorcenter.ieee.org/"
             "your-publishing-schedule/setting-up-your-document/",
)

IEEE_DOUBLE_COLUMN = PaperProfile(
    name="ieee_double_column",
    full_name="IEEE Transactions, double column (\\figure*)",
    width_mm=181.0,               # 7.16 in = 181.86 mm; 181 mm is the
                                  #   canonical IEEEtran column-pair total
    column="double",
    font_pt=10.0,                 # SAME 10 pt as single column — the
                                  #   font is ABSOLUTE, not relative to
                                  #   figure size.  Wider axes, same
                                  #   readable text.
    legend_pt=10.0,
    tick_pt=10.0,
    linewidth_pt=1.0,
    axes_linewidth_pt=0.7,
    marker_size_pt=4.5,
    # 10 pt absolute on a 181 mm wide figure = ~3.5 mm = 0.019 of width.
    # That means margins as FRACTION OF FIGURE can be much tighter:
    margin_left=0.080,
    margin_right=0.99,
    margin_top=0.97,
    margin_bottom=0.155,          # 10 pt xlabel @ 18 cm = 0.05 of height
                                  #   (aspect 0.38) so bottom ~0.15
    wspace=0.20, hspace=0.30,
    default_aspect=0.38,
    spec_url="IEEE Editorial Style Manual",
)

# --- IEEJ D-section / B-section ------------------------------------------

IEEJ_SINGLE_COLUMN = PaperProfile(
    name="ieej_single_column",
    full_name="IEEJ Trans D / B, single column",
    width_mm=88.0,                # IEEJ規定: 単欄88mm
    column="single",
    font_pt=10.0,                 # LAB STANDARD: 10 pt @ 8 cm column.
    legend_pt=10.0,
    tick_pt=10.0,
    linewidth_pt=1.0,
    axes_linewidth_pt=0.7,
    marker_size_pt=4.5,
    margin_left=0.165,
    margin_right=0.98,
    margin_top=0.97,
    margin_bottom=0.200,
    wspace=0.20, hspace=0.32,
    default_aspect=0.72,
    spec_url="IEEJ 論文投稿の手引き",
)

IEEJ_DOUBLE_COLUMN = PaperProfile(
    name="ieej_double_column",
    full_name="IEEJ Trans D / B, double column (両欄)",
    width_mm=180.0,
    column="double",
    font_pt=10.0,                 # SAME 10 pt; font is absolute.
    legend_pt=10.0,
    tick_pt=10.0,
    linewidth_pt=1.0,
    axes_linewidth_pt=0.7,
    marker_size_pt=4.5,
    margin_left=0.080,
    margin_right=0.99,
    margin_top=0.97,
    margin_bottom=0.155,
    wspace=0.20, hspace=0.30,
    default_aspect=0.38,
    spec_url="IEEJ 論文投稿の手引き",
)

# --- IGTE / Compumag digest ----------------------------------------------

IGTE_DIGEST_DOUBLE = PaperProfile(
    name="igte_digest_double",
    full_name="IGTE / Compumag digest, double-column A4",
    width_mm=170.0,
    column="page",
    font_pt=10.0,                 # SAME 10 pt absolute.
    legend_pt=10.0,
    tick_pt=10.0,
    linewidth_pt=1.0,
    axes_linewidth_pt=0.7,
    marker_size_pt=4.5,
    margin_left=0.085,
    margin_right=0.99,
    margin_top=0.97,
    margin_bottom=0.160,
    wspace=0.20, hspace=0.30,
    default_aspect=0.40,
    spec_url="IGTE Symposium digest template",
)

IGTE_DIGEST_SINGLE = PaperProfile(
    name="igte_digest_single",
    full_name="IGTE / Compumag digest, single-column",
    width_mm=82.0,                # 8.2 cm, close to the 8 cm LAB anchor.
    column="single",
    font_pt=10.0,                 # SAME 10 pt absolute.
    legend_pt=10.0,
    tick_pt=10.0,
    linewidth_pt=1.0,
    axes_linewidth_pt=0.7,
    marker_size_pt=4.5,
    margin_left=0.175,            # 8.2 cm is 8% NARROWER than 8.89, so
                                  #   labels eat slightly MORE figure
                                  #   fraction
    margin_right=0.98,
    margin_top=0.97,
    margin_bottom=0.210,
    wspace=0.20, hspace=0.32,
    default_aspect=0.72,
    spec_url="IGTE Symposium digest template",
)

# --- Beamer 16:9 talk slides (CEFC / IGTE oral) --------------------------
# A beamer ``aspectratio=169`` paper is 160 x 90 mm; the usable text width
# is ~150 mm (full) or ~72 mm (a half-width column inside a two-column
# frame).  Slide figures use a viewing-distance rule distinct from papers:
#   * Times New Roman, authored at 24 pt so every visible label, tick,
#     legend, and annotation remains >=20 pt after actual slide scaling.
#   * NO in-figure title.  The beamer frametitle / slide caption carries
#     it; emit_paper_figure()'s no-title gate fires for these profiles too.
#   * Strokes / markers are heavier than the paper default for projector
#     legibility.
# IMPORTANT: author the figure at the width it will OCCUPY on the slide
# (apply_lab_style(embed_width_cm=...) or these profiles' width_mm) and
# \includegraphics[width=<that width>] at 100%.  Do NOT let
# \includegraphics scale a 150 mm figure into an 8 cm column -- that
# can shrink the 24 pt authored target below the 20 pt displayed floor.

BEAMER_169_FULL = PaperProfile.from_base(
    name="beamer_169_full",
    full_name="Beamer 16:9 slide, full text width",
    width_mm=150.0,
    column="page",
    base_pt=24.0,                 # projected-slide readability floor
    small_offset=0.0,             # ticks must not drop below 24 pt
    linewidth_pt=1.3,             # heavier strokes read better on screen
    axes_linewidth_pt=0.9,
    marker_size_pt=5.0,
    margin_left=0.085, margin_right=0.99,
    margin_top=0.97, margin_bottom=0.16,
    wspace=0.20, hspace=0.30,
    default_aspect=0.50,          # a 16:9 slide favours WIDE figures
    spec_url="beamer aspectratio=169 (160 x 90 mm)",
)

BEAMER_169_HALF = PaperProfile.from_base(
    name="beamer_169_half",
    full_name="Beamer 16:9 slide, half column (two-column frame)",
    width_mm=72.0,
    column="single",
    base_pt=24.0,
    small_offset=0.0,
    linewidth_pt=1.3,
    axes_linewidth_pt=0.9,
    marker_size_pt=5.0,
    margin_left=0.175, margin_right=0.98,
    margin_top=0.97, margin_bottom=0.20,
    wspace=0.20, hspace=0.32,
    default_aspect=0.80,
    spec_url="beamer aspectratio=169 (160 x 90 mm)",
)


PROFILES: dict[str, PaperProfile] = {
    p.name: p for p in [
        IEEE_SINGLE_COLUMN,
        IEEE_DOUBLE_COLUMN,
        IEEJ_SINGLE_COLUMN,
        IEEJ_DOUBLE_COLUMN,
        IGTE_DIGEST_DOUBLE,
        IGTE_DIGEST_SINGLE,
        BEAMER_169_FULL,
        BEAMER_169_HALF,
    ]
}


def get_profile(name_or_profile) -> PaperProfile:
    """Accept either a string key or a PaperProfile object."""
    if isinstance(name_or_profile, PaperProfile):
        return name_or_profile
    key = str(name_or_profile).strip().lower()
    if key not in PROFILES:
        raise ValueError(
            f"Unknown profile {name_or_profile!r}.  Available: "
            f"{', '.join(sorted(PROFILES))}"
        )
    return PROFILES[key]


# ============================================================
# Layout scaffold
# ============================================================

# Per-panel-count delta-corrections to the profile's default margins.
# Multi-panel layouts can take MORE of the figure for axes because
# (a) shared internal labels are dropped / smaller, (b) we're willing
# to push exterior margins tighter because more axes amortize.
_MARGIN_DELTAS = {
    # (nrows, ncols) -> dict of subplots_adjust corrections (added to
    # profile defaults).  Positive = looser, negative = tighter.
    (1, 1): dict(),
    (1, 2): dict(margin_left=-0.020, margin_right=+0.005,
                 wspace=-0.02),
    (1, 3): dict(margin_left=-0.030, margin_right=+0.005,
                 wspace=-0.04),
    (2, 1): dict(margin_top=+0.005, margin_bottom=-0.020,
                 hspace=-0.05),
    (2, 2): dict(margin_left=-0.020, margin_right=+0.005,
                 margin_top=+0.005, margin_bottom=-0.020,
                 wspace=-0.02, hspace=-0.05),
}


def paper_figure(
    profile,
    nrows: int = 1,
    ncols: int = 1,
    aspect: Optional[float] = None,
    rel_width: float = 1.0,
    figsize_override_in: Optional[tuple[float, float]] = None,
    sharex: bool = False,
    sharey: bool = False,
    use_times_roman: bool = True,
    panel_labels="auto",
    panel_label_kwargs: Optional[dict] = None,
):
    """Create a paper-quality (fig, axes_2d) at the journal's exact width.

    Wraps:
        apply_lab_style() -> exact-width figure -> plt.subplots ->
        per-layout subplots_adjust -> optional add_panel_labels.

    Args:
        profile: PaperProfile or string key (e.g. 'ieee_double_column').
        nrows, ncols: grid shape.
        aspect: height / width ratio of the WHOLE figure (default uses
            profile.default_aspect; for nrows > 1 it is multiplied by
            nrows to keep each row at the per-axes aspect).
        rel_width: fraction of the profile's journal width.  Default
            1.0 = full column.  Use 1.5 for a 1.5-column figure
            (ieee_single_column * 1.5 ~ 13 cm), or 0.5 for a
            half-column inset.  Allows custom column widths without
            forking a new profile.  Tueplots pattern.
        figsize_override_in: explicit (width_in, height_in); overrides
            both rel_width AND the journal-width derivation.  Use for
            truly non-standard layouts only.
        sharex, sharey: forwarded to plt.subplots.
        use_times_roman: enforce Times New Roman, the Sugahara Lab
            standard.  False is only for deliberate non-public drafts.
        panel_labels: "auto" (default) places (a), (b), (c) labels
            automatically iff nrows*ncols > 1.  True forces labels
            even on single-panel; False suppresses them.  Mirrors the
            ultraplot `abc=True` ergonomics.
        panel_label_kwargs: extra kwargs for add_panel_labels.

    Returns:
        (fig, axes_2d).  axes_2d is ALWAYS a (nrows, ncols) ndarray of
        Axes objects, even for nrows=ncols=1 (you get axes_2d[0,0]).
        This is intentional — it lets the same downstream loop work
        regardless of layout.
    """
    import numpy as np
    plt, matplotlib = _get_mpl()
    prof = get_profile(profile)

    # --- Figure size ---
    if figsize_override_in is not None:
        figsize_in = figsize_override_in
    else:
        w_in = prof.width_in * float(rel_width)
        a = aspect if aspect is not None else prof.default_aspect
        # For nrows > 1, scale the height up so each ROW is the
        # per-axes aspect.  This avoids the "2x1 looks squashed because
        # default aspect is for 1x1" trap.
        h_in = w_in * a * (nrows if nrows > 1 else 1)
        figsize_in = (w_in, h_in)

    # --- rcParams ---
    rc = {
        "font.size": prof.font_pt,
        "axes.labelsize": prof.font_pt,
        "axes.titlesize": prof.font_pt,
        "xtick.labelsize": prof.tick_pt,
        "ytick.labelsize": prof.tick_pt,
        "legend.fontsize": prof.legend_pt,
        "axes.linewidth": prof.axes_linewidth_pt,
        "lines.linewidth": prof.linewidth_pt,
        "lines.markersize": prof.marker_size_pt,
        "axes.grid": True,
        "grid.linestyle": ":",
        "grid.alpha": 0.4,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "legend.frameon": False,
        "pdf.fonttype": 42,       # TrueType — IEEE/Elsevier require this
        "ps.fonttype": 42,
        # Colorblind-safe palette (Okabe-Ito): lab default for every
        # profile.  Override via prof.color_cycle.
        "axes.prop_cycle": _make_axes_prop_cycle(list(prof.color_cycle)),
    }
    if use_times_roman:
        _assert_times_new_roman_available()
        rc["font.family"] = "serif"
        rc["font.serif"] = ["Times New Roman"]
        rc["mathtext.fontset"] = "stix"
        rc["mathtext.default"] = "it"
    plt.rcParams.update(rc)

    # --- Figure + subplots ---
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=figsize_in,
        sharex=sharex, sharey=sharey,
        squeeze=False,            # always return 2D ndarray
        dpi=600,
    )
    fig._lab_medium = (
        "presentation" if prof.name.startswith("beamer_169_") else "paper"
    )
    fig._lab_authored_width_cm = float(figsize_in[0]) * 2.54
    fig._lab_embed_width_cm = float(figsize_in[0]) * 2.54
    # axes is already 2D thanks to squeeze=False.

    # --- subplots_adjust per profile + layout delta ---
    deltas = _MARGIN_DELTAS.get((nrows, ncols), {})
    fig.subplots_adjust(
        left=prof.margin_left + deltas.get("margin_left", 0),
        right=prof.margin_right + deltas.get("margin_right", 0),
        top=prof.margin_top + deltas.get("margin_top", 0),
        bottom=prof.margin_bottom + deltas.get("margin_bottom", 0),
        wspace=prof.wspace + deltas.get("wspace", 0),
        hspace=prof.hspace + deltas.get("hspace", 0),
    )

    # --- Optional panel labels ---
    # `"auto"` -> apply when multi-panel.  True -> always.  False -> never.
    auto_apply = (panel_labels == "auto" and (nrows * ncols) > 1)
    forced = (panel_labels is True)
    if auto_apply or forced:
        add_panel_labels(axes, **(panel_label_kwargs or {}))

    return fig, axes


# ============================================================
# add_panel_labels — (a), (b), (c) for multi-panel figures
# ============================================================


def add_panel_labels(
    axes,
    labels: Optional[list[str]] = None,
    position: Literal["top-left", "top-right", "bottom-left", "bottom-right"]
        = "top-left",
    offset_axes: tuple[float, float] = (0.04, 0.96),
    fontsize: Optional[float] = None,
    fontweight: str = "bold",
    paren: str = "()",
    color: str = "black",
    bbox_facecolor: Optional[str] = None,
    bbox_pad: float = 0.15,
):
    """Place (a), (b), (c)... at consistent positions across axes.

    Args:
        axes: 1D or 2D iterable of Axes (e.g. the array returned by
            paper_figure).  Order: row-major.
        labels: Override the auto (a), (b)... sequence with arbitrary
            strings (e.g. ["(i)", "(ii)", "(iii)"]).
        position: Anchor corner inside each axes.
        offset_axes: (x_axes_frac, y_axes_frac) inside the axes.
            Default (0.04, 0.96) = 4% from the left edge, 96% from
            bottom = effectively top-left corner with a small inset.
        fontsize: Override; defaults to rcParams['axes.labelsize'].
        fontweight: 'bold' is the IEEE convention.
        paren: Pair of characters wrapping the letter.  Default '()'.
            Use '[]' for some IEEJ styles, '' for no parens.
        color: Text color.
        bbox_facecolor: If not None, draw a colored bbox behind the
            label (useful when overlapping data).  Default None = no
            bbox (cleanest look).
        bbox_pad: Padding inside the bbox.

    Returns:
        list of Text objects, in axes-order.
    """
    plt, _ = _get_mpl()
    import numpy as np
    arr = np.asarray(axes).ravel()
    if labels is None:
        labels = [f"{paren[:1] if paren else ''}{chr(ord('a') + i)}"
                  f"{paren[1:] if paren else ''}"
                  for i in range(len(arr))]
    if len(labels) < len(arr):
        raise ValueError(
            f"add_panel_labels: {len(labels)} labels provided but "
            f"{len(arr)} axes — pass `labels` of matching length or "
            f"None for auto-generation."
        )

    ha_va = {
        "top-left":     ("left", "top",
                         offset_axes[0],         offset_axes[1]),
        "top-right":    ("right", "top",
                         1.0 - offset_axes[0],   offset_axes[1]),
        "bottom-left":  ("left", "bottom",
                         offset_axes[0],         1.0 - offset_axes[1]),
        "bottom-right": ("right", "bottom",
                         1.0 - offset_axes[0],   1.0 - offset_axes[1]),
    }
    if position not in ha_va:
        raise ValueError(f"position must be one of {sorted(ha_va)}, got "
                         f"{position!r}")
    ha, va, xa, ya = ha_va[position]
    bbox_kw = None
    if bbox_facecolor is not None:
        bbox_kw = dict(facecolor=bbox_facecolor, edgecolor="none",
                       pad=bbox_pad)

    if fontsize is None:
        fontsize = plt.rcParams["axes.labelsize"]

    texts = []
    for ax, lbl in zip(arr, labels):
        t = ax.text(xa, ya, lbl, transform=ax.transAxes,
                    ha=ha, va=va,
                    fontsize=fontsize, fontweight=fontweight,
                    color=color, bbox=bbox_kw)
        texts.append(t)
    return texts


# ============================================================
# measure_figure_efficiency — gate input metric
# ============================================================


def measure_figure_efficiency(fig) -> dict:
    """Compute axes-area / total-area + per-margin breakdown.

    The key metric is ``axes_area_fraction``: the sum of the axes
    bounding boxes' area divided by the figure area.  A well-designed
    paper figure is typically 0.78-0.88 efficient (less than 0.70
    typically signals wasted margin).

    Includes detail useful for the auto-tighten loop and for human
    diagnostic:
      - per-margin fraction of figure NOT used by axes (top/bottom/
        left/right of the BOUNDING box of all axes)
      - estimated wspace / hspace (median fractional gap between
        adjacent axes columns/rows)
      - per-axes bbox + area
    """
    plt, _ = _get_mpl()
    fig.canvas.draw()

    fw, fh = fig.get_size_inches()
    fig_area = fw * fh
    axes = fig.get_axes()
    if not axes:
        return dict(
            fig_size_inches=(fw, fh),
            fig_area_in2=fig_area,
            n_axes=0,
            axes_area_total_in2=0.0,
            axes_area_fraction=0.0,
            margin_breakdown={"left": 1.0, "right": 0.0,
                                "top": 1.0, "bottom": 0.0},
            wspace_estimated=None,
            hspace_estimated=None,
            axes=[],
        )

    # Collect axes positions in figure-fraction coordinates.
    positions = []   # list of (x0, y0, x1, y1) in figure fraction
    axes_info = []
    for i, ax in enumerate(axes):
        bbox = ax.get_position()
        x0, y0, w, h = bbox.x0, bbox.y0, bbox.width, bbox.height
        positions.append((x0, y0, x0 + w, y0 + h))
        axes_info.append({
            "index": i,
            "bbox_fig_frac": (x0, y0, x0 + w, y0 + h),
            "area_in2": w * h * fig_area,
            "width_in": w * fw,
            "height_in": h * fh,
        })

    axes_area = sum(a["area_in2"] for a in axes_info)
    axes_fraction = axes_area / fig_area if fig_area > 0 else 0.0

    # Per-margin (figure fraction) — relative to the BOUNDING box of
    # all axes (so for nrows*ncols > 1, "left" is the leftmost x0).
    xs0 = [p[0] for p in positions]
    ys0 = [p[1] for p in positions]
    xs1 = [p[2] for p in positions]
    ys1 = [p[3] for p in positions]
    margin = {
        "left":   min(xs0),
        "right":  1.0 - max(xs1),
        "top":    1.0 - max(ys1),
        "bottom": min(ys0),
    }

    # Estimate wspace, hspace by looking at neighbour columns/rows.
    # Sort axes-columns by x0; the gap between column N's x1 and
    # column N+1's x0 is the horizontal subplot gap.  Convert to
    # fraction-of-mean-axes-width (matplotlib's wspace convention).
    wspace_est = None
    if len(positions) > 1:
        col_centers = sorted({round(p[0], 3) for p in positions})
        if len(col_centers) > 1:
            # Pair adjacent column starts; subtract the column's width
            x0s_sorted = sorted({p[0] for p in positions})
            gaps = []
            widths = []
            for i in range(len(x0s_sorted) - 1):
                # Find the rightmost x1 among axes whose x0 == x0s_sorted[i]
                xs_right = [p[2] for p in positions
                            if abs(p[0] - x0s_sorted[i]) < 1e-6]
                if not xs_right:
                    continue
                xr = max(xs_right)
                gap = x0s_sorted[i + 1] - xr
                if gap > 0:
                    gaps.append(gap)
                    widths.append(xr - x0s_sorted[i])
            if gaps and widths:
                mean_w = sum(widths) / len(widths)
                if mean_w > 0:
                    wspace_est = (sum(gaps) / len(gaps)) / mean_w

    hspace_est = None
    if len(positions) > 1:
        row_centers = sorted({round(p[1], 3) for p in positions})
        if len(row_centers) > 1:
            y0s_sorted = sorted({p[1] for p in positions})
            gaps = []
            heights = []
            for i in range(len(y0s_sorted) - 1):
                ys_top = [p[3] for p in positions
                          if abs(p[1] - y0s_sorted[i]) < 1e-6]
                if not ys_top:
                    continue
                yt = max(ys_top)
                gap = y0s_sorted[i + 1] - yt
                if gap > 0:
                    gaps.append(gap)
                    heights.append(yt - y0s_sorted[i])
            if gaps and heights:
                mean_h = sum(heights) / len(heights)
                if mean_h > 0:
                    hspace_est = (sum(gaps) / len(gaps)) / mean_h

    return dict(
        fig_size_inches=(fw, fh),
        fig_area_in2=fig_area,
        n_axes=len(axes_info),
        axes_area_total_in2=axes_area,
        axes_area_fraction=axes_fraction,
        margin_breakdown=margin,
        wspace_estimated=wspace_est,
        hspace_estimated=hspace_est,
        axes=axes_info,
    )


# ============================================================
# auto_tighten — iterative margin shrinker
# ============================================================


def _compute_overhang_px(fig) -> dict:
    """Per-side maximum pixel overhang of any text artist past fig.bbox.

    Used by auto_tighten() as a **delta** check.  Matplotlib reports
    tick-label bboxes that routinely extend a few dozen pixels past
    the figure boundary even on an untouched figure (because the bbox
    accounts for descender / ascender / centering quirks).  An absolute
    "no text outside fig.bbox" rule would therefore reject every
    tightening on a fresh paper_figure().  Comparing to baseline
    overhang plus a small tolerance fixes that.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fb = fig.bbox
    over = {"left": 0.0, "right": 0.0, "top": 0.0, "bottom": 0.0}
    for ax in fig.get_axes():
        artists = [ax.xaxis.label, ax.yaxis.label, ax.title]
        artists += list(ax.get_xticklabels())
        artists += list(ax.get_yticklabels())
        for t in artists:
            if not t.get_text():
                continue
            try:
                bb = t.get_window_extent(renderer=renderer)
            except RuntimeError:
                continue
            over["left"]   = max(over["left"],   fb.x0 - bb.x0)
            over["right"]  = max(over["right"],  bb.x1 - fb.x1)
            over["bottom"] = max(over["bottom"], fb.y0 - bb.y0)
            over["top"]    = max(over["top"],    bb.y1 - fb.y1)
    return over


# Per-side overhang delta tolerances for the auto_tighten safety check.
# Use FIGURE-FRACTION tolerances rather than absolute pixels so the
# behaviour is dpi-independent (a 100-dpi vs 600-dpi figure produces
# 6x more pixels for the same physical layout, but the "label is too
# close to the edge" decision should be the same).
#
# The tolerance absorbs sub-millimetre bbox drift on top of the
# BASELINE snapshot (which already accounts for matplotlib reporting
# tick-label bboxes a few dozen pixels outside an untouched figure).
# It MUST stay well below the size of a text label: a tolerance wider
# than the label itself lets auto_tighten push a whole axis label off
# the canvas and still call it "within tolerance".
#
# MEASURED 2026-08-10, and the reason these were cut from 0.02/0.03:
# on a 181 mm 2x1 ieee_double_column with y labels, the old 2%-of-width
# tolerance was 85 px @600 dpi = 3.6 mm -- LARGER than a 10 pt label --
# so auto_tighten(0.78) shrank left 0.080 -> 0.040 and left BOTH y
# labels 46-65 px OUTSIDE the canvas while every gate reported clean.
# 0.4% of width = 0.72 mm on that figure: enough for bbox rounding,
# far too little to lose a label.
_OVERHANG_TOL_X_FRAC = 0.004
_OVERHANG_TOL_Y_FRAC = 0.006


def auto_tighten(
    fig,
    target_axes_fraction: float = 0.80,
    max_iter: int = 24,
    step: float = 0.005,
    sides: tuple[str, ...] = ("left", "right", "top", "bottom",
                                "wspace", "hspace"),
    min_margin: float = 0.005,
    max_wspace_hspace_reduction: float = 0.10,
    overhang_tol_x_frac: float = _OVERHANG_TOL_X_FRAC,
    overhang_tol_y_frac: float = _OVERHANG_TOL_Y_FRAC,
) -> dict:
    """Iteratively shrink subplots_adjust margins to hit target efficiency.

    Algorithm:
      0. Snapshot baseline per-side **overhang** of text artists past
         fig.bbox (matplotlib reports these even on untouched figures).
      1. Measure current efficiency.  If >= target, stop.
      2. For each side in `sides`, try shrinking by `step`; accept iff
         the resulting per-side overhang did NOT increase past
         (baseline + overhang_tolerance_px).  Multiple sides per
         iteration are accepted (not first-success-only) so a single
         iter can shrink left + bottom + wspace together.

    Stops when no side can tighten safely, target reached, or max_iter.

    Args:
        fig: matplotlib Figure (mutated in place).
        target_axes_fraction: stop when reached.
        max_iter: hard cap on iterations.
        step: per-iteration shrinkage per side (figure fraction).
        sides: which margins to consider.
        min_margin: floor on outer margins (avoid pushing labels off
            the page entirely).
        max_wspace_hspace_reduction: cumulative cap on wspace/hspace
            shrinkage (default 0.10, e.g. wspace 0.20 -> 0.10 max).
        overhang_tolerance_px: per-side overhang delta tolerance (px).

    Returns:
        dict with:
          - final_efficiency (float)
          - reached_target (bool)
          - trace: list of {iter, efficiency, applied_sides:[str],
                            rejected:[(side, reason)]}
          - final_adjust (dict)
          - initial_adjust (dict)
          - baseline_overhang_px (dict per side)
    """
    plt, _ = _get_mpl()

    def _current_adjust():
        sp = fig.subplotpars
        return dict(
            left=sp.left, right=sp.right,
            top=sp.top, bottom=sp.bottom,
            wspace=sp.wspace, hspace=sp.hspace,
        )

    initial = _current_adjust()
    baseline_overhang = _compute_overhang_px(fig)
    eff0 = measure_figure_efficiency(fig)["axes_area_fraction"]
    trace = [{"iter": 0, "efficiency": eff0,
              "applied_sides": [], "rejected": []}]
    if eff0 >= target_axes_fraction:
        return dict(
            final_efficiency=eff0,
            reached_target=True,
            trace=trace,
            final_adjust=_current_adjust(),
            initial_adjust=initial,
            baseline_overhang_px=baseline_overhang,
        )

    wspace_reduced = 0.0
    hspace_reduced = 0.0

    fb_w_px = fig.bbox.width
    fb_h_px = fig.bbox.height
    tol_x_px = overhang_tol_x_frac * fb_w_px
    tol_y_px = overhang_tol_y_frac * fb_h_px

    def _safe_delta() -> tuple[bool, str]:
        """(True, '') iff overhang within baseline + per-axis tolerance."""
        cur = _compute_overhang_px(fig)
        for s in ("left", "right"):
            if cur[s] > baseline_overhang[s] + tol_x_px:
                return False, (
                    f"{s} overhang {cur[s]:.0f}px > baseline "
                    f"{baseline_overhang[s]:.0f}px + {tol_x_px:.0f}px tol "
                    f"({overhang_tol_x_frac*100:.1f}% of fig width)"
                )
        for s in ("top", "bottom"):
            if cur[s] > baseline_overhang[s] + tol_y_px:
                return False, (
                    f"{s} overhang {cur[s]:.0f}px > baseline "
                    f"{baseline_overhang[s]:.0f}px + {tol_y_px:.0f}px tol "
                    f"({overhang_tol_y_frac*100:.1f}% of fig height)"
                )
        return True, ""

    for it in range(1, max_iter + 1):
        applied_sides = []
        rejected = []
        for side in sides:
            cur = _current_adjust()
            saved = dict(cur)
            adj = dict(cur)
            # IMPORTANT: subplots_adjust(left=X) means "axes START at X
            # fraction from the LEFT" -- so INCREASING left WIDENS the
            # left margin (less axes area).  To tighten, DECREASE left.
            # Symmetric reasoning for right (=axes-end x-fraction),
            # top (=axes-end y-fraction from BOTTOM), and bottom
            # (=axes-start y-fraction from bottom).  We tighten by
            # DEcreasing left / bottom and INcreasing right / top.
            if side == "left":
                adj["left"] = max(adj["left"] - step, min_margin)
                if adj["left"] >= cur["left"]:
                    continue   # already at floor, no shrink possible
            elif side == "right":
                adj["right"] = min(adj["right"] + step, 1.0 - min_margin)
                if adj["right"] <= cur["right"]:
                    continue
            elif side == "top":
                adj["top"] = min(adj["top"] + step, 1.0 - min_margin)
                if adj["top"] <= cur["top"]:
                    continue
            elif side == "bottom":
                adj["bottom"] = max(adj["bottom"] - step, min_margin)
                if adj["bottom"] >= cur["bottom"]:
                    continue
            elif side == "wspace":
                if wspace_reduced + step > max_wspace_hspace_reduction:
                    rejected.append((side, "wspace_cap_reached"))
                    continue
                adj["wspace"] = max(adj["wspace"] - step, 0.0)
            elif side == "hspace":
                if hspace_reduced + step > max_wspace_hspace_reduction:
                    rejected.append((side, "hspace_cap_reached"))
                    continue
                adj["hspace"] = max(adj["hspace"] - step, 0.0)
            else:
                continue

            fig.subplots_adjust(**adj)
            ok, why = _safe_delta()
            if ok:
                if side == "wspace":
                    wspace_reduced += step
                if side == "hspace":
                    hspace_reduced += step
                applied_sides.append(side)
            else:
                fig.subplots_adjust(**saved)
                rejected.append((side, why))

        eff = measure_figure_efficiency(fig)["axes_area_fraction"]
        trace.append({"iter": it, "efficiency": eff,
                       "applied_sides": applied_sides,
                       "rejected": rejected})
        if not applied_sides:
            break
        if eff >= target_axes_fraction:
            break

    final_eff = measure_figure_efficiency(fig)["axes_area_fraction"]
    return dict(
        final_efficiency=final_eff,
        reached_target=final_eff >= target_axes_fraction,
        trace=trace,
        final_adjust=_current_adjust(),
        initial_adjust=initial,
        baseline_overhang_px=baseline_overhang,
    )


# ============================================================
# emit_paper_figure — the validation gate
# ============================================================


def _check_colors_are_cvd_safe(fig, allowed_palette=None) -> list[dict]:
    """Lint Line2D / Patch / Collection colors against an
    allowed colorblind-vision-deficient (CVD) safe palette.

    Any artist whose color falls OUTSIDE the allowed palette (within
    a small euclidean RGB tolerance) is flagged.  Returns a list of
    {axis, label, color_hex, reason}.  Empty list = clean.

    The default palette is Okabe-Ito 8-color (Wong 2011, Nature
    Methods 8:441).  Pass `allowed_palette=[...]` to override (e.g.
    a project's brand palette).

    Notes:
      - Matplotlib's default `tab10` cycle is NOT CVD-safe (orange and
        red are confusable in deuteranopia).  Lab default is Okabe-Ito.
      - Greyscale colors (any rgb where r==g==b) are always allowed
        because they survive every CVD condition AND print well.
    """
    plt, _ = _get_mpl()
    import matplotlib.colors as mcolors
    if allowed_palette is None:
        allowed_palette = OKABE_ITO
    allowed_rgb = [mcolors.to_rgb(c) for c in allowed_palette]
    tol = 0.04   # ~10/255 per channel — generous to allow alpha-tinted
                 # variants of the palette
    out = []

    def _is_grey(rgb):
        return abs(rgb[0] - rgb[1]) < 0.02 and abs(rgb[1] - rgb[2]) < 0.02

    def _is_in_palette(rgb):
        if _is_grey(rgb):
            return True
        for ar in allowed_rgb:
            if (abs(rgb[0] - ar[0]) < tol
                and abs(rgb[1] - ar[1]) < tol
                and abs(rgb[2] - ar[2]) < tol):
                return True
        return False

    for i, ax in enumerate(fig.get_axes()):
        for line in ax.get_lines():
            try:
                rgb = mcolors.to_rgb(line.get_color())
            except (ValueError, TypeError):
                continue
            label = line.get_label() or ""
            if label.startswith("_"):
                continue
            if not _is_in_palette(rgb):
                out.append({
                    "axis": i, "label": label,
                    "color": mcolors.to_hex(rgb),
                    "reason": (
                        f"color {mcolors.to_hex(rgb)} not in "
                        f"Okabe-Ito CVD-safe palette and is not greyscale"
                    ),
                })
    return out


def _check_pdf_fonts_embedded(pdf_path) -> list[str]:
    """Open the saved PDF and verify Type-1/Type-3 fonts are NOT present.

    Type-1 and Type-3 fonts are rejection-worthy for IEEE / Elsevier
    (Type-3 is raster glyph; Type-1 has poor unicode support).  Only
    TrueType (Type-42) is accepted.  Returns a list of violation
    descriptions; empty list = clean.

    Implementation: scan the raw PDF stream for `/Subtype /Type1` or
    `/Subtype /Type3`.  This avoids depending on pdffonts/poppler
    being installed.
    """
    from pathlib import Path
    p = Path(pdf_path)
    if not p.exists() or p.suffix.lower() != ".pdf":
        return []  # nothing to check
    violations = []
    try:
        with open(p, "rb") as f:
            blob = f.read()
    except OSError as e:
        violations.append(f"could not read {p}: {e}")
        return violations
    if b"/Subtype /Type3" in blob or b"/Subtype/Type3" in blob:
        violations.append(
            f"{p.name} contains Type-3 (raster) fonts — IEEE / Elsevier"
            f" pre-flight will reject. Confirm rcParams['pdf.fonttype']"
            f" == 42 BEFORE savefig()."
        )
    # Type-1 (PostScript) — accept only if matched with TrueType
    # variant; pure Type1-only embed is rejection-worthy at modern
    # publishers.  We just flag it.
    has_type1 = b"/Subtype /Type1" in blob or b"/Subtype/Type1" in blob
    has_truetype = b"/Subtype /TrueType" in blob or b"/Subtype/TrueType" in blob
    if has_type1 and not has_truetype:
        violations.append(
            f"{p.name} contains Type-1 fonts but no TrueType — some"
            f" publishers (Elsevier, IEEE) prefer Type-42 TrueType."
        )
    return violations


_BAD_SERIF_FALLBACK_FONT_FRAGMENTS = (
    b"DejaVuSerif", b"DejaVu Serif",
    b"LiberationSerif", b"Liberation Serif",
    b"NimbusRomNo9", b"Nimbus Roman",
)


def _check_pdf_times_new_roman(pdf_path) -> list[str]:
    """Post-save check for common serif fallbacks in the final PDF."""
    from pathlib import Path
    p = Path(pdf_path)
    if not p.exists() or p.suffix.lower() != ".pdf":
        return []
    try:
        blob = p.read_bytes()
    except OSError as e:
        return [f"could not read {p}: {e}"]
    hits = sorted({
        frag.decode("ascii", "replace")
        for frag in _BAD_SERIF_FALLBACK_FONT_FRAGMENTS
        if frag in blob
    })
    if hits:
        return [
            f"{p.name} embeds serif fallback font(s): {', '.join(hits)}. "
            "Sugahara Lab figures must be generated with Times New Roman, "
            "not a fallback family."
        ]
    return []


def _contains_cjk(s: str) -> bool:
    """True if the string contains any Japanese / CJK character.

    Covers Hiragana, Katakana (incl. halfwidth), CJK Unified Ideographs
    (+ Extension A), CJK symbols/punctuation, and fullwidth forms.  Used
    to detect Japanese text in figure labels — a paper-figure must use
    ENGLISH labels only (Japanese text forces the PDF to embed a CJK
    font, which English-language publishers reject / cannot render).
    """
    for ch in s:
        cp = ord(ch)
        if (0x3040 <= cp <= 0x30FF        # Hiragana + Katakana
                or 0x3400 <= cp <= 0x4DBF  # CJK Ext A
                or 0x4E00 <= cp <= 0x9FFF  # CJK Unified Ideographs
                or 0x3000 <= cp <= 0x303F  # CJK symbols & punctuation
                or 0xFF00 <= cp <= 0xFFEF):  # fullwidth/halfwidth forms
            return True
    return False


def _check_no_japanese_text(fig) -> list[dict]:
    """Lab rule: paper figures carry NO Japanese text.

    Scans every text artist (axis labels, tick labels, titles, legend
    entries, annotations / ax.text) for CJK characters.  Any hit is a
    violation: the figure must be re-labelled in English BEFORE it goes
    into an IEEE / IEEJ-English / IGTE / Compumag paper.

    Why this is enforced on the GRAPH side (not just paper-writing):
      - The root cause is the matplotlib figure having Japanese labels.
      - If we let it through to savefig, the PDF embeds a Japanese font
        (MS Gothic / Yu Gothic / IPAex / Noto CJK), bloating the file
        and failing the publisher's font pre-flight.
      - Catching it at emit time (before the PDF exists) gives a 1-line
        fix: re-label the axis in English.

    Returns a list of {kind, index, text} dicts; empty = clean.
    """
    out = []
    for i, ax in enumerate(fig.get_axes()):
        # axis labels
        for kind, getter in (("xlabel", ax.get_xlabel),
                              ("ylabel", ax.get_ylabel),
                              ("title", ax.get_title)):
            txt = getter()
            if txt and _contains_cjk(txt):
                out.append({"kind": kind, "axis": i, "text": txt})
        # tick labels (only those explicitly set with text; auto numeric
        # ticks have no CJK)
        for kind, ticks in (("xticklabel", ax.get_xticklabels()),
                            ("yticklabel", ax.get_yticklabels())):
            for t in ticks:
                s = t.get_text()
                if s and _contains_cjk(s):
                    out.append({"kind": kind, "axis": i, "text": s})
        # legend entries
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                s = t.get_text()
                if s and _contains_cjk(s):
                    out.append({"kind": "legend", "axis": i, "text": s})
        # free-standing ax.text / annotate
        for t in ax.texts:
            s = t.get_text()
            if s and _contains_cjk(s):
                out.append({"kind": "text", "axis": i, "text": s})
    # figure-level suptitle
    sup = getattr(fig, "_suptitle", None)
    if sup is not None and _contains_cjk(sup.get_text()):
        out.append({"kind": "suptitle", "axis": -1, "text": sup.get_text()})
    return out


def _check_text_overflow(fig, tol_pt: float = 0.5) -> list[dict]:
    """Return free-text artists whose drawn bbox exits the figure canvas.

    This complements the axis-label clipping checks in ``_lab_api``.
    Concept diagrams often use ``ax.text(...)`` / ``annotate(...)``
    instead of axis labels; those labels can be cut off while ordinary
    plot gates still pass.  We scan only free diagram text here:
    axis labels, tick labels, titles, and legend entries are handled by
    their own dedicated gates.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    dpi = float(fig.dpi)
    tol_px = tol_pt * dpi / 72.0
    fb = fig.bbox
    out = []

    def _record(artist, *, axis, kind):
        if artist is None or not artist.get_visible():
            return
        text = artist.get_text()
        if not text:
            return
        try:
            bb = artist.get_window_extent(renderer=renderer)
        except Exception:
            return
        for side, over_px in (
            ("left", fb.x0 - bb.x0),
            ("bottom", fb.y0 - bb.y0),
            ("right", bb.x1 - fb.x1),
            ("top", bb.y1 - fb.y1),
        ):
            if over_px > tol_px:
                out.append({
                    "axis": axis,
                    "kind": kind,
                    "text": text,
                    "side": side,
                    "overhang_pt": round(over_px * 72.0 / dpi, 2),
                })

    for i, ax in enumerate(fig.get_axes()):
        excluded = {
            id(ax.xaxis.label),
            id(ax.yaxis.label),
            id(ax.title),
            id(getattr(ax, "_left_title", None)),
            id(getattr(ax, "_right_title", None)),
        }
        excluded.update(id(t) for t in ax.get_xticklabels())
        excluded.update(id(t) for t in ax.get_yticklabels())
        leg = ax.get_legend()
        if leg is not None:
            excluded.update(id(t) for t in leg.get_texts())
        for t in ax.texts:
            if id(t) not in excluded:
                _record(t, axis=i, kind="text")

    for t in fig.texts:
        if t is getattr(fig, "_suptitle", None):
            continue
        _record(t, axis=-1, kind="figure_text")
    return out


def audit_text_overflow(fig, tol_pt: float = 0.5) -> list[dict]:
    """Lint free diagram text for canvas overflow.

    Returns ``[{"axis", "kind", "text", "side", "overhang_pt"}, ...]``;
    an empty list means every ``ax.text`` / ``fig.text`` label fits
    inside the fixed figure canvas.
    """
    return _check_text_overflow(fig, tol_pt=tol_pt)


# Japanese / CJK font-name fragments that, if present in a saved PDF's
# font dictionary, mean a CJK glyph set got embedded.  Curated to avoid
# false positives on Latin families (e.g. "Century Gothic" is NOT here;
# we match the Japanese-specific families only).
_JAPANESE_FONT_FRAGMENTS = (
    b"MS-Gothic", b"MSGothic", b"MS Gothic",
    b"MS-Mincho", b"MSMincho", b"MS Mincho",
    b"MS-PGothic", b"MS-PMincho",
    b"Meiryo", b"YuGothic", b"Yu Gothic", b"YuMincho", b"Yu Mincho",
    b"IPAGothic", b"IPAMincho", b"IPAexGothic", b"IPAexMincho", b"IPAex",
    b"NotoSansCJK", b"NotoSerifCJK", b"Noto Sans CJK", b"Noto Serif CJK",
    b"NotoSansJP", b"NotoSerifJP",
    b"Hiragino", b"HiraKaku", b"HiraMin", b"HiraginoSans",
    b"TakaoGothic", b"TakaoMincho", b"TakaoPGothic",
    b"Osaka", b"Kozuka", b"KozGoPr", b"KozMinPr",
    b"SourceHanSans", b"SourceHanSerif", b"Source Han",
    b"GenShinGothic", b"Migu", b"Sazanami",
)


def _check_no_japanese_font_in_pdf(pdf_path) -> list[str]:
    """Belt-and-suspenders: scan a saved PDF for embedded Japanese fonts.

    Complements ``_check_no_japanese_text`` (which scans the live matplotlib
    text artists).  Even when artists look clean, a stray glyph or a
    mathtext fallback could pull a CJK font; this confirms the FINAL
    bytes carry no Japanese family.  Returns violation strings; empty =
    clean.
    """
    from pathlib import Path
    p = Path(pdf_path)
    if not p.exists() or p.suffix.lower() != ".pdf":
        return []
    try:
        with open(p, "rb") as f:
            blob = f.read()
    except OSError as e:
        return [f"could not read {p}: {e}"]
    hits = sorted({frag.decode("ascii", "replace")
                   for frag in _JAPANESE_FONT_FRAGMENTS if frag in blob})
    if hits:
        return [
            f"{p.name} embeds Japanese font(s): {', '.join(hits)}. "
            f"Paper figures must use English labels only — a CJK font in "
            f"the figure indicates Japanese axis/legend/annotation text. "
            f"Re-label in English and re-save."
        ]
    return []


def _check_no_in_figure_title(fig) -> list[str]:
    """Lab rule: the title goes in the LaTeX \\caption{} (or, for a talk
    slide, the beamer frametitle / slide caption), NEVER inside the figure
    -- on ANY figure, paper OR slide.

    Returns a list of violations (axis index + title text).  An empty
    list means clean.  emit_paper_figure() escalates this to a raise on
    `on_fail='raise'` since shipping a figure with an in-figure title is a
    style break.  The ``beamer_169_*`` slide profiles are NOT exempt: a
    slide's title is the beamer frametitle, not baked into the PNG.
    """
    bad = []
    # Per-axes titles
    for i, ax in enumerate(fig.get_axes()):
        ttl = ax.get_title()
        if ttl and ttl.strip():
            bad.append(f"axis {i}: ax.set_title({ttl!r})")
    # Figure-level suptitle
    sup = getattr(fig, "_suptitle", None)
    if sup is not None and sup.get_text().strip():
        bad.append(f"figure: fig.suptitle({sup.get_text()!r})")
    return bad


def _is_axis_reference_line(line, ax) -> bool:
    """Return whether *line* is an ``axhline``/``axvline`` reference line.

    Those helpers span the whole axes at a constant coordinate, so any legend
    placed inside the axes necessarily crosses them.  They are background
    furniture -- a zero line, a threshold, a spec limit -- not plotted data,
    and the lab rule is that a legend must not cover *data*.  The discriminator
    is the transform: ``ax.plot`` draws in ``ax.transData`` while ``axhline`` /
    ``axvline`` use the blended axis transforms, so an unlabelled but genuine
    data curve is still audited.
    """
    try:
        transform = line.get_transform()
        return (
            transform is ax.get_xaxis_transform()
            or transform is ax.get_yaxis_transform()
        )
    except Exception:
        return False


def _line_style_draws(linestyle) -> bool:
    """Return whether a Line2D linestyle renders a stroked path."""
    if linestyle is None:
        return False
    if isinstance(linestyle, tuple):        # (offset, onoffseq) dash spec
        offset, onoffseq = linestyle
        return bool(onoffseq) and any(float(x) > 0 for x in onoffseq)
    return str(linestyle).strip().lower() not in {"", "none"}


def _marker_draws(marker) -> bool:
    """Return whether a Line2D marker renders anything."""
    if marker is None:
        return False
    if isinstance(marker, str):
        return marker.strip().lower() not in {"", "none"}
    return True                             # Path / MarkerStyle / int style


def _check_legend_no_overlap(fig) -> list[dict]:
    """Lab rule: legends must not overlap plotted data.

    The check runs in display coordinates after the canvas is drawn.  It
    tests the rendered path of every ``Line2D`` (so a sparse two-point line
    crossing a legend is still caught), its vertices/markers, and collection
    offsets such as scatter points.  The legend bbox is padded by 1 pt to
    account for stroke and marker extent.
    """
    plt, _ = _get_mpl()
    import numpy as np
    from matplotlib.transforms import Bbox

    fig.canvas.draw()
    out = []
    for i, ax in enumerate(fig.get_axes()):
        legend = ax.get_legend()
        if legend is None:
            continue
        try:
            leg_bbox = legend.get_window_extent()
        except RuntimeError:
            continue
        pad_px = fig.dpi / 72.0
        hit_bbox = Bbox.from_extents(
            leg_bbox.x0 - pad_px,
            leg_bbox.y0 - pad_px,
            leg_bbox.x1 + pad_px,
            leg_bbox.y1 + pad_px,
        )
        for line_index, line in enumerate(ax.lines):
            if not line.get_visible():
                continue
            if _is_axis_reference_line(line, ax):
                continue
            # Only artists that actually render can overlap the legend.  A
            # marker-only Line2D (linestyle='None', e.g. ax.plot(x, y, 'o'))
            # still carries connecting segments in its path, but none of them
            # are stroked, so testing those segments reports an overlap the
            # reader can never see.
            draws_stroke = (
                _line_style_draws(line.get_linestyle())
                and float(line.get_linewidth()) > 0.0
            )
            draws_marker = _marker_draws(line.get_marker())
            if not draws_stroke and not draws_marker:
                continue
            label = line.get_label()
            if not label or label.startswith("_"):
                label = f"Line2D[{line_index}]"
            path_disp = line.get_path().transformed(line.get_transform())
            pts_disp = np.asarray(path_disp.vertices, dtype=float)
            finite = np.isfinite(pts_disp).all(axis=1)
            pts_disp = pts_disp[finite]
            inside = (
                (pts_disp[:, 0] >= hit_bbox.x0)
                & (pts_disp[:, 0] <= hit_bbox.x1)
                & (pts_disp[:, 1] >= hit_bbox.y0)
                & (pts_disp[:, 1] <= hit_bbox.y1)
            )
            # Vertices lie on the rendered artist either way: they are the
            # marker positions for a marker-only line and points on the curve
            # for a stroked one.
            n_over = int(inside.sum())
            segment_hit = bool(
                draws_stroke
                and path_disp.intersects_bbox(hit_bbox, filled=False)
            )
            if n_over > 0 or segment_hit:
                out.append({
                    "axis": i,
                    "label": label,
                    "artist": "Line2D",
                    "n_overlap_points": n_over,
                    "segment_intersection": segment_hit,
                    "fraction": n_over / max(pts_disp.shape[0], 1),
                })

        for collection_index, collection in enumerate(ax.collections):
            if not collection.get_visible():
                continue
            label = collection.get_label()
            if not label or label.startswith("_"):
                label = f"{type(collection).__name__}[{collection_index}]"
            get_offsets = getattr(collection, "get_offsets", None)
            get_offset_transform = getattr(collection, "get_offset_transform", None)
            if get_offsets is None or get_offset_transform is None:
                continue
            offsets = np.asarray(get_offsets(), dtype=float)
            if offsets.size == 0:
                continue
            offsets = np.atleast_2d(offsets)
            pts_disp = get_offset_transform().transform(offsets)
            finite = np.isfinite(pts_disp).all(axis=1)
            pts_disp = pts_disp[finite]
            inside = (
                (pts_disp[:, 0] >= hit_bbox.x0)
                & (pts_disp[:, 0] <= hit_bbox.x1)
                & (pts_disp[:, 1] >= hit_bbox.y0)
                & (pts_disp[:, 1] <= hit_bbox.y1)
            )
            n_over = int(inside.sum())
            if n_over > 0:
                out.append({
                    "axis": i,
                    "label": label,
                    "artist": type(collection).__name__,
                    "n_overlap_points": n_over,
                    "segment_intersection": False,
                    "fraction": n_over / max(pts_disp.shape[0], 1),
                })
    return out


def emit_paper_figure(
    fig,
    path,
    profile,
    min_axes_fraction: float = 0.72,
    on_fail: Literal["raise", "warn", "auto_tighten"] = "raise",
    check_title_in_figure: bool = True,
    check_legend_overlap: bool = True,
    check_colorblind_safe: bool = True,
    check_font_embedding: bool = True,
    check_no_japanese: bool = True,
    check_times_new_roman: bool = True,
    check_min_font_size: bool = True,
    min_font_pt: Optional[float] = None,
    max_font_pt: Optional[float] = None,
    embed_width_cm: Optional[float] = None,
    check_text_overflow: bool = True,
    text_overflow_tol_pt: float = 0.5,
    dpi: int = 600,
    save_pdf: bool = True,
    save_png: bool = True,
    verbose: bool = True,
) -> dict:
    """Save a figure at the journal's exact width with an efficiency gate.

    The pipeline:
      1. Ensure the figure's width matches the profile (within 1 mm).
         If not, resize the figure (preserving aspect).
      2. Measure efficiency.
      3. If efficiency >= min_axes_fraction: save.
      4. Else, behave per `on_fail`:
         'raise' (default):
            ValueError with a per-margin suggestion of WHICH margin is
            the biggest waste.
         'warn':
            warnings.warn() but proceed to save.
         'auto_tighten':
            run auto_tighten() once; re-measure; if still failing
            switch to 'warn' and save.

    Args:
        fig: matplotlib Figure.
        path: Output path WITHOUT extension (we add .pdf and .png).
        profile: PaperProfile or string key.
        min_axes_fraction: gate floor.  0.72 is the empirical default
            (achievable for most layouts AFTER auto_tighten); raise to
            0.78-0.82 for ultra-tight, lower to 0.65 for first drafts.
            Note: 8 pt IEEE x/y labels + tick labels physically require
            ~13-18% of figure for axis frame, so 0.80+ is reachable
            only for layouts with shared axes (e.g. 1x2 with sharey).
        on_fail: gate behaviour (see above).
        check_no_japanese: when True (default), reject figures that
            contain Japanese (CJK) text in any artist (axis labels,
            ticks, legend, annotations, title) AND reject saved PDFs
            that embed a Japanese font.  Paper figures must use English
            labels only.  Set False only for a figure deliberately
            targeting a Japanese-language venue.
        check_times_new_roman: when True (default), require matplotlib
            rcParams to request Times New Roman and reject saved PDFs
            that contain common serif fallback fonts.
        check_min_font_size: when True (default), inspect every visible text
            artist after scaling to the final paper/slide width.
        min_font_pt: optional final-size font floor.  When omitted, paper
            profiles require 10 pt and 16:9 slide profiles require 20 pt.
            Slide profiles still author at 24 pt by default to leave scaling
            margin before the 20 pt final-size gate.
        max_font_pt: optional final-size font ceiling.  When omitted, ordinary
            paper profiles use 10.5 pt, keeping all plot text at the 10 pt body
            target without making a figure look typographically oversized.
            Slide profiles have no default ceiling.
        embed_width_cm: actual width at which the output will be embedded in
            the paper or pasted into PowerPoint.  Defaults to the profile
            width (100% scale).  Passing the real pasted width is mandatory
            when the slide layout scales the generated figure.
        check_text_overflow: when True (default), reject free diagram
            labels created with ``ax.text`` / ``fig.text`` if their
            rendered bbox extends past the fixed figure canvas.  Axis
            labels, tick labels, legends, and titles are handled by
            separate gates.
        dpi: PNG/JPEG resolution.  600 for line-art with thin strokes,
            300 for general use, 1200 for IEEE camera-ready.
        save_pdf, save_png: which formats to emit.
        verbose: print pre/post measurement to stdout.

    Returns:
        Same dict shape as measure_figure_efficiency, plus:
          - 'wrote': list of paths written
          - 'gate_passed': bool
          - 'auto_tightened': bool
    """
    plt, _ = _get_mpl()
    prof = get_profile(profile)
    from pathlib import Path
    p = Path(path)
    if p.suffix:
        # User passed an extension — drop it, we'll add back
        p = p.with_suffix("")

    # --- Pre-flight (1): no titles inside the figure ---
    # Lab rule: every title goes in the LaTeX \caption{}, never in the
    # figure itself.  This check is FIRST because if it fires the user
    # has a 1-line fix (delete the set_title call) before we even start
    # measuring efficiency.
    # No in-figure title on ANY figure -- paper OR slide.  The title goes
    # in the LaTeX \caption / beamer frametitle, never inside the plot
    # (the beamer_169_* slide profiles are NOT exempt).
    if check_title_in_figure:
        title_violations = _check_no_in_figure_title(fig)
        if title_violations:
            msg = (
                "emit_paper_figure: titles detected inside the figure.\n"
                "Lab style: titles go in the LaTeX \\caption{...}, NOT in\n"
                "the figure.  Remove the following calls and write the\n"
                "title in the caption:\n  "
                + "\n  ".join(title_violations)
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                # auto_tighten cannot fix a title, so escalate to raise.
                raise ValueError(msg)

    # --- Pre-flight (1b): colorblind-safe palette ---
    # Lab rule: every artist's color must be either greyscale or in
    # the Okabe-Ito 8-color palette (Wong 2011, Nature Methods).
    # Matplotlib's default tab10 cycle confuses deuteranopia/protanopia
    # readers (orange vs red); we set Okabe-Ito by default in
    # paper_figure(), but a user could override with raw plt.plot(...,
    # color='red') -- this check catches that.
    if check_colorblind_safe:
        cvd_violations = _check_colors_are_cvd_safe(fig)
        if cvd_violations:
            lines = [
                f"  axis {v['axis']}: '{v['label']}' uses {v['color']}: "
                f"{v['reason']}"
                for v in cvd_violations
            ]
            msg = (
                "emit_paper_figure: non-CVD-safe colors detected.\n"
                + "\n".join(lines) + "\n"
                "Fix:\n"
                "  - use the profile's default color cycle (don't pass\n"
                "    color=... explicitly), OR\n"
                "  - pick from Okabe-Ito: black, orange (#E69F00), sky-blue\n"
                "    (#56B4E9), bluish-green (#009E73), yellow (#F0E442),\n"
                "    blue (#0072B2), vermillion (#D55E00), reddish-purple\n"
                "    (#CC79A7), OR\n"
                "  - any pure greyscale (r==g==b) -- always CVD-safe."
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                raise ValueError(msg)

    # --- Pre-flight (1c): no Japanese text in the figure ---
    # Lab rule: paper figures use ENGLISH labels only.  Japanese text
    # forces a CJK font into the PDF, which English-language publishers
    # (IEEE / Elsevier) reject in pre-flight and which bloats the file.
    # Caught BEFORE savefig so the fix is a 1-line re-label.
    if check_no_japanese:
        jp_violations = _check_no_japanese_text(fig)
        if jp_violations:
            lines = [
                f"  axis {v['axis']} {v['kind']}: {v['text']!r}"
                for v in jp_violations
            ]
            msg = (
                "emit_paper_figure: Japanese (CJK) text detected in the "
                "figure.\n"
                + "\n".join(lines) + "\n"
                "Lab rule: paper figures carry ENGLISH labels only.\n"
                "Fix:\n"
                "  - Re-label axes / legend / annotations in English.\n"
                "  - Keep the Japanese version (if any) as a SEPARATE\n"
                "    figure for Japanese-language venues; the IEEE / IEEJ-\n"
                "    English / IGTE / Compumag figure must be English.\n"
                "Why: Japanese text embeds a CJK font (MS Gothic / Yu\n"
                "Gothic / IPAex / Noto CJK) in the PDF — publishers reject\n"
                "it in font pre-flight and it bloats the file."
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                # auto_tighten cannot fix Japanese text, so escalate.
                raise ValueError(msg)

    # --- Pre-flight (1d): Times New Roman is the lab standard ---
    if check_times_new_roman:
        _assert_times_new_roman_available()
        if not _rcparams_request_times_new_roman():
            msg = (
                "emit_paper_figure: Times New Roman is not requested in "
                "matplotlib rcParams.\n"
                "Sugahara Lab standard: publication figures are generated "
                "with Times New Roman.  Build the figure with "
                "paper_figure(..., use_times_roman=True) or set:\n"
                "  plt.rcParams['font.family'] = 'serif'\n"
                "  plt.rcParams['font.serif'] = ['Times New Roman']"
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                raise ValueError(msg)

    # --- Pre-flight (1e): audit every Text artist at FINAL embed scale ---
    # rcParams only describe authored text.  PowerPoint / LaTeX may scale the
    # generated figure, so a 24 pt source label can become <20 pt on a slide.
    # Paper figures similarly require every item to remain >=10 pt at the
    # actual printed width, including ticks, legends, and annotations.
    effective_min_font_pt = min_font_pt
    is_slide_profile = prof.name.startswith("beamer_169_")
    if effective_min_font_pt is None:
        effective_min_font_pt = 20.0 if is_slide_profile else 10.0
    effective_max_font_pt = max_font_pt
    if effective_max_font_pt is None and not is_slide_profile:
        effective_max_font_pt = 10.5
    intended_embed_width_cm = (
        float(embed_width_cm)
        if embed_width_cm is not None
        else float(prof.width_mm) / 10.0
    )
    if intended_embed_width_cm <= 0:
        raise ValueError("embed_width_cm must be positive.")
    authored_width_cm = float(prof.width_mm) / 10.0
    embed_scale = intended_embed_width_cm / authored_width_cm
    font_size_violations = []
    if check_min_font_size and effective_min_font_pt is not None:
        from .tools import check_min_font
        font_size_violations = check_min_font(
            fig,
            min_pt=float(effective_min_font_pt),
            max_pt=(
                None
                if effective_max_font_pt is None
                else float(effective_max_font_pt)
            ),
            embed_width_cm=intended_embed_width_cm,
            authored_width_cm=authored_width_cm,
        )
        if font_size_violations:
            lines = [
                f"  {v['kind']} {v['text']!r}: {v['render_pt']:.1f} pt "
                f"source -> {v['visible_pt']:.1f} pt displayed "
                f"({v['reason']})"
                for v in font_size_violations
            ]
            has_above = any(
                v["reason"] == "above maximum" for v in font_size_violations
            )
            msg = (
                "emit_paper_figure: visible text outside the final-size font "
                "band.\n" if has_above else
                "emit_paper_figure: visible text below the final-size font floor.\n"
            ) + (
                "\n".join(lines) + "\n"
                f"Source width {authored_width_cm:g} cm -> final width "
                f"{intended_embed_width_cm:g} cm ({embed_scale:.3f}x).\n"
                f"Allowed final-size range: {float(effective_min_font_pt):g} pt"
                + (
                    " or larger"
                    if effective_max_font_pt is None
                    else f" to {float(effective_max_font_pt):g} pt"
                )
                + ".\nAdjust source font sizes or the final embed width; "
                "simplify crowded content rather than shrinking or enlarging "
                "text away from the paper body size."
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                raise ValueError(msg)

    # --- Pre-flight (1d): free diagram text must fit the canvas ---
    # This catches concept diagrams where ax.text() / annotate() labels
    # are clipped by tight margins.  It deliberately excludes axis
    # labels, tick labels, titles, and legends; those have dedicated
    # checks elsewhere.
    if check_text_overflow:
        text_overflow = _check_text_overflow(fig, tol_pt=text_overflow_tol_pt)
        if text_overflow:
            lines = [
                f"  axis {v['axis']} {v['kind']} {v['text']!r} "
                f"extends {v['overhang_pt']:.1f} pt past the {v['side']} edge"
                for v in text_overflow
            ]
            msg = (
                "emit_paper_figure: diagram text extends past the figure "
                "canvas.\n"
                + "\n".join(lines) + "\n"
                "Fix:\n"
                "  - move the text inside the axes / figure,\n"
                "  - shorten the label,\n"
                "  - increase aspect / width, or\n"
                "  - reserve margin with fig.subplots_adjust(...)."
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                # auto_tighten cannot know which diagram label to move.
                raise ValueError(msg)

    # --- Pre-flight (2): no legend overlapping data lines ---
    # Lab rule: legends must not overlap data.  An overlapping legend
    # is the single most common reviewer-visible obvious flaw.
    if check_legend_overlap:
        overlaps = _check_legend_no_overlap(fig)
        if overlaps:
            lines = [
                f"  axis {o['axis']}: legend overlaps '{o['label']}' "
                + (
                    f"({o['n_overlap_points']} points = "
                    f"{o['fraction']*100:.0f}% of curve)"
                    if o["n_overlap_points"]
                    else "(a drawn segment crosses the legend box; no "
                         "sampled point lands inside it)"
                )
                for o in overlaps
            ]
            msg = (
                "emit_paper_figure: legend overlaps data line(s):\n"
                + "\n".join(lines) + "\n"
                "Fix:\n"
                "  - call find_best_legend_loc(ax) before ax.legend(loc=...)\n"
                "  - OR use label_curve_endpoints(ax, [...]) for direct\n"
                "    inline labels in the right margin (lab preferred)"
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                raise ValueError(msg)

    # --- Pre-flight (3): enforce width ---
    fw, fh = fig.get_size_inches()
    target_w_in = prof.width_in
    if abs(fw - target_w_in) > (1.0 / 25.4):     # > 1 mm difference
        # Preserve aspect ratio when resizing.
        aspect = fh / fw if fw > 0 else 1.0
        fig.set_size_inches(target_w_in, target_w_in * aspect)
        if verbose:
            print(f"[emit_paper_figure] Resized fig from "
                  f"{fw:.2f}x{fh:.2f} in to "
                  f"{target_w_in:.2f}x{target_w_in * aspect:.2f} in "
                  f"(target = {prof.full_name})")

    # --- Measure ---
    pre = measure_figure_efficiency(fig)
    if verbose:
        print(f"[emit_paper_figure] {prof.full_name}: efficiency "
              f"{pre['axes_area_fraction']:.3f} (target >= "
              f"{min_axes_fraction:.2f})")

    auto_tightened = False
    if pre["axes_area_fraction"] < min_axes_fraction:
        if on_fail == "auto_tighten":
            if verbose:
                print(f"[emit_paper_figure] auto_tighten -> target "
                      f"{min_axes_fraction:.2f}")
            auto_tighten(fig, target_axes_fraction=min_axes_fraction)
            auto_tightened = True
            post = measure_figure_efficiency(fig)
            if post["axes_area_fraction"] < min_axes_fraction:
                import warnings
                warnings.warn(
                    f"emit_paper_figure: still {post['axes_area_fraction']:.3f}"
                    f" after auto_tighten (target {min_axes_fraction:.2f}). "
                    f"Saving anyway because on_fail='auto_tighten' degrades "
                    f"to warn after one attempt."
                )
        elif on_fail == "warn":
            import warnings
            warnings.warn(_efficiency_hint_message(pre, min_axes_fraction))
        elif on_fail == "raise":
            raise ValueError(_efficiency_hint_message(pre, min_axes_fraction))
        else:
            raise ValueError(
                f"on_fail must be 'raise', 'warn', or 'auto_tighten'; "
                f"got {on_fail!r}"
            )

    # --- Save ---
    wrote = []
    if save_pdf:
        pdf_path = p.with_suffix(".pdf")
        fig.savefig(pdf_path, dpi=dpi, bbox_inches=None)   # exact size, no tight
        wrote.append(str(pdf_path))
    if save_png:
        png_path = p.with_suffix(".png")
        fig.savefig(png_path, dpi=dpi, bbox_inches=None)
        wrote.append(str(png_path))

    # --- Post-save: verify PDF font embedding ---
    # paper_figure() sets rcParams['pdf.fonttype']=42, but a user could
    # override that between paper_figure() and emit_paper_figure().  The
    # only way to know whether the FINAL pdf has Type-42 (TrueType) is
    # to inspect the file bytes.  Cheap, runs after save, surfaces the
    # IEEE/Elsevier pre-flight reject reason BEFORE submission.
    font_violations = []
    if check_font_embedding and save_pdf:
        pdf_path = p.with_suffix(".pdf")
        font_violations = _check_pdf_fonts_embedded(pdf_path)
        if font_violations:
            msg = (
                "emit_paper_figure: saved PDF has font-embedding "
                "issues that publishers reject:\n  "
                + "\n  ".join(font_violations) + "\n"
                "Fix:\n"
                "  - Ensure plt.rcParams['pdf.fonttype'] == 42 BEFORE\n"
                "    savefig().  paper_figure() sets this; do not\n"
                "    override.\n"
                "  - If you must use a custom font, use a TrueType\n"
                "    (.ttf / .otf) family installed on the system."
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                # File is already written; raise AFTER so the user
                # can still inspect what was produced.
                raise ValueError(msg)

    # --- Post-save: verify Times New Roman did not silently fall back ---
    times_new_roman_violations = []
    if check_times_new_roman and save_pdf:
        pdf_path = p.with_suffix(".pdf")
        times_new_roman_violations = _check_pdf_times_new_roman(pdf_path)
        if times_new_roman_violations:
            msg = (
                "emit_paper_figure: saved PDF appears to use a serif "
                "fallback instead of Times New Roman:\n  "
                + "\n  ".join(times_new_roman_violations) + "\n"
                "Fix: rebuild the figure with paper_figure() after "
                "installing Times New Roman / rebuilding the matplotlib "
                "font cache."
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                raise ValueError(msg)

    # --- Post-save: verify NO Japanese font embedded in the PDF ---
    # The live-artist scan above (_check_no_japanese_text) catches the
    # common case; this confirms the FINAL bytes also carry no CJK font
    # (e.g. caught a mathtext fallback or a stray glyph the artist scan
    # missed).
    japanese_font_violations = []
    if check_no_japanese and save_pdf:
        pdf_path = p.with_suffix(".pdf")
        japanese_font_violations = _check_no_japanese_font_in_pdf(pdf_path)
        if japanese_font_violations:
            msg = (
                "emit_paper_figure: saved PDF embeds a Japanese font:\n  "
                + "\n  ".join(japanese_font_violations) + "\n"
                "Fix: re-label all figure text in English (the figure\n"
                "must contain no Japanese characters), then re-save."
            )
            if on_fail == "warn":
                import warnings
                warnings.warn(msg)
            elif on_fail in ("raise", "auto_tighten"):
                raise ValueError(msg)

    final = measure_figure_efficiency(fig)
    final["wrote"] = wrote
    final["gate_passed"] = final["axes_area_fraction"] >= min_axes_fraction
    final["auto_tightened"] = auto_tightened
    final["profile"] = prof.name
    final["font_violations"] = font_violations
    final["font_size_violations"] = font_size_violations
    final["min_font_pt"] = effective_min_font_pt
    final["max_font_pt"] = effective_max_font_pt
    final["authored_width_cm"] = authored_width_cm
    final["embed_width_cm"] = intended_embed_width_cm
    final["embed_scale"] = embed_scale
    final["times_new_roman_violations"] = times_new_roman_violations
    final["japanese_font_violations"] = japanese_font_violations

    if verbose:
        print(f"[emit_paper_figure] wrote {len(wrote)} file(s): "
              f"{', '.join(wrote)}")
        print(f"[emit_paper_figure] final efficiency: "
              f"{final['axes_area_fraction']:.3f}  "
              f"gate_passed={final['gate_passed']}")
        if font_violations:
            print(f"[emit_paper_figure] font issues: {len(font_violations)}")
    return final


def _efficiency_hint_message(measurement: dict, target: float) -> str:
    """Construct the suggestion message for the gate."""
    margin = measurement["margin_breakdown"]
    # Sort margins by how much room they have, descending
    sorted_m = sorted(margin.items(), key=lambda kv: -kv[1])
    suggestion_lines = ["Suggested actions (largest waste first):"]
    for side, frac in sorted_m:
        if frac > 0.10:
            suggestion_lines.append(
                f"  - {side:6s} margin = {frac:.3f} of figure; "
                f"consider fig.subplots_adjust({side}="
                f"{(frac - 0.04 if side in ('left', 'bottom') else 1.0 - frac + 0.04):.3f}"
                f")  [shave ~0.04]"
            )
    wspace = measurement.get("wspace_estimated")
    hspace = measurement.get("hspace_estimated")
    if wspace and wspace > 0.18:
        suggestion_lines.append(
            f"  - wspace ~ {wspace:.2f}; "
            f"fig.subplots_adjust(wspace={max(wspace - 0.05, 0.10):.2f}) "
            f"likely tightens 1x2/2x2 layouts safely"
        )
    if hspace and hspace > 0.25:
        suggestion_lines.append(
            f"  - hspace ~ {hspace:.2f}; "
            f"fig.subplots_adjust(hspace={max(hspace - 0.05, 0.18):.2f}) "
            f"likely tightens 2x1/2x2 layouts safely"
        )
    suggestion_lines.append(
        "  - OR call auto_tighten(fig, target_axes_fraction="
        f"{target:.2f}) to apply the above iteratively."
    )
    return (
        f"emit_paper_figure: axes_area_fraction = "
        f"{measurement['axes_area_fraction']:.3f} < "
        f"{target:.2f} (paper-quality floor).\n"
        + "\n".join(suggestion_lines)
    )
