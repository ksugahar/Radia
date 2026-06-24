"""Golden tests for radia_mcp.figure.paper_figure (v0.78.0).

Lock in:
  - Per-profile baseline efficiency (no profile drift)
  - auto_tighten gains >= 5 efficiency points on a representative layout
  - auto_tighten does NOT introduce clipping past the per-side
    tolerance vs. baseline
  - emit_paper_figure() gate behaviour: raise / warn / auto_tighten
  - Output files actually appear on disk at the correct dimensions

These tests REQUIRE matplotlib at runtime.  If matplotlib is not
available they are skipped (not failed); the MCP server tools above
remain testable without matplotlib via test_each_server_selftest.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest


# All tests need matplotlib; skip the whole module if absent.
matplotlib = pytest.importorskip("matplotlib")
import matplotlib  # noqa: E402
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.font_manager as fm  # noqa: E402


def _tnr_available():
    try:
        p = fm.findfont("Times New Roman", fallback_to_default=False)
        return "times" in Path(str(p)).name.lower()
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not _tnr_available(),
    reason="Times New Roman is required by the Sugahara Lab figure standard",
)

from radia_mcp.figure import (  # noqa: E402
    PROFILES,
    PaperProfile,
    OKABE_ITO,
    paper_figure,
    measure_figure_efficiency,
    auto_tighten,
    emit_paper_figure,
    add_panel_labels,
    audit_text_overflow,
)


# --------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------


def _populate(axes):
    """Put representative data + labels in every panel."""
    import numpy as np
    for i, ax in enumerate(axes.flat):
        x = np.linspace(0, 10, 50)
        y = x ** 1.5 + 2 * i
        ax.plot(x, y, label=f"data {i}")
        ax.set_xlabel(r"$f$ (Hz)")
        ax.set_ylabel(r"$|Z|$ ($\Omega$)")


# --------------------------------------------------------------------
# Baseline efficiency lock per profile
# --------------------------------------------------------------------


@pytest.mark.parametrize("prof_name", sorted(PROFILES.keys()))
def test_baseline_efficiency_in_band(prof_name):
    """Each profile's 1x1 baseline must land in a sensible 0.50..0.90 band.

    Catches profile drift where someone tweaks margins so loose that the
    default falls below 0.50 (waste) or so tight that there is no room
    for the user's labels (above 0.90).  Either extreme indicates the
    profile is broken.
    """
    fig, axes = paper_figure(prof_name, nrows=1, ncols=1)
    _populate(axes)
    m = measure_figure_efficiency(fig)
    plt.close(fig)
    eff = m["axes_area_fraction"]
    assert 0.50 <= eff <= 0.90, (
        f"profile {prof_name!r} baseline 1x1 efficiency {eff:.3f} "
        f"is outside the 0.50..0.90 sanity band — profile margins drifted"
    )


# --------------------------------------------------------------------
# auto_tighten convergence + safety
# --------------------------------------------------------------------


def test_auto_tighten_gains_at_least_5pt_ieee_1x2():
    """The canonical use-case (IEEE 2-column 1x2) gains >= 5 efficiency
    points after auto_tighten.

    If this drops, either the per-layout margin defaults already drove
    efficiency too high to gain (in which case the default deltas
    should be loosened so the gate hint is more visible), or
    auto_tighten regressed.  Either way: investigate.
    """
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=2)
    _populate(axes)
    pre = measure_figure_efficiency(fig)["axes_area_fraction"]
    res = auto_tighten(fig, target_axes_fraction=0.80, max_iter=24)
    plt.close(fig)
    gain_pts = (res["final_efficiency"] - pre) * 100
    assert gain_pts >= 5.0, (
        f"IEEE 1x2 auto_tighten gained only {gain_pts:.1f} pts "
        f"({pre:.3f} -> {res['final_efficiency']:.3f}); expected >= 5"
    )


def test_auto_tighten_no_new_clipping():
    """auto_tighten must not push any side overhang past baseline + tol.

    The per-side overhang tolerance is 2% of figure width / 3% of
    figure height (see _OVERHANG_TOL_* in _paper_figure.py).  After
    auto_tighten, the final overhang on each side must be <=
    baseline + tolerance.  Locks in the safety invariant.
    """
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=2)
    _populate(axes)
    from radia_mcp.figure._paper_figure import _compute_overhang_px
    baseline = _compute_overhang_px(fig)
    res = auto_tighten(fig, target_axes_fraction=0.85, max_iter=30)
    final_over = _compute_overhang_px(fig)
    fb_w = fig.bbox.width
    fb_h = fig.bbox.height
    tol_x = 0.02 * fb_w
    tol_y = 0.03 * fb_h
    plt.close(fig)
    for side, tol in [("left", tol_x), ("right", tol_x),
                       ("top", tol_y), ("bottom", tol_y)]:
        assert final_over[side] <= baseline[side] + tol + 0.5, (
            f"auto_tighten pushed {side} overhang past safety:\n"
            f"  baseline = {baseline[side]:.1f} px\n"
            f"  final    = {final_over[side]:.1f} px\n"
            f"  tol      = {tol:.1f} px (= {tol/(fb_w if side in ('left','right') else fb_h)*100:.1f}% of fig {side[0]}-dim)"
        )


def test_auto_tighten_converges_within_3_passes():
    """Repeated auto_tighten calls converge — each pass either makes
    real progress (the baseline-overhang is re-snapshotted so a 2nd
    pass can squeeze a few more pixels), or does nothing.  Pass 3
    must gain no more than 1 pt over pass 2.

    Without this lock, an auto_tighten regression that oscillates
    (alternately tightening then over-correcting) would silently leak
    into shipped figures.
    """
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=2)
    _populate(axes)
    auto_tighten(fig, target_axes_fraction=0.85, max_iter=24)
    eff1 = measure_figure_efficiency(fig)["axes_area_fraction"]
    auto_tighten(fig, target_axes_fraction=0.85, max_iter=24)
    eff2 = measure_figure_efficiency(fig)["axes_area_fraction"]
    auto_tighten(fig, target_axes_fraction=0.85, max_iter=24)
    eff3 = measure_figure_efficiency(fig)["axes_area_fraction"]
    plt.close(fig)
    # Monotonically non-decreasing
    assert eff2 >= eff1 - 0.005, (
        f"auto_tighten regressed pass 1->2: {eff1:.3f} -> {eff2:.3f}"
    )
    assert eff3 >= eff2 - 0.005, (
        f"auto_tighten regressed pass 2->3: {eff2:.3f} -> {eff3:.3f}"
    )
    # Converged: pass 3 gains <= 1 pt vs pass 2
    assert (eff3 - eff2) <= 0.01, (
        f"auto_tighten not converged after 3 passes: pass 2 -> "
        f"{eff2:.3f}, pass 3 -> {eff3:.3f} (delta {eff3-eff2:.3f})"
    )


# --------------------------------------------------------------------
# emit_paper_figure gate
# --------------------------------------------------------------------


def test_emit_paper_figure_raises_below_floor(tmp_path):
    """on_fail='raise' must raise ValueError when efficiency < floor."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=2)
    _populate(axes)
    out = tmp_path / "fig"
    with pytest.raises(ValueError, match=r"axes_area_fraction"):
        emit_paper_figure(
            fig, str(out), "ieee_double_column",
            min_axes_fraction=0.95,     # impossibly high — guaranteed fail
            on_fail="raise",
            verbose=False,
        )
    plt.close(fig)
    # Files MUST NOT exist (raise short-circuits the save).
    assert not out.with_suffix(".pdf").exists()
    assert not out.with_suffix(".png").exists()


def test_emit_paper_figure_auto_tighten_recovers(tmp_path):
    """on_fail='auto_tighten' should pull efficiency above the floor
    (or save anyway with a warn after one attempt)."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=2)
    _populate(axes)
    out = tmp_path / "fig_at"
    res = emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.72,
        on_fail="auto_tighten",
        verbose=False,
    )
    plt.close(fig)
    # Final efficiency should be at least equal to (probably > than) 0.72.
    assert res["axes_area_fraction"] >= 0.70, (
        f"auto_tighten gate left efficiency at "
        f"{res['axes_area_fraction']:.3f}"
    )
    assert out.with_suffix(".pdf").exists()
    assert out.with_suffix(".png").exists()


def test_emit_paper_figure_warn_saves(tmp_path):
    """on_fail='warn' must save even when below floor."""
    import warnings
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=2)
    _populate(axes)
    out = tmp_path / "fig_warn"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        emit_paper_figure(
            fig, str(out), "ieee_double_column",
            min_axes_fraction=0.99, on_fail="warn",
            verbose=False,
        )
    plt.close(fig)
    assert out.with_suffix(".pdf").exists()
    assert out.with_suffix(".png").exists()
    assert any("axes_area_fraction" in str(w.message) for w in caught), (
        f"no efficiency-warn surfaced; got {[str(w.message) for w in caught]}"
    )


def test_emit_paper_figure_resizes_to_profile_width(tmp_path):
    """emit_paper_figure must enforce the journal's exact width when the
    user passes a misshapen figure."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    _populate(axes)
    # Deliberately mess up the size — emit should fix it.
    fig.set_size_inches(3.0, 2.0)
    out = tmp_path / "fig_resized"
    emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.50,   # forgiving so we test resize, not gate
        on_fail="warn",
        verbose=False,
    )
    # After emit, fig must be at profile width (181 mm = 7.126 in)
    w_in, _ = fig.get_size_inches()
    plt.close(fig)
    target_w_in = 181.0 / 25.4
    assert abs(w_in - target_w_in) < (1.0 / 25.4), (
        f"emit_paper_figure did not resize to profile width; "
        f"got {w_in:.3f} in, expected {target_w_in:.3f} in"
    )


# --------------------------------------------------------------------
# add_panel_labels
# --------------------------------------------------------------------


def test_add_panel_labels_writes_one_per_axis():
    """(a), (b), (c) ... labels — one per axis, in the right order."""
    fig, axes = paper_figure("ieee_double_column", nrows=2, ncols=2)
    _populate(axes)
    texts = add_panel_labels(axes)
    plt.close(fig)
    assert len(texts) == 4
    expected = ["(a)", "(b)", "(c)", "(d)"]
    for t, exp in zip(texts, expected):
        assert t.get_text() == exp, (
            f"expected {exp}, got {t.get_text()}"
        )


def test_add_panel_labels_custom_labels():
    """Allow custom labels (e.g. roman numerals)."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=3)
    _populate(axes)
    texts = add_panel_labels(axes, labels=["(i)", "(ii)", "(iii)"])
    plt.close(fig)
    assert [t.get_text() for t in texts] == ["(i)", "(ii)", "(iii)"]


def test_add_panel_labels_count_mismatch_raises():
    """Wrong number of labels -> ValueError, not silent truncation."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=3)
    _populate(axes)
    with pytest.raises(ValueError, match=r"labels provided"):
        add_panel_labels(axes, labels=["only", "two"])
    plt.close(fig)


# --------------------------------------------------------------------
# paper_figure() layout invariants
# --------------------------------------------------------------------


@pytest.mark.parametrize("nrows,ncols", [(1, 1), (1, 2), (2, 1), (2, 2), (1, 3)])
def test_paper_figure_returns_2d_axes(nrows, ncols):
    """Even for 1x1, axes is a (1,1) ndarray (consistent indexing)."""
    fig, axes = paper_figure("ieee_double_column", nrows=nrows, ncols=ncols)
    plt.close(fig)
    assert axes.shape == (nrows, ncols), (
        f"axes shape {axes.shape} != ({nrows},{ncols})"
    )


def test_paper_figure_width_matches_profile():
    """The figure created by paper_figure() must be exactly the profile's
    width in inches (within 0.1 mm)."""
    for name, prof in PROFILES.items():
        fig, _ = paper_figure(name)
        w_in, _ = fig.get_size_inches()
        plt.close(fig)
        target_w_in = prof.width_mm / 25.4
        assert abs(w_in - target_w_in) < (0.1 / 25.4), (
            f"profile {name}: figure width {w_in:.3f} in vs target "
            f"{target_w_in:.3f} in"
        )


def test_paper_figure_rcparams_have_type42_font():
    """Verify the Type-42 font embed is set (IEEE / Elsevier requirement)."""
    fig, _ = paper_figure("ieee_double_column")
    plt.close(fig)
    assert plt.rcParams["pdf.fonttype"] == 42, (
        "pdf.fonttype must be 42 (TrueType) for IEEE-compliant figures; "
        f"got {plt.rcParams['pdf.fonttype']}"
    )
    assert plt.rcParams["ps.fonttype"] == 42


def test_paper_figure_rcparams_use_times_new_roman():
    """Sugahara Lab standard: paper_figure() generates TNR figures."""
    fig, _ = paper_figure("ieee_double_column")
    plt.close(fig)
    fam = plt.rcParams["font.family"]
    fam = [fam] if isinstance(fam, str) else list(fam)
    assert "serif" in [str(x).lower() for x in fam]
    assert plt.rcParams["font.serif"][0] == "Times New Roman"
    assert plt.rcParams["mathtext.fontset"] == "stix"


def test_emit_raises_when_times_new_roman_not_requested(tmp_path):
    """emit_paper_figure() must reject non-TNR rcParams."""
    old_family = plt.rcParams["font.family"]
    old_serif = list(plt.rcParams["font.serif"])
    old_math = plt.rcParams["mathtext.fontset"]
    try:
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = ["DejaVu Serif"]
        plt.rcParams["mathtext.fontset"] = "dejavuserif"
        fig, ax = plt.subplots(figsize=(3.5, 2.4), dpi=100)
        ax.plot([0, 1], [0, 1], color="0.4")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        with pytest.raises(ValueError, match=r"Times New Roman"):
            emit_paper_figure(
                fig, str(tmp_path / "bad_font"), "ieee_single_column",
                min_axes_fraction=0.20,
                on_fail="raise",
                verbose=False,
                check_colorblind_safe=False,
                check_legend_overlap=False,
            )
    finally:
        plt.rcParams["font.family"] = old_family
        plt.rcParams["font.serif"] = old_serif
        plt.rcParams["mathtext.fontset"] = old_math
        plt.close("all")


# --------------------------------------------------------------------
# LAB FONT RULE: 10 pt @ 8 cm, absolute (not relative to column width)
# --------------------------------------------------------------------


@pytest.mark.parametrize("prof_name", sorted(PROFILES.keys()))
def test_profile_uses_10pt_body_font(prof_name):
    """Sugahara Lab absolute font rule: 10 pt body, 10 pt legend, 9 pt
    ticks for every profile regardless of column width.

    Catches regressions where someone tries to shrink the font to fit
    more axes or to bump it up for wider columns.  The 10 pt is an
    ABSOLUTE physical-readability standard, not relative to figure
    size.
    """
    prof = PROFILES[prof_name]
    assert prof.font_pt == 10.0, (
        f"{prof_name}.font_pt = {prof.font_pt} != 10.0 (lab rule)"
    )
    assert prof.legend_pt == 10.0, (
        f"{prof_name}.legend_pt = {prof.legend_pt} != 10.0 "
        f"(legends carry data, not afterthoughts -- must match body)"
    )
    assert prof.tick_pt == 9.0, (
        f"{prof_name}.tick_pt = {prof.tick_pt} != 9.0 "
        f"(IEEE/IEEJ convention: 1 pt below body)"
    )


def test_paper_figure_rcparams_have_10pt_body_font():
    """When paper_figure() runs, plt.rcParams['font.size'] must be 10.0."""
    fig, _ = paper_figure("ieee_single_column")
    plt.close(fig)
    assert plt.rcParams["font.size"] == 10.0, (
        f"font.size = {plt.rcParams['font.size']} != 10.0 after "
        f"paper_figure('ieee_single_column')"
    )


# --------------------------------------------------------------------
# Title-in-figure prohibition (lab rule: titles go in LaTeX caption)
# --------------------------------------------------------------------


def test_emit_raises_on_ax_set_title(tmp_path):
    """ax.set_title('...') must trigger the gate with a clear message."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    _populate(axes)
    axes[0, 0].set_title("Impedance vs frequency")   # forbidden
    # pytest.raises's `match` is re.search WITHOUT re.DOTALL, so we
    # use a single-token pattern that survives newlines in the message.
    with pytest.raises(ValueError, match=r"titles detected inside"):
        emit_paper_figure(
            fig, str(tmp_path / "fig"), "ieee_double_column",
            on_fail="raise", verbose=False,
        )
    plt.close(fig)


def test_emit_raises_on_fig_suptitle(tmp_path):
    """fig.suptitle('...') is equally forbidden."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=2)
    _populate(axes)
    fig.suptitle("AC analysis sweep")               # forbidden
    with pytest.raises(ValueError, match=r"titles detected inside"):
        emit_paper_figure(
            fig, str(tmp_path / "fig"), "ieee_double_column",
            on_fail="raise", verbose=False,
        )
    plt.close(fig)


def test_emit_allows_empty_title(tmp_path):
    """An axis with no title (default) must pass the title check."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    _populate(axes)
    # No set_title -- default empty
    out = tmp_path / "fig_clean"
    emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.40,             # relax efficiency for unit test
        on_fail="raise", verbose=False,
        check_legend_overlap=False,         # no legend -> overlap N/A
    )
    plt.close(fig)
    assert out.with_suffix(".pdf").exists()


def test_emit_title_check_can_be_disabled(tmp_path):
    """`check_title_in_figure=False` lets the title through (slide deck use)."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    _populate(axes)
    axes[0, 0].set_title("Slide deck title is OK with override")
    out = tmp_path / "fig_title_ok"
    emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.40,
        on_fail="raise", verbose=False,
        check_title_in_figure=False,
        check_legend_overlap=False,
    )
    plt.close(fig)
    assert out.with_suffix(".pdf").exists()


# --------------------------------------------------------------------
# Diagram text overflow gate
# --------------------------------------------------------------------


def test_audit_text_overflow_detects_ax_text_outside_canvas():
    """Free diagram labels made with ax.text() must be lintable."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    ax = axes[0, 0]
    ax.text(-0.20, 0.5, "outside label", transform=ax.transAxes)
    rep = audit_text_overflow(fig)
    plt.close(fig)
    assert rep, "expected free-text overflow to be reported"
    assert rep[0]["kind"] == "text"
    assert rep[0]["side"] == "left"
    assert rep[0]["overhang_pt"] > 0


def test_emit_raises_on_diagram_text_overflow(tmp_path):
    """emit_paper_figure must not ship clipped diagram labels."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    ax = axes[0, 0]
    ax.plot([0, 1], [0, 1])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.text(-0.20, 0.5, "outside label", transform=ax.transAxes)
    with pytest.raises(ValueError, match=r"diagram text extends"):
        emit_paper_figure(
            fig, str(tmp_path / "fig"), "ieee_double_column",
            min_axes_fraction=0.40,
            on_fail="raise",
            verbose=False,
            check_legend_overlap=False,
        )
    plt.close(fig)


def test_emit_text_overflow_check_can_be_disabled(tmp_path):
    """Known intentional outside text can be allowed explicitly."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    ax = axes[0, 0]
    ax.plot([0, 1], [0, 1])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.text(-0.20, 0.5, "outside label", transform=ax.transAxes)
    out = tmp_path / "fig_text_overflow_off"
    emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.40,
        on_fail="raise",
        verbose=False,
        check_legend_overlap=False,
        check_text_overflow=False,
    )
    plt.close(fig)
    assert out.with_suffix(".pdf").exists()


# --------------------------------------------------------------------
# Legend-overlap prohibition
# --------------------------------------------------------------------


def test_emit_raises_on_legend_overlap(tmp_path):
    """A legend placed ON TOP of data must trigger the gate.

    Use a strictly-ascending line so the upper-right corner of the
    axes (matplotlib's default `loc='upper right'`) is guaranteed to
    contain data points.  A sin wave can have its data dip away from
    the upper-right corner exactly when the legend lands there, which
    would falsely PASS the check.
    """
    import numpy as np
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    ax = axes[0, 0]
    # Strictly ascending line: y = x; upper-right corner is where the
    # data lives, so loc='upper right' guarantees overlap.
    x = np.linspace(0, 10, 200)
    y = x.copy()
    ax.plot(x, y, label="ascending line")
    ax.set_xlabel(r"$f$ (Hz)")
    ax.set_ylabel(r"$|Z|$ ($\Omega$)")
    ax.legend(loc="upper right")
    with pytest.raises(ValueError, match=r"legend overlaps data"):
        emit_paper_figure(
            fig, str(tmp_path / "fig"), "ieee_double_column",
            on_fail="raise", verbose=False,
        )
    plt.close(fig)


def test_emit_passes_when_legend_in_safe_corner(tmp_path):
    """A legend pushed into an empty corner must NOT trigger the gate."""
    import numpy as np
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    ax = axes[0, 0]
    # Monotone-increasing curve from bottom-left to top-right --
    # bottom-right corner is empty.
    x = np.linspace(0, 10, 200)
    y = x ** 1.5
    ax.plot(x, y, label="data")
    ax.set_xlabel(r"$f$ (Hz)")
    ax.set_ylabel(r"$|Z|$ ($\Omega$)")
    ax.legend(loc="upper left", frameon=False)
    out = tmp_path / "fig_legend_ok"
    emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.40,
        on_fail="raise", verbose=False,
    )
    plt.close(fig)
    assert out.with_suffix(".pdf").exists()


def test_emit_legend_check_can_be_disabled(tmp_path):
    """`check_legend_overlap=False` lets known overlaps through."""
    import numpy as np
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    ax = axes[0, 0]
    x = np.linspace(0, 10, 200)
    ax.plot(x, 5 + 4 * np.sin(x), label="data")
    ax.set_xlabel(r"$f$ (Hz)")
    ax.set_ylabel(r"$|Z|$ ($\Omega$)")
    ax.legend(loc="upper right")
    out = tmp_path / "fig_overlap_ok"
    emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.40,
        on_fail="raise", verbose=False,
        check_legend_overlap=False,
    )
    plt.close(fig)
    assert out.with_suffix(".pdf").exists()


# --------------------------------------------------------------------
# v0.80.0: CVD-safe color palette (Okabe-Ito / Wong 2011)
# --------------------------------------------------------------------


def test_default_color_cycle_is_okabe_ito():
    """Every profile must default to the Okabe-Ito CVD-safe cycle.

    Wong 2011 (Nature Methods 8:441) is the de facto colorblind-safe
    standard in scientific publishing.  Matplotlib's default `tab10`
    confuses deuteranopia readers (orange vs red).  This locks the
    lab default so no profile drifts to tab10.
    """
    from radia_mcp.figure import OKABE_ITO
    for name, prof in PROFILES.items():
        assert tuple(prof.color_cycle) == tuple(OKABE_ITO), (
            f"profile {name}: color_cycle != Okabe-Ito"
        )


def test_paper_figure_sets_okabe_ito_rcparams():
    """After paper_figure(), plt.rcParams['axes.prop_cycle'] yields
    Okabe-Ito colors in order (no tab10 leak)."""
    from radia_mcp.figure import OKABE_ITO
    fig, _ = paper_figure("ieee_single_column")
    plt.close(fig)
    actual = plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    # Compare first three (full equality breaks if matplotlib normalises
    # hex to lowercase or vice-versa)
    actual_low = [c.lower() for c in actual[:3]]
    expected_low = [c.lower() for c in OKABE_ITO[:3]]
    assert actual_low == expected_low, (
        f"rcParams cycle = {actual[:3]} != Okabe-Ito {OKABE_ITO[:3]}"
    )


def test_emit_raises_on_non_cvd_safe_color(tmp_path):
    """A line drawn with `color='red'` (NOT in Okabe-Ito) triggers
    the CVD-safe gate."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    ax = axes[0, 0]
    ax.plot([0, 1, 2], [0, 1, 4], color="red", label="bad")   # forbidden
    ax.set_xlabel(r"$f$ (Hz)")
    ax.set_ylabel(r"$|Z|$ ($\Omega$)")
    with pytest.raises(ValueError, match=r"non-CVD-safe colors"):
        emit_paper_figure(
            fig, str(tmp_path / "fig"), "ieee_double_column",
            on_fail="raise", verbose=False,
            check_legend_overlap=False,  # legend not present
        )
    plt.close(fig)


def test_emit_greyscale_passes(tmp_path):
    """Pure greyscale colors are always CVD-safe."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    ax = axes[0, 0]
    ax.plot([0, 1, 2], [0, 1, 4], color="0.5", label="grey")   # OK
    ax.set_xlabel(r"$f$ (Hz)")
    ax.set_ylabel(r"$|Z|$ ($\Omega$)")
    out = tmp_path / "fig_grey"
    emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.40,
        on_fail="raise", verbose=False,
        check_legend_overlap=False,
    )
    plt.close(fig)
    assert out.with_suffix(".pdf").exists()


def test_emit_cvd_check_can_be_disabled(tmp_path):
    """`check_colorblind_safe=False` lets non-safe colors through."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    ax = axes[0, 0]
    ax.plot([0, 1, 2], [0, 1, 4], color="red")
    ax.set_xlabel(r"$f$ (Hz)")
    ax.set_ylabel(r"$|Z|$ ($\Omega$)")
    out = tmp_path / "fig_cvd_off"
    emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.40,
        on_fail="raise", verbose=False,
        check_colorblind_safe=False,
        check_legend_overlap=False,
    )
    plt.close(fig)
    assert out.with_suffix(".pdf").exists()


# --------------------------------------------------------------------
# v0.80.0: PaperProfile.from_base (tueplots-style font derivation)
# --------------------------------------------------------------------


def test_from_base_derives_font_sizes_correctly():
    """from_base(base_pt=10, small_offset=1) gives
    font=10, legend=10, tick=9."""
    p = PaperProfile.from_base(
        name="test", full_name="test", width_mm=88.9, column="single",
        base_pt=10.0, small_offset=1.0,
    )
    assert p.font_pt == 10.0
    assert p.legend_pt == 10.0       # legend ALWAYS = body
    assert p.tick_pt == 9.0          # body - small_offset


def test_from_base_with_different_base():
    """Smaller base (slide deck case) derives proportionally."""
    p = PaperProfile.from_base(
        name="t", full_name="t", width_mm=82, column="single",
        base_pt=8.0, small_offset=1.0,
    )
    assert p.font_pt == 8.0
    assert p.legend_pt == 8.0
    assert p.tick_pt == 7.0


# --------------------------------------------------------------------
# v0.80.0: rel_width (tueplots-style fractional width)
# --------------------------------------------------------------------


def test_rel_width_scales_figure():
    """rel_width=0.5 of a double-column gives half-width figure."""
    fig_full, _ = paper_figure("ieee_double_column", rel_width=1.0)
    full_w = fig_full.get_size_inches()[0]
    plt.close(fig_full)
    fig_half, _ = paper_figure("ieee_double_column", rel_width=0.5)
    half_w = fig_half.get_size_inches()[0]
    plt.close(fig_half)
    assert abs(half_w * 2 - full_w) < 0.01, (
        f"rel_width=0.5 produced {half_w:.3f} in, expected "
        f"{full_w/2:.3f} in (half of {full_w:.3f})"
    )


def test_rel_width_1_5_for_1_5_column_figure():
    """rel_width=1.5 of single-column = 1.5-column figure."""
    fig, _ = paper_figure("ieee_single_column", rel_width=1.5)
    w = fig.get_size_inches()[0]
    plt.close(fig)
    expected_w = (88.9 / 25.4) * 1.5
    assert abs(w - expected_w) < 0.01, (
        f"rel_width=1.5 produced {w:.3f} in, expected {expected_w:.3f} in"
    )


# --------------------------------------------------------------------
# v0.80.0: auto panel labels on multi-panel by default
# --------------------------------------------------------------------


def test_panel_labels_auto_applied_on_multipanel():
    """By default, multi-panel layouts get (a)(b)(c) automatically."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=2)
    # Without explicit panel_labels=True, auto should apply for 1x2
    text_ax0 = axes[0, 0].texts[0].get_text() if axes[0, 0].texts else ""
    text_ax1 = axes[0, 1].texts[0].get_text() if axes[0, 1].texts else ""
    plt.close(fig)
    assert text_ax0 == "(a)", f"ax0 text = {text_ax0!r}, expected '(a)'"
    assert text_ax1 == "(b)", f"ax1 text = {text_ax1!r}, expected '(b)'"


def test_panel_labels_not_applied_on_single_panel():
    """1x1 figures get NO auto-label (would clutter a single panel)."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    has_label = bool(axes[0, 0].texts)
    plt.close(fig)
    assert not has_label, (
        "single-panel figure should not auto-receive (a) label"
    )


def test_panel_labels_can_be_forced_on_single_panel():
    """Explicit panel_labels=True applies even on 1x1."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1,
                              panel_labels=True)
    has_label = (axes[0, 0].texts[0].get_text() == "(a)"
                 if axes[0, 0].texts else False)
    plt.close(fig)
    assert has_label, "panel_labels=True should force (a) on 1x1"


def test_panel_labels_false_suppresses_on_multipanel():
    """Explicit panel_labels=False keeps multi-panel clean."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=2,
                              panel_labels=False)
    has_any = any(ax.texts for ax in axes.flat)
    plt.close(fig)
    assert not has_any, "panel_labels=False should suppress all labels"


# --------------------------------------------------------------------
# v0.80.0: PDF font embedding verification
# --------------------------------------------------------------------


def test_emit_pdf_has_no_type3_fonts(tmp_path):
    """After emit_paper_figure(), the saved PDF must not contain
    /Subtype /Type3 (raster-glyph) entries."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    _populate(axes)
    out = tmp_path / "fig_pdfcheck"
    res = emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.40, on_fail="raise", verbose=False,
        check_legend_overlap=False,
    )
    plt.close(fig)
    assert out.with_suffix(".pdf").exists()
    # No font violations should be reported
    assert res["font_violations"] == [], (
        f"font_violations should be empty for default Type-42 setup; "
        f"got {res['font_violations']}"
    )


def test_emit_pdf_check_can_be_disabled(tmp_path):
    """`check_font_embedding=False` skips the post-save scan."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    _populate(axes)
    out = tmp_path / "fig_pdf_nocheck"
    res = emit_paper_figure(
        fig, str(out), "ieee_double_column",
        min_axes_fraction=0.40, on_fail="raise", verbose=False,
        check_legend_overlap=False,
        check_font_embedding=False,
    )
    plt.close(fig)
    assert res["font_violations"] == []   # not run -> empty


def test_emit_raises_on_type3_fonts(tmp_path):
    """If the user overrides pdf.fonttype to 3, the gate must catch it."""
    fig, axes = paper_figure("ieee_double_column", nrows=1, ncols=1)
    _populate(axes)
    # Sabotage: force Type-3 fonts on the SAVE itself.
    plt.rcParams["pdf.fonttype"] = 3
    out = tmp_path / "fig_type3"
    try:
        with pytest.raises(ValueError, match=r"font-embedding"):
            emit_paper_figure(
                fig, str(out), "ieee_double_column",
                min_axes_fraction=0.40, on_fail="raise", verbose=False,
                check_legend_overlap=False,
            )
    finally:
        # Restore to avoid leaking sabotaged rcParam to other tests
        plt.rcParams["pdf.fonttype"] = 42
        plt.close(fig)
