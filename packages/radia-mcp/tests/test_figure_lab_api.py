"""Tests for the misuse-proof radia_mcp.figure lab API + builders + audit.

These lock the fail-loud behaviour that was missing when the CEFC 2026 oral
figures slipped through with in-figure titles + \\linewidth downscaling.
"""
import os
import pytest

# Skip the whole module if matplotlib is not installed (CI minimal env).
# These tests need matplotlib at import time (Agg backend must be set before
# any figure submodule uses it).  The radia-mcp-matrix CI installs only mcp;
# this guard prevents collection-time ImportError.
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from radia_mcp.figure import (
    lab_figure, save_lab_figure, legend_no_overlap, audit_tex_figures,
    audit_label_overflow, scaling_loglog, bh_curve, check_min_font,
    lab_savefig,
)


def _tnr_available():
    import matplotlib.font_manager as fm
    try:
        p = fm.findfont("Times New Roman", fallback_to_default=False)
        return "times" in os.path.basename(str(p)).lower()
    except Exception:
        return False


requires_tnr = pytest.mark.skipif(
    not _tnr_available(), reason="Times New Roman not installed in matplotlib")


# ------------------------------------------------------------------
# fail-loud gates
# ------------------------------------------------------------------
def test_font_gate_rejects_paper_annotation_below_10pt():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.0 / 2.54, 2.0))
    ax.text(0.5, 0.5, "small annotation", fontsize=9.9)
    bad = check_min_font(fig, min_pt=10.0, embed_width_cm=8.0)
    assert len(bad) == 1
    assert bad[0]["text"] == "small annotation"
    assert bad[0]["visible_pt"] == pytest.approx(9.9)
    plt.close(fig)


def test_font_gate_rejects_oversized_paper_annotation():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.0 / 2.54, 2.0))
    ax.text(0.5, 0.5, "oversized annotation", fontsize=11.0)
    bad = check_min_font(
        fig, min_pt=10.0, max_pt=10.5, embed_width_cm=8.0
    )
    assert len(bad) == 1
    assert bad[0]["reason"] == "above maximum"
    assert bad[0]["visible_pt"] == pytest.approx(11.0)
    plt.close(fig)


def test_font_gate_uses_actual_powerpoint_paste_width():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(25.0 / 2.54, 4.0))
    ax.set_axis_off()
    ax.text(0.5, 0.5, "source 24 pt", fontsize=24.0)
    bad = check_min_font(fig, min_pt=20.0, embed_width_cm=20.0)
    assert len(bad) == 1
    assert bad[0]["embed_scale"] == pytest.approx(0.8)
    assert bad[0]["visible_pt"] == pytest.approx(19.2)
    plt.close(fig)


def test_font_gate_accepts_exact_20pt_after_powerpoint_scaling():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(25.0 / 2.54, 4.0))
    ax.set_axis_off()
    ax.text(0.5, 0.5, "source 25 pt", fontsize=25.0)
    assert check_min_font(fig, min_pt=20.0, embed_width_cm=20.0) == []
    plt.close(fig)


def test_lab_savefig_rejects_downscaled_slide_text(tmp_path):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(25.0 / 2.54, 4.0))
    ax.set_axis_off()
    ax.text(0.5, 0.5, "source 24 pt", fontsize=24.0)
    with pytest.raises(ValueError, match="Required for presentation: >= 20 pt"):
        lab_savefig(
            fig,
            str(tmp_path / "bad_slide.png"),
            medium="presentation",
            embed_width_cm=20.0,
        )
    plt.close(fig)


def test_save_lab_figure_rejects_downscaled_slide_text(tmp_path):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(25.0 / 2.54, 4.0))
    ax.set_axis_off()
    ax.text(0.5, 0.5, "source 24 pt", fontsize=24.0)
    with pytest.raises(ValueError, match="Required for presentation: >= 20 pt"):
        save_lab_figure(
            fig,
            str(tmp_path / "bad_slide"),
            embed_width_cm=20.0,
            medium="presentation",
            save_pdf=False,
            save_png=False,
            check_times_new_roman=False,
            check_label_overflow=False,
        )
    plt.close(fig)


@requires_tnr
def test_lab_figure_authors_at_embed_width():
    import matplotlib.pyplot as plt
    fig, ax = lab_figure(8.0, aspect=0.6)
    # 8.0 cm = 3.150 in authored width (apply_lab_style returns inches)
    assert abs(fig.get_size_inches()[0] - 8.0/2.54) < 0.05
    assert getattr(fig, "_lab_embed_width_cm", None) == 8.0
    plt.close(fig)


@requires_tnr
def test_save_raises_on_in_figure_title(tmp_path):
    fig, ax = lab_figure(8.0)
    ax.plot([0, 1], [0, 1])
    ax.set_title("this belongs in the caption")
    with pytest.raises(ValueError):
        save_lab_figure(fig, str(tmp_path / "bad"), 8.0, save_pdf=False)


@requires_tnr
def test_save_ok_emits_exact_latex_width(tmp_path):
    fig, ax = lab_figure(8.0)
    ax.plot([0, 1], [0, 1]); ax.set_xlabel("x"); ax.set_ylabel("y")
    info = save_lab_figure(fig, str(tmp_path / "good"), 8.0, save_pdf=True, save_png=True)
    assert info["gates"] == "passed"
    assert r"width=8.00cm" in info["latex"]
    assert os.path.isfile(str(tmp_path / "good.png"))
    assert os.path.isfile(str(tmp_path / "good.pdf"))


@requires_tnr
def test_builders_round_trip(tmp_path):
    fig, ax = scaling_loglog([1e3, 1e4, 1e5],
                             {"dense": [1.0, 10.0, None], "H": [0.5, 2.0, 6.0]},
                             "memory (GB)", ram_line_gb=128)
    info = save_lab_figure(fig, str(tmp_path / "scal"), 7.6, save_pdf=False)
    assert "width=7.60cm" in info["latex"]
    fig2, ax2 = bh_curve([0, 100, 300], [0.0, 1.8, 2.5])
    info2 = save_lab_figure(fig2, str(tmp_path / "bh"), 6.0, save_pdf=False)
    assert "width=6.00cm" in info2["latex"]


# ------------------------------------------------------------------
# .tex embed audit (no matplotlib font dependency)
# ------------------------------------------------------------------
def test_audit_flags_height_linewidth_and_missing(tmp_path):
    tex = tmp_path / "deck.tex"
    tex.write_text(
        "\\graphicspath{{figures/}}\n"
        "\\includegraphics[height=3cm]{a}\n"
        "\\includegraphics[width=\\linewidth]{b}\n"
        "\\includegraphics[width=8cm]{c}\n",
        encoding="utf-8")
    rep = audit_tex_figures(str(tex))
    assert rep["n_figures"] == 3
    by = {r["figure"]: r for r in rep["figures"]}
    assert any("HEIGHT" in x for x in by["a"]["risks"])
    assert any("linewidth" in x.lower() for x in by["b"]["risks"])
    assert by["c"]["fixed_cm_width"] is True
    # all three referenced files are absent -> each flagged (FILE NOT FOUND)
    assert rep["n_flagged"] == 3


# ------------------------------------------------------------------
# label-overflow lint (margins too tight -> xlabel / ylabel clipped)
# ------------------------------------------------------------------
def test_audit_label_overflow_detects_clipped_axis_labels():
    """A full-bleed axes (position [0,0,1,1]) pushes the x- and y-labels off the canvas; the lint flags
    both with a positive overhang -- no Times New Roman needed (pure matplotlib)."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.set_position([0.0, 0.0, 1.0, 1.0])
    ax.plot([0, 1], [0, 1]); ax.set_xlabel("a long x axis label"); ax.set_ylabel("a long y axis label")
    rep = audit_label_overflow(fig, include_ticklabels=False)
    wheres = {r["where"] for r in rep}
    assert "xlabel" in wheres and "ylabel" in wheres, rep
    assert all(r["overhang_pt"] > 0 for r in rep)
    plt.close(fig)


def test_audit_label_overflow_clean_when_roomy():
    """With a normal layout that reserves room for the labels, the axis-label lint is clean."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot([0, 1], [0, 1]); ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.subplots_adjust(left=0.2, bottom=0.2, right=0.95, top=0.95)
    assert audit_label_overflow(fig, include_ticklabels=False) == []
    plt.close(fig)


def test_save_lab_figure_raises_on_unfixable_clip(tmp_path):
    """With tighten=False (no back-off) a hand-forced full-bleed axes clips the labels -> the save gate
    fails LOUD rather than ship a clipped figure (pure matplotlib, no TNR request)."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.set_position([0.0, 0.0, 1.0, 1.0])
    ax.plot([0, 1], [0, 1]); ax.set_xlabel("clipped x label"); ax.set_ylabel("clipped y label")
    with pytest.raises(ValueError, match="overflow the figure canvas"):
        save_lab_figure(
            fig, str(tmp_path / "clip"), 8.0, save_pdf=False, tighten=False,
            check_times_new_roman=False,
        )
    plt.close(fig)


@requires_tnr
def test_save_lab_figure_raises_when_times_new_roman_not_requested(tmp_path):
    """save_lab_figure() is a lab-standard save path, so bare DejaVu rcParams fail."""
    import matplotlib.pyplot as plt
    old_family = plt.rcParams["font.family"]
    old_serif = list(plt.rcParams["font.serif"])
    try:
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = ["DejaVu Serif"]
        fig, ax = plt.subplots(figsize=(3, 2))
        ax.plot([0, 1], [0, 1])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        with pytest.raises(RuntimeError, match="Times New Roman"):
            save_lab_figure(
                fig, str(tmp_path / "badfont"), 8.0, save_pdf=False,
                check_label_overflow=False,
            )
    finally:
        plt.rcParams["font.family"] = old_family
        plt.rcParams["font.serif"] = old_serif
        plt.close("all")


@requires_tnr
def test_save_lab_figure_backs_off_clipped_label(tmp_path):
    """The exact Rac/Rdc figure whose xlabel auto_tighten(0.80) clipped at the bottom edge: the back-off
    (re-run tight_layout, figure size unchanged) reserves the room so the gate passes -- the
    '余白を攻めすぎ' bottom clip never ships."""
    import numpy as np
    fig, ax = lab_figure(8.5, aspect=0.66)
    ax.semilogx(np.logspace(3, 6, 60), np.linspace(1.0, 2.5, 60), label="solid")
    ax.set_xlabel("frequency  $f$  (Hz)")
    ax.set_ylabel(r"$R_\mathrm{AC}/R_\mathrm{DC}$  (per strand)")
    ax.legend()
    info = save_lab_figure(fig, str(tmp_path / "rac"), 8.5, save_pdf=False, save_png=True)
    assert info["gates"] == "passed" and os.path.isfile(str(tmp_path / "rac.png"))


# ------------------------------------------------------------------
# legend placement + axes-to-the-limit
# ------------------------------------------------------------------
@requires_tnr
def test_legend_no_overlap_returns_a_placement():
    import numpy as np
    import matplotlib.pyplot as plt
    fig, ax = lab_figure(8.0)
    x = np.linspace(0, 1, 50)
    ax.plot(x, x, label="rising")
    ax.plot(x, 1 - x, label="falling")
    placed = legend_no_overlap(ax)
    assert placed.startswith("inside") or placed == "outside-right"
    assert ax.get_legend() is not None
    plt.close(fig)


@requires_tnr
def test_save_reports_axes_fraction(tmp_path):
    fig, ax = lab_figure(8.0)
    ax.plot([0, 1], [0, 1]); ax.set_xlabel("x"); ax.set_ylabel("y")
    info = save_lab_figure(fig, str(tmp_path / "eff"), 8.0, save_pdf=False)
    # auto_tighten ran (default) -> the axes fill a healthy fraction
    assert info["axes_fraction"] is not None
    assert 0.4 < info["axes_fraction"] <= 1.0
