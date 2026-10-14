"""Golden tests for radia_mcp.chart2d (22 paper-quality 2D chart types).

Locks in:
  - Every chart_type in CHARTS resolves through get_spec (incl. aliases)
  - Every chart_type produces a non-empty Python recipe
  - Representative chart types render to valid PNG bytes server-side
  - render_chart inherits the graph profile's Okabe-Ito palette
  - chart2d_catalog tool shape

Skipped (with reason) if matplotlib is unavailable.
"""
from __future__ import annotations

import pytest

matplotlib = pytest.importorskip("matplotlib")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from radia_mcp.chart2d import (  # noqa: E402
    CHARTS, ALIASES, get_spec, chart_groups, chart_recipe,
    render_chart,
)


# --------------------------------------------------------------------
# Catalog invariants
# --------------------------------------------------------------------


def test_catalog_has_22_types():
    """Pin the count so a silent removal is caught."""
    assert len(CHARTS) == 22, (
        f"chart catalog has {len(CHARTS)} types; expected 22"
    )


def test_catalog_groups_are_complete():
    """Every chart belongs to exactly one of 4 canonical groups."""
    expected_groups = {"line", "dist", "field", "point"}
    g = chart_groups()
    assert set(g.keys()) == expected_groups, (
        f"groups = {set(g.keys())} != {expected_groups}"
    )
    # Every chart is in exactly one group
    seen = []
    for grp, names in g.items():
        seen.extend(names)
    assert sorted(seen) == sorted(CHARTS.keys()), (
        "chart-name <-> group asymmetry"
    )


def test_aliases_resolve_via_get_spec():
    """get_spec('nyquist') resolves to the 'phase' canonical."""
    for alias, canonical in ALIASES.items():
        spec = get_spec(alias)
        assert spec.name == canonical, (
            f"get_spec({alias!r}) returned {spec.name!r}, "
            f"expected {canonical!r}"
        )


def test_get_spec_unknown_raises_with_hint():
    """Unknown chart type raises ValueError with the available list."""
    with pytest.raises(ValueError, match=r"Unknown chart type"):
        get_spec("not_a_real_chart_type")


# --------------------------------------------------------------------
# Recipe generation: 22 + every alias
# --------------------------------------------------------------------


@pytest.mark.parametrize("chart_type", sorted(CHARTS.keys()))
def test_recipe_for_every_chart_type(chart_type):
    """Every chart type's recipe must contain paper_figure(...) +
    emit_paper_figure(...) (lab style gates included)."""
    recipe = chart_recipe(chart_type)
    assert "paper_figure(" in recipe, (
        f"recipe[{chart_type}] missing paper_figure() call"
    )
    assert "emit_paper_figure(" in recipe, (
        f"recipe[{chart_type}] missing emit_paper_figure() gate"
    )
    # Recipe must be long enough to be useful
    assert len(recipe) >= 600, (
        f"recipe[{chart_type}] only {len(recipe)} chars (suspiciously short)"
    )


def test_recipe_respects_profile_arg():
    """Passing profile='ieej_double_column' must appear in the recipe."""
    recipe = chart_recipe("line", profile="ieej_double_column")
    assert "'ieej_double_column'" in recipe, (
        "profile= kwarg not threaded into recipe"
    )


# --------------------------------------------------------------------
# Server-side render: representative chart types
# --------------------------------------------------------------------
#
# Picking one per group + the Bode pair (multi-axes) + polar (special).
# Running all 22 server-side would 4x the test time without proportional
# coverage; the representatives exercise every code path in _render.


def _png_header_ok(b: bytes) -> bool:
    return b[:8] == b"\x89PNG\r\n\x1a\n"


def test_render_line_produces_valid_fig():
    fig, axes = render_chart(
        "line",
        x=np.linspace(0, 10, 50),
        y=np.linspace(0, 10, 50) ** 1.5,
        xlabel=r"$f$ (Hz)", ylabel=r"$|Z|$ ($\Omega$)",
    )
    assert axes.shape == (1, 1)
    # Save and verify PNG header
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100)
    plt.close(fig)
    assert _png_header_ok(buf.getvalue())


def test_render_bode_produces_2x1_layout():
    """Bode pair must produce a 2x1 panel layout."""
    f = np.logspace(0, 6, 100)
    omega = 2 * np.pi * f
    H = 1.0 / (1.0 + 1j * omega / (2 * np.pi * 1e3))
    fig, axes = render_chart(
        "bode",
        freq=f, magnitude=20 * np.log10(np.abs(H)),
        phase=np.angle(H, deg=True),
    )
    assert axes.shape == (2, 1), (
        f"bode layout {axes.shape}, expected (2, 1)"
    )
    plt.close(fig)


def test_render_contourf_produces_valid_fig():
    X, Y = np.meshgrid(np.linspace(-2, 2, 30), np.linspace(-2, 2, 30))
    Z = np.exp(-(X ** 2 + Y ** 2))
    fig, axes = render_chart("contourf", X=X, Y=Y, Z=Z,
                              xlabel="x (mm)", ylabel="y (mm)")
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100)
    plt.close(fig)
    assert _png_header_ok(buf.getvalue())


def test_render_histogram_produces_valid_fig():
    fig, axes = render_chart("histogram",
                              data=np.random.normal(size=300),
                              bins=20)
    plt.close(fig)


def test_render_phase_sets_equal_aspect():
    """Nyquist / impedance locus must use equal aspect (circles look round)."""
    z = (1 + 1j * np.linspace(0, 5, 50)) / (1 - 0.3j * np.linspace(0, 5, 50))
    fig, axes = render_chart("phase", z=z)
    ax = axes[0, 0]
    # matplotlib stores 'equal' as 1.0 (float)
    assert ax.get_aspect() == 1.0, (
        f"phase plot aspect = {ax.get_aspect()}, expected 1.0 ('equal')"
    )
    plt.close(fig)


def test_render_polar_uses_polar_axes():
    fig, axes = render_chart("polar",
                              theta=np.linspace(0, 2 * np.pi, 100),
                              r=1 + 0.5 * np.cos(3 * np.linspace(0, 2 * np.pi, 100)))
    ax = axes[0, 0]
    assert ax.name == "polar", (
        f"polar axes name = {ax.name!r}, expected 'polar'"
    )
    plt.close(fig)


def test_render_inherits_okabe_ito_palette():
    """Server-side render must use the lab CVD-safe palette by default."""
    from radia_mcp.figure import OKABE_ITO
    fig, axes = render_chart("line",
                              x=np.linspace(0, 1, 10),
                              y=np.linspace(0, 1, 10))
    # First line's color should be OKABE_ITO[0] = black '#000000'
    color = axes[0, 0].lines[0].get_color()
    plt.close(fig)
    import matplotlib.colors as mcolors
    rgb = mcolors.to_rgb(color)
    assert rgb == (0.0, 0.0, 0.0), (
        f"first line color = {color}, expected black "
        f"(Okabe-Ito index 0)"
    )


def test_render_unknown_raises():
    """Unknown chart type rejection."""
    with pytest.raises(ValueError, match=r"Unknown chart type"):
        render_chart("definitely_not_a_chart")


# --------------------------------------------------------------------
# chart2d_catalog tool
# --------------------------------------------------------------------


def test_catalog_tool_shape():
    """chart2d_catalog() returns the documented dict shape."""
    from radia_mcp.chart2d.server import chart2d_catalog
    cat = chart2d_catalog("all")
    assert cat["n_types"] == 22
    assert set(cat["groups"]) == {"line", "dist", "field", "point"}
    assert isinstance(cat["types"], list)
    assert len(cat["types"]) == 22
    for t in cat["types"]:
        for k in ("name", "group", "description"):
            assert k in t, f"catalog type missing {k}"
    assert "ieee_single_column" in cat["available_profiles"]
    assert set(cat["return_modes"]) == {"recipe", "image", "both"}


def test_catalog_tool_filters_by_group():
    from radia_mcp.chart2d.server import chart2d_catalog
    cat = chart2d_catalog("field")
    field_charts = {"contour", "contourf", "pcolormesh", "quiver",
                    "streamplot", "imshow", "polar"}
    assert set(t["name"] for t in cat["types"]) == field_charts
