"""Surrogate contracts for balanced (bisection) equal-flux seeding and the charts.

Assumes closed-form cell fluxes:
- g = 1 + xi * eta on the unit square: a rectangle [x0, x1] x [y0, y1]
  carries (x1 - x0)(y1 - y0) + (x1^2 - x0^2)(y1^2 - y0^2) / 4;
- g = xi - 1/2 carries 1/8 on each side of xi = 1/2, and a cell [x0, x1]
  on one side carries |(x1 - 1/2)^2 - (x0 - 1/2)^2| (y1 - y0) / 2;
- B = z_hat through a planar annulus r in [a, b]: a cell carries
  (r1^2 - r0^2) dtheta / 2, total pi (b^2 - a^2);
- B = (x, y, 0) (1 + z / H) through the cylinder r = R, z in [0, H]: a cell
  carries R^2 dtheta H (dz' + (z1'^2 - z0'^2) / 2) with z' = z / H.
The table reproduces each density exactly because it is bilinear in the chart.
"""

import numpy as np
import pytest

from radia.flux_seeding import (annular_patch, balanced_cell_flux_quadrature,
                                balanced_equal_flux_seeds, cylindrical_patch,
                                equal_flux_seeds, planar_patch)

UNIT = planar_patch([0, 0, 0], [1, 0, 0], [0, 1, 0])


def _seed(density, n, nodes=5, chart=UNIT, **kwargs):
    grid = np.linspace(0, 1, nodes)
    return balanced_equal_flux_seeds(density, xi_nodes=grid, eta_nodes=grid, n_lines=n,
                                     chart=chart.chart, metric=chart.metric, **kwargs)


@pytest.mark.parametrize("n", [1, 7, 12])
def test_every_cell_carries_the_same_closed_form_flux(n):
    """Guards: bisection fractions, exact marginal inversion and full coverage."""
    seeds = _seed(lambda xi, eta: 1 + xi * eta, n, nodes=3)
    x0, x1, y0, y1 = seeds.cells.T
    flux = (x1 - x0) * (y1 - y0) + (x1**2 - x0**2) * (y1**2 - y0**2) / 4
    np.testing.assert_allclose(flux, np.full(n, 1.25 / n), rtol=1e-13)
    assert np.sum((x1 - x0) * (y1 - y0)) == pytest.approx(1.0, abs=1e-14)
    assert (seeds.total_flux_Wb, seeds.target_cell_flux_Wb) == pytest.approx((1.25, 1.25 / n))
    # Each seed lies inside its own cell.
    xi, eta = seeds.chart_coords.T
    assert np.all((x0 < xi) & (xi < x1) & (y0 < eta) & (eta < y1))
    np.testing.assert_array_equal(seeds.ids, np.arange(n))


@pytest.mark.parametrize("sign, expected", [("positive", 1), ("negative", -1)])
def test_one_sign_of_a_sign_changing_surface_is_seeded(sign, expected):
    """Guards: the clipped density, both flux parts and the seeded side."""
    density = lambda xi, eta: xi - 0.5 + 0 * eta
    seeds = _seed(density, 5, sign=sign)
    assert seeds.orientation == expected
    assert seeds.total_flux_Wb == pytest.approx(expected * 0.125, abs=1e-15)
    assert (seeds.positive_flux_Wb, seeds.negative_flux_Wb) == pytest.approx((0.125, -0.125))
    xi = seeds.chart_coords[:, 0]
    assert np.all(xi > 0.5) if expected > 0 else np.all(xi < 0.5), xi
    # Clipped g is linear in each cell on the side it covers: (x - 1/2) dx integrates exactly.
    x0, x1, y0, y1 = seeds.cells.T
    lo, hi = (0.5, 1.0) if expected > 0 else (0.0, 0.5)
    a, b = np.clip(x0, lo, hi), np.clip(x1, lo, hi)
    flux = np.abs(((b - 0.5) ** 2 - (a - 0.5) ** 2) / 2) * (y1 - y0)
    np.testing.assert_allclose(flux, np.full(5, 0.125 / 5), rtol=1e-12)
    check = balanced_cell_flux_quadrature(density, seeds, xi_nodes=np.linspace(0, 1, 5),
                                          eta_nodes=np.linspace(0, 1, 5))
    assert check["max_relative_cell_deviation"] < 1e-12 and check["coverage"] == pytest.approx(1)


def test_sign_change_requires_an_explicit_side():
    """Guards: auto never silently picks one side of a mixed surface."""
    with pytest.raises(ValueError, match="changes sign.*sign='positive'"):
        _seed(lambda xi, eta: xi - 0.5 + 0 * eta, 4)


@pytest.mark.parametrize("density, kwargs, message", [
    (lambda xi, eta: 1 + 0 * xi * eta, {"sign": "negative"}, "no negative flux"),
    (lambda xi, eta: 0 * xi * eta, {}, "identically zero"),
    (lambda xi, eta: 1 + 0 * xi * eta, {"sign": "both"}, "sign must be one of"),
    (lambda xi, eta: np.where(xi > 0.6, np.inf, 1.0) + 0 * eta, {}, "finite"),
])
def test_invalid_balanced_inputs_fail_loudly(density, kwargs, message):
    """Guards: absent, zero, nonfinite and unknown-sign requests raise."""
    with pytest.raises(ValueError, match=message):
        _seed(density, 4, **kwargs)
    with pytest.raises(ValueError, match="n_lines must be positive"):
        _seed(lambda xi, eta: 1 + 0 * xi * eta, 0)


def test_cells_stay_near_square_where_the_grid_gives_strips():
    """Guards: cuts chosen by the children's attainable aspect, not a fixed axis.

    A concentrated flux spot stretches the 6 x 6 grid's cells; the bisection
    of the same 36 seeds keeps them markedly squarer, and both layouts still
    carry equal flux to the table's accuracy.
    """
    spot = lambda xi, eta: 1 + 50 * np.exp(-((xi - 0.3) ** 2 + (eta - 0.6) ** 2) / 0.01)
    nodes = np.linspace(0, 1, 65)
    balanced = balanced_equal_flux_seeds(spot, xi_nodes=nodes, eta_nodes=nodes, n_lines=36,
                                         chart=UNIT.chart, metric=UNIT.metric)
    grid = equal_flux_seeds(spot, xi_nodes=nodes, eta_nodes=nodes, n_u=6, n_v=6,
                            chart=UNIT.chart, metric=UNIT.metric)
    b, g = balanced.diagnostics, grid.diagnostics
    assert b["max_cell_aspect"] < 6 < 2 * 6 < g["max_cell_aspect"], (b, g)
    assert b["median_cell_aspect"] < g["median_cell_aspect"], (b, g)
    check = balanced_cell_flux_quadrature(spot, balanced, xi_nodes=nodes, eta_nodes=nodes)
    assert check["max_relative_cell_deviation"] < 0.02 and check["coverage"] == pytest.approx(1)


def test_annulus_cells_carry_equal_flux_and_locate_inverts_the_chart():
    """Guards: polar chart, its area element, metric and locate()."""
    a, b = 0.01, 0.03
    patch = annular_patch([0.1, 0.2, 0.3], [0, 0, 2], b, r_inner_m=a)
    field = lambda p: np.broadcast_to([0.0, 0.0, 1.0], np.shape(p))
    seeds = _seed(patch.flux_density(field), 12, nodes=9, chart=patch)
    assert seeds.total_flux_Wb == pytest.approx(np.pi * (b**2 - a**2), rel=1e-13)
    x0, x1, y0, y1 = seeds.cells.T
    r0, r1 = a + x0 * (b - a), a + x1 * (b - a)
    np.testing.assert_allclose((r1**2 - r0**2) * (y1 - y0) * 2 * np.pi / 2,
                               np.full(12, seeds.target_cell_flux_Wb), rtol=1e-12)
    s, xi, eta = patch.locate(seeds.points_m)
    np.testing.assert_allclose(s, 0, atol=1e-15)
    np.testing.assert_allclose(np.column_stack((xi, eta)), seeds.chart_coords, atol=1e-12)
    np.testing.assert_allclose(seeds.points_m[:, 2], 0.3, atol=1e-15)


def test_cylinder_flux_points_outward_and_cells_are_equal():
    """Guards: cylinder chart orientation (outward), height scaling and locate()."""
    radius, height = 0.05, 0.2
    patch = cylindrical_patch([0, 0, 0], [0, 0, height], radius)
    field = lambda p: np.asarray(p) * [1.0, 1.0, 0.0] * (1 + np.asarray(p)[..., 2:3] / height)
    seeds = _seed(patch.flux_density(field), 10, nodes=9, chart=patch)
    total = radius**2 * 2 * np.pi * height * 1.5
    assert (seeds.orientation, seeds.total_flux_Wb) == (1, pytest.approx(total, rel=1e-13))
    x0, x1, y0, y1 = seeds.cells.T
    flux = radius**2 * 2 * np.pi * (x1 - x0) * height * ((y1 - y0) + (y1**2 - y0**2) / 2)
    np.testing.assert_allclose(flux, np.full(10, total / 10), rtol=1e-12)
    s, xi, eta = patch.locate(seeds.points_m)
    np.testing.assert_allclose(s, 0, atol=1e-15)
    np.testing.assert_allclose(np.column_stack((xi, eta)), seeds.chart_coords, atol=1e-12)


@pytest.mark.parametrize("build, message", [
    (lambda: annular_patch([0, 0, 0], [0, 0, 1], 0.01, r_inner_m=0.02), "radii"),
    (lambda: annular_patch([0, 0, 0], [0, 0, 0], 0.01), "nonzero"),
    (lambda: annular_patch([0, 0, 0], [0, 0, 1], 0.01, theta_range=(0, 7)), "theta_range"),
    (lambda: cylindrical_patch([0, 0, 0], [0, 0, 1], -1.0), "radius_m"),
    (lambda: cylindrical_patch([0, 0, 0], [0, 0, 1], 1.0, reference=[0, 0, 2]), "parallel"),
])
def test_degenerate_charts_are_rejected(build, message):
    """Guards: inverted radii, zero axes, overlong sectors and parallel references."""
    with pytest.raises(ValueError, match=message):
        build()
