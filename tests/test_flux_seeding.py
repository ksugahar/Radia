"""Surrogate contracts for two-stage equal-flux seeding on a planar chart.

Assumes: the closed-form marginal and conditional inverses of the
nonseparable density g = 1 + xi * eta on the unit square,
    V(eta) = (eta + eta**2 / 4) / (5 / 4),
    U(xi | eta) = (xi + eta * xi**2 / 2) / (1 + eta / 2),
which a piecewise-linear table reproduces exactly because g is bilinear.
"""

import numpy as np
import pytest

from radia.flux_seeding import cell_flux_quadrature, equal_flux_seeds, planar_patch

UNIT = planar_patch([0, 0, 0], [1, 0, 0], [0, 1, 0])


def bilinear(xi, eta):
    return 1 + xi * eta


def wavy(xi, eta):
    return 1 + xi * eta + 0.5 * np.sin(np.pi * xi) * np.sin(np.pi * eta)


def analytic_inverse(u, v):
    eta = 2 * (np.sqrt(1 + 5 * v / 4) - 1)
    xi = (np.sqrt(1 + 2 * eta * u * (1 + eta / 2)) - 1) / eta
    return xi, eta


def test_nonseparable_density_matches_closed_form_inverse():
    """Guards: marginal-then-conditional order, quadratic inversion and seed IDs."""
    nodes = np.linspace(0, 1, 3)
    seeds = equal_flux_seeds(bilinear, xi_nodes=nodes, eta_nodes=nodes, n_u=4, n_v=5,
                             chart=UNIT.chart)
    u, v = seeds.flux_coords.T
    xi, eta = analytic_inverse(u, v)
    np.testing.assert_allclose(seeds.chart_coords, np.column_stack((xi, eta)), atol=1e-14)
    np.testing.assert_array_equal(seeds.ids, np.arange(20))
    np.testing.assert_allclose(u[:4], [0.125, 0.375, 0.625, 0.875])
    assert seeds.total_flux_Wb == pytest.approx(1.25)
    assert seeds.target_cell_flux_Wb == pytest.approx(1.25 / 20)
    assert seeds.diagnostics["forward_residual_u_max"] < 1e-14
    assert seeds.diagnostics["forward_residual_v_max"] < 1e-14


def test_conditional_is_solved_at_the_exact_seed_eta():
    """Guards: each seed row is re-sampled at its own eta, not at a table row.

    g = 1 + 3 xi eta**2 is linear in xi, so the conditional at any eta is
    exact, while the coarse eta table makes the marginal approximate.
    """
    g = lambda xi, eta: 1 + 3 * xi * eta**2
    seeds = equal_flux_seeds(g, xi_nodes=[0, 1], eta_nodes=[0, 0.5, 1], n_u=3, n_v=4,
                             chart=UNIT.chart)
    xi, eta = seeds.chart_coords.T
    u = seeds.flux_coords[:, 0]
    exact_u = (xi + 1.5 * eta**2 * xi**2) / (1 + 1.5 * eta**2)
    np.testing.assert_allclose(exact_u, u, atol=1e-14)
    assert len(np.unique(eta)) == 4


def test_negative_uniform_field_keeps_its_orientation():
    """Guards: a globally negative flux is placed on |g| and reported with its sign."""
    patch = planar_patch([0.1, -0.2, 0.3], [0.02, 0, 0], [0, 0.04, 0])
    field = lambda p: np.broadcast_to([0.0, 0.0, -0.5], np.shape(p))
    seeds = equal_flux_seeds(patch.flux_density(field), xi_nodes=[0, 1], eta_nodes=[0, 1],
                             n_u=2, n_v=2, chart=patch.chart)
    assert seeds.orientation == -1
    assert seeds.total_flux_Wb == pytest.approx(-0.5 * 0.02 * 0.04)
    assert seeds.target_cell_flux_Wb == pytest.approx(-0.5 * 0.02 * 0.04 / 4)
    np.testing.assert_allclose(seeds.points_m[:, :2],
                               [[0.105, -0.19], [0.115, -0.19], [0.105, -0.17], [0.115, -0.17]])


def test_endpoint_singularity_is_accepted():
    """Guards: a density vanishing on a chart edge is a coordinate singularity."""
    seeds = equal_flux_seeds(lambda xi, eta: xi * eta, xi_nodes=np.linspace(0, 1, 5),
                             eta_nodes=np.linspace(0, 1, 5), n_u=2, n_v=2, chart=UNIT.chart)
    # g = xi eta is bilinear: V = eta**2 and U = xi**2.
    np.testing.assert_allclose(seeds.chart_coords, np.sqrt(seeds.flux_coords), atol=1e-14)


@pytest.mark.parametrize("density, message", [
    (lambda xi, eta: xi - 0.5 + 0 * eta, "changes sign"),
    (lambda xi, eta: 0 * xi * eta, "identically zero"),
    (lambda xi, eta: np.where(xi > 0.6, np.nan, 1.0) + 0 * eta, "finite"),
    (lambda xi, eta: np.abs(xi - 0.5) + 0 * eta, "interior xi sample"),
])
def test_invalid_densities_fail_loudly(density, message):
    """Guards: sign-changing, zero, nonfinite and interior-degenerate inputs raise."""
    nodes = np.linspace(0, 1, 5)
    with pytest.raises(ValueError, match=message):
        equal_flux_seeds(density, xi_nodes=nodes, eta_nodes=nodes, n_u=2, n_v=2,
                         chart=UNIT.chart)


def test_noninvertible_planar_chart_is_rejected():
    """Guards: parallel chart vectors have no area element."""
    with pytest.raises(ValueError, match="not invertible"):
        planar_patch([0, 0, 0], [1, 0, 0], [2, 0, 0])


def test_curved_cell_quadrature_separates_table_and_quadrature_error():
    """Guards: independent cell integration over curved-in-chart xi edges.

    The table error falls at second order with refinement; raising the
    Gauss order on a fixed table changes the deviation only at quadrature
    level, and every cell set covers the whole patch flux.
    """
    deviations = []
    for n in (9, 17, 33):
        nodes = np.linspace(0, 1, n)
        seeds = equal_flux_seeds(wavy, xi_nodes=nodes, eta_nodes=nodes, n_u=6, n_v=6,
                                 chart=UNIT.chart)
        check = cell_flux_quadrature(wavy, seeds, xi_nodes=nodes, eta_nodes=nodes, order=8)
        deviations.append(check["max_relative_cell_deviation"])
        assert check["coverage"] == pytest.approx(1, abs=1e-12)
        assert check["independent_total_flux_Wb"] == pytest.approx(1.25 + 2 / np.pi**2)
    assert deviations[0] / deviations[1] > 3.5 and deviations[1] / deviations[2] > 3.5, deviations

    nodes = np.linspace(0, 1, 17)
    seeds = equal_flux_seeds(wavy, xi_nodes=nodes, eta_nodes=nodes, n_u=6, n_v=6,
                             chart=UNIT.chart)
    by_order = [cell_flux_quadrature(wavy, seeds, xi_nodes=nodes, eta_nodes=nodes,
                                     order=k)["max_relative_cell_deviation"] for k in (6, 10)]
    assert by_order[0] == pytest.approx(by_order[1], rel=1e-6), by_order
