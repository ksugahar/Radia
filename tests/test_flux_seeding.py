"""Surrogate contracts for equal-flux seeding on prescribed charts."""

import numpy as np
import pytest

from radia.flux_seeding import equal_flux_seeds


def sector_chart(radius, v):
    phi = np.pi * v / 2
    return np.stack((radius * np.cos(phi), radius * np.sin(phi), np.zeros_like(phi)), axis=-1)


@pytest.mark.parametrize("sign", [1, -1])
def test_annular_sector_seeds_match_analytic_inverse(sign):
    """Guards: signed cell weights, chart inversion and covered patch flux."""
    u_edges = np.linspace(0, 0.5, 9)
    if sign < 0:
        u_edges = -u_edges[::-1]
    result = equal_flux_seeds(
        sector_chart, lambda r, v: sign * np.pi * r / 2,
        s_nodes=np.linspace(1, 2, 1001), v_nodes=[0, 0.5, 1],
        u_edges=u_edges, v_edges=np.linspace(0, 1, 5), s_reference=1,
    )
    centres_u = result.cells_uv[:, :2].mean(axis=1)
    centres_v = result.cells_uv[:, 2:].mean(axis=1)
    expected_radius = np.sqrt(1 + 4 * sign * centres_u / np.pi)
    np.testing.assert_allclose(result.points_m, sector_chart(expected_radius, centres_v), atol=2e-7)
    assert len(result.points_m) == 32
    assert result.signed_cell_flux_Wb == pytest.approx(sign / 64)
    assert result.signed_patch_flux_Wb == pytest.approx(sign * 3 * np.pi / 4)
    assert result.coverage == pytest.approx(2 / (3 * np.pi))


def test_complete_cells_require_interior_v_knots():
    """Guards: a row can contain a boundary minimum between its two edges."""
    result = equal_flux_seeds(
        lambda s, v: np.column_stack((s, v, np.zeros_like(s))),
        lambda s, v: 0.5 + np.abs(v - 0.5),
        s_nodes=[0, 0.5, 1], v_nodes=[0, 0.5, 1],
        u_edges=[0, 0.25, 0.5, 0.75, 1], v_edges=[0, 1], s_reference=0,
    )
    np.testing.assert_array_equal(result.cells_uv, [[0, 0.25, 0, 1], [0.25, 0.5, 0, 1]])
    np.testing.assert_allclose(result.points_m[:, 0], [0.25, 0.75])
    assert result.coverage == pytest.approx(2 / 3)


def test_sign_change_is_rejected():
    """Guards: inversion fails loudly on a patch with changing flux sign."""
    with pytest.raises(ValueError, match="one nonzero sign"):
        equal_flux_seeds(
            sector_chart, lambda r, v: v - 0.5,
            s_nodes=[1, 1.5, 2], v_nodes=[0, 0.5, 1],
            u_edges=[0, 0.25], v_edges=[0, 1], s_reference=1,
        )
