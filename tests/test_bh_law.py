"""The shared soft-iron B(H) interpolant keeps dB/dH >= mu0 whenever the table's M rises."""
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from radia.bh_law import (FLUX_DENSITY_PCHIP, MAGNETIZATION_PCHIP, MU_0, bh_interpolation_kind,
                          monotone_bh_pchip)

# Rising-M rows on which PCHIP of B(H) itself dips to dB/dH = 0 between rows
KNEE = np.array([[0, 0], [100, 0.2], [300, 0.9], [600, 1.3], [1000, 1.45], [3000, 1.6],
                 [10000, 1.75], [50000, 1.9], [200000, 2.1]], dtype=float)
# M = B/mu0 - H falls above 3e4 A/m (the saturation-audit stress fixture shape)
FALLING = np.array([[0, 0], [100, 0.6], [1e3, 1.7], [1e4, 2.05], [3e4, 2.1], [1e5, 2.15]], dtype=float)


def _grid(table):
    return np.unique(np.concatenate([np.linspace(0.0, table[-1, 0], 40001),
                                     np.geomspace(1e-2, table[-1, 0], 40001)]))


def test_rising_magnetization_keeps_vacuum_slope_floor():
    H, B = KNEE.T
    assert bh_interpolation_kind(H, B) == MAGNETIZATION_PCHIP
    h = _grid(KNEE)
    assert np.min(PchipInterpolator(H, B).derivative()(h)) < 0.5*MU_0   # the defect this law removes
    law = monotone_bh_pchip(H, B)
    assert np.min(law.derivative()(h)) >= MU_0*(1.0 - 1e-12)
    assert np.all(np.diff(law(h)/MU_0 - h) >= -1e-9*np.max(B/MU_0))
    assert np.allclose(law(H), B, rtol=0.0, atol=1e-12)


def test_law_is_c1_on_the_table_knots_and_a_ppoly_drop_in():
    H, B = KNEE.T
    law = monotone_bh_pchip(H, B, extrapolate=False)
    assert law.c.shape == (4, len(H) - 1) and np.array_equal(law.x, H)
    assert np.allclose(law.c[-1], B[:-1], rtol=0.0, atol=1e-12)            # value at each left knot
    slope = law.derivative()
    for knot in H[1:-1]:
        left, right = slope(knot*(1 - 1e-12)), slope(knot*(1 + 1e-12))
        assert left == pytest.approx(right, rel=1e-6)
    assert np.isnan(law(H[-1]*1.5))                                         # caller owns the tail
    primitive = law.antiderivative()
    grid = np.linspace(0.0, H[-1], 200001)
    trapezoid = float(np.sum(0.5*(law(grid[1:]) + law(grid[:-1]))*np.diff(grid)))
    assert primitive(H[-1]) - primitive(0.0) == pytest.approx(trapezoid, rel=1e-8)


def test_falling_magnetization_keeps_monotone_flux_density():
    H, B = FALLING.T
    assert bh_interpolation_kind(H, B) == FLUX_DENSITY_PCHIP
    law = monotone_bh_pchip(H, B)
    reference = PchipInterpolator(H, B)
    h = _grid(FALLING)
    assert np.array_equal(law(h), reference(h))
    assert np.min(law.derivative()(h)) >= 0.0
