"""Tests for radia.em_material.EMMaterial surface-impedance helpers.

dowell_Zs must return the correct THICK-conductor asymptotic limit for a large
xi = R / delta rather than silently returning nan: np.tanh overflows past
xi ~ 710 emitting a RuntimeWarning (NOT an OverflowError), so the old
`except OverflowError` never fired and dowell_Zs silently returned nan for thick
conductors (e.g. an induction-heating steel workpiece at high frequency).
"""
import math

import numpy as np

from radia.em_material import EMMaterial


def test_dowell_Zs_thick_conductor_no_nan():
    """xi = R/delta > 710 (the overflow regime where the old code silently
    returned nan): dowell_Zs returns the finite asymptotic surface impedance
    Z_s -> (1+1j) rho / delta."""
    m = EMMaterial("steel", sigma=1.0e6, mu_r=100.0)
    R = 0.01
    for f in (1.0e8, 1.0e9, 1.0e10):          # xi ~ 2000, 6300, 20000 (>> 710)
        delta = m.skin_depth(f)
        xi = R / delta
        assert xi > 710.0                      # the np.tanh overflow regime
        Zs = m.dowell_Zs(f, R)
        assert np.isfinite(Zs.real) and np.isfinite(Zs.imag), (xi, Zs)
        Z_limit = complex(1, 1) * m.rho / delta
        assert abs(Zs - Z_limit) / abs(Z_limit) < 1.0e-6, (xi, Zs, Z_limit)


def test_dowell_Zs_continuous_at_threshold():
    """The tanh branch and the xi > 30 asymptotic branch agree at the xi = 30
    threshold (no discontinuity in Z_s)."""
    m = EMMaterial("cu", sigma=5.8e7, mu_r=1.0)
    R = 0.01
    MU0 = 4.0e-7 * math.pi

    def f_for_xi(xi):
        return (xi / R) ** 2 * 2.0 * m.rho / (MU0 * m.mu_r) / (2.0 * math.pi)

    Z_lo = m.dowell_Zs(f_for_xi(29.99), R)     # tanh branch
    Z_hi = m.dowell_Zs(f_for_xi(30.01), R)     # asymptotic branch
    assert abs(Z_hi - Z_lo) / abs(Z_lo) < 1.0e-3


def test_dowell_Zs_thin_conductor_matches_tanh():
    """For a thin conductor (small xi) dowell_Zs uses the full tanh formula,
    not the asymptotic limit."""
    m = EMMaterial("cu", sigma=5.8e7, mu_r=1.0)
    R = 1.0e-4                                  # thin -> small xi
    f = 1.0e4
    delta = m.skin_depth(f)
    xi = R / delta
    assert xi < 30.0                           # the tanh-branch regime
    gamma = complex(1, 1) * xi
    expect = (m.rho / R) * gamma * np.tanh(gamma)
    Zs = m.dowell_Zs(f, R)
    assert abs(Zs - expect) / abs(expect) < 1.0e-12


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
