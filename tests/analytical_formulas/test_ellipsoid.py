"""Tests for radia.analytical_formulas.ellipsoid (Part 5, eq 38-44)."""

import math

import pytest

from radia.analytical_formulas import (
    demag_factor_oblate,
    demag_factor_prolate,
    demag_factor_rotational,
    ellipsoid_internal_field,
    ellipsoid_torque,
)
from radia.analytical_formulas.ellipsoid import MU_0


# ---------------------------------------------------------------------------
# Sphere limit
# ---------------------------------------------------------------------------


def test_sphere_demag_factors():
    Nx, Ny, Nz = demag_factor_rotational(1.0, 1.0)
    assert Nx == pytest.approx(1.0 / 3.0, abs=1e-12)
    assert Ny == pytest.approx(1.0 / 3.0, abs=1e-12)
    assert Nz == pytest.approx(1.0 / 3.0, abs=1e-12)


def test_sphere_via_prolate_branch():
    # Sphere is the boundary of the prolate range.
    assert demag_factor_prolate(1.0, 1.0) == pytest.approx(1.0 / 3.0, abs=1e-12)


def test_sphere_via_oblate_branch():
    assert demag_factor_oblate(1.0, 1.0) == pytest.approx(1.0 / 3.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Geometric limits
# ---------------------------------------------------------------------------


def test_long_needle_polar_factor_vanishes():
    # c >> a -- prolate, N_z -> 0.
    Nz = demag_factor_prolate(1000.0, 1.0)
    assert 0.0 < Nz < 1e-4


def test_thin_disk_polar_factor_unity():
    # c << a -- oblate, N_z -> 1.
    Nz = demag_factor_oblate(0.001, 1.0)
    assert 0.998 < Nz < 1.0


def test_demag_factors_sum_to_unity():
    for c, a in [(0.3, 1.0), (0.7, 1.0), (1.5, 1.0), (3.0, 1.0)]:
        Nx, Ny, Nz = demag_factor_rotational(c, a)
        assert Nx + Ny + Nz == pytest.approx(1.0, abs=1e-12)


def test_equatorial_factors_equal():
    for c, a in [(0.3, 1.0), (1.5, 1.0), (5.0, 1.0)]:
        Nx, Ny, _ = demag_factor_rotational(c, a)
        assert Nx == Ny


# ---------------------------------------------------------------------------
# Series-vs-direct continuity at sphere
# ---------------------------------------------------------------------------


def test_demag_factor_continuous_through_sphere():
    # |c/a - 1| = 1e-3 should match the direct branch within Taylor accuracy.
    Nz_eps_prolate = demag_factor_prolate(1.0 + 1e-3, 1.0)
    Nz_eps_oblate = demag_factor_oblate(1.0 - 1e-3, 1.0)
    # Both must straddle 1/3.
    assert Nz_eps_prolate < 1.0 / 3.0
    assert Nz_eps_oblate > 1.0 / 3.0
    # Relative deviation from 1/3 should be O((c/a - 1)).
    assert abs(Nz_eps_prolate - 1.0 / 3.0) < 5e-4
    assert abs(Nz_eps_oblate - 1.0 / 3.0) < 5e-4


def test_demag_factor_argument_validation():
    with pytest.raises(ValueError):
        demag_factor_prolate(0.5, 1.0)        # c < a is oblate
    with pytest.raises(ValueError):
        demag_factor_oblate(2.0, 1.0)         # c > a is prolate
    with pytest.raises(ValueError):
        demag_factor_prolate(-1.0, 1.0)
    with pytest.raises(ValueError):
        demag_factor_prolate(1.0, 0.0)


# ---------------------------------------------------------------------------
# Internal field
# ---------------------------------------------------------------------------


def test_internal_field_no_material():
    # chi_r = 0 means no magnetic response -> H_inside = H_0
    assert ellipsoid_internal_field(1234.5, chi_r=0.0, N=0.4) == pytest.approx(1234.5)


def test_internal_field_sphere_high_mu():
    # H_inside = H_0 / (1 + chi_r * 1/3)
    H_in = ellipsoid_internal_field(1000.0, chi_r=999.0, N=1.0 / 3.0)
    assert H_in == pytest.approx(1000.0 / (1.0 + 999.0 / 3.0), rel=1e-12)


# ---------------------------------------------------------------------------
# Torque
# ---------------------------------------------------------------------------


def test_torque_zero_when_aligned():
    for alpha in (0.0, math.pi):
        T = ellipsoid_torque(1000.0, alpha, chi_r=999.0, a=0.01, c=0.05)
        assert T == pytest.approx(0.0, abs=1e-15)


def test_torque_zero_when_perpendicular():
    # alpha = pi/2 also gives sin(2 alpha) = 0
    T = ellipsoid_torque(1000.0, math.pi / 2, chi_r=999.0, a=0.01, c=0.05)
    assert T == pytest.approx(0.0, abs=1e-15)


def test_prolate_torque_is_restoring():
    # Prolate (c > a) has N_z < 1/3 -> 1 - 3 N_z > 0 -> T_z < 0 at alpha in (0, pi/2)
    T = ellipsoid_torque(1000.0, math.pi / 4, chi_r=999.0, a=0.01, c=0.1)
    assert T < 0.0


def test_oblate_torque_is_repulsive():
    # Oblate (c < a) has N_z > 1/3 -> 1 - 3 N_z < 0 -> T_z > 0 at alpha in (0, pi/2)
    T = ellipsoid_torque(1000.0, math.pi / 4, chi_r=999.0, a=0.1, c=0.01)
    assert T > 0.0


def test_torque_quadratic_in_field_amplitude():
    # T_z scales as H_0**2.
    T1 = ellipsoid_torque(100.0, math.pi / 4, chi_r=10.0, a=0.01, c=0.02)
    T2 = ellipsoid_torque(200.0, math.pi / 4, chi_r=10.0, a=0.01, c=0.02)
    assert T2 / T1 == pytest.approx(4.0, rel=1e-12)


def test_torque_volume_scaling():
    # T_z scales as V = (4/3) pi a**2 c
    T1 = ellipsoid_torque(100.0, math.pi / 4, chi_r=10.0, a=0.01, c=0.02)
    T2 = ellipsoid_torque(100.0, math.pi / 4, chi_r=10.0, a=0.02, c=0.04)
    # Doubling both a and c gives V * 8
    assert T2 / T1 == pytest.approx(8.0, rel=1e-12)
