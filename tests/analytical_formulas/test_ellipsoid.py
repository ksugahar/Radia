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


# ---------------------------------------------------------------------------
# Osborn 1945 reference values (Phys. Rev. 67, 351, Tables I and II)
# ---------------------------------------------------------------------------

# Polar demag factor along the symmetry axis, computed from the closed forms
# in Osborn, J.A., "Demagnetizing Factors of the General Ellipsoid",
# Phys. Rev. 67 (1945) 351-357 with high-precision SymPy. Aspect ratio
# m = c / a (polar / equatorial).
OSBORN_REFERENCE = {
    # (m, expected N_z) -- 18-digit values produced by an independent
    # mpmath script (mp.dps = 30) using the same closed forms.
    "prolate, m=2":   (2.0,  0.17356399753396423),
    "prolate, m=5":   (5.0,  0.0558209698024552),
    "prolate, m=10":  (10.0, 0.02028588030156383),
    "prolate, m=100": (100.0, 0.000429898719881890),
    "oblate, m=0.5":  (0.5,  0.5272002825625699),
    "oblate, m=0.1":  (0.1,  0.8608042765278006),
    "oblate, m=0.01": (0.01, 0.9844897069128692),
}


@pytest.mark.parametrize("label", list(OSBORN_REFERENCE))
def test_against_osborn_reference(label):
    m, expected = OSBORN_REFERENCE[label]
    Nz = demag_factor_rotational(m, 1.0)[2]
    # double precision input is the bottleneck; allow ~1e-14 absolute drift.
    assert Nz == pytest.approx(expected, rel=1e-13, abs=1e-14)


# ---------------------------------------------------------------------------
# Limit consistency: the half-difference between branches at fixed |c/a - 1|
# ---------------------------------------------------------------------------


def test_branches_agree_with_taylor_series():
    # For very small (c/a - 1) the two branches must match the Taylor series
    # N_z = 1/3 - 2 eps / 15 + 8 eps**2 / 105 - O(eps**3)  (eps = (c/a)**2 - 1).
    # The next-order coefficient is ~ -8/315, so the truncation residual is
    # bounded by ~ |eps|**3 / 30.
    eps = 5e-3
    Nz_pro = demag_factor_prolate(1.0 + eps, 1.0)
    Nz_obl = demag_factor_oblate(1.0 - eps, 1.0)
    eps_pro = (1.0 + eps) ** 2 - 1.0
    eps_obl = (1.0 - eps) ** 2 - 1.0
    Nz_pro_series = 1.0 / 3.0 - 2 * eps_pro / 15 + 8 * eps_pro ** 2 / 105
    Nz_obl_series = 1.0 / 3.0 - 2 * eps_obl / 15 + 8 * eps_obl ** 2 / 105
    assert abs(Nz_pro - Nz_pro_series) < abs(eps_pro) ** 3
    assert abs(Nz_obl - Nz_obl_series) < abs(eps_obl) ** 3
