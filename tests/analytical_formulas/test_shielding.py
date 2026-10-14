"""Tests for radia.analytical_formulas.shielding (Part 1, eq 23-24)."""

import pytest

from radia.analytical_formulas import (
    shielding_factor_cylinder,
    shielding_factor_sphere,
)


# ---------------------------------------------------------------------------
# No-shield limit
# ---------------------------------------------------------------------------


def test_cylinder_no_shielding_when_mu_unity():
    assert shielding_factor_cylinder(0.5, 1.0, mu_r=1.0) == pytest.approx(1.0, abs=1e-12)


def test_sphere_no_shielding_when_mu_unity():
    assert shielding_factor_sphere(0.5, 1.0, mu_r=1.0) == pytest.approx(1.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Asymptotes
# ---------------------------------------------------------------------------


def test_cylinder_high_mu_thin_shell_asymptote():
    # S_cyl ~ 4 b**2 / [(mu - 1)**2 (b**2 - a**2)] ~ 4 b / [mu * (b - a)] for thin shell
    a, b, mu_r = 0.99, 1.0, 100000.0
    S = shielding_factor_cylinder(a, b, mu_r)
    # The leading order is 4 b**2 / (mu_r**2 (b**2 - a**2)) when (a/b)**2 close to 1.
    # numerically just check S << 1.
    assert 0.0 < S < 1e-2


def test_sphere_high_mu_thin_shell_asymptote():
    a, b, mu_r = 0.99, 1.0, 100000.0
    S = shielding_factor_sphere(a, b, mu_r)
    assert 0.0 < S < 1e-2


def test_cylinder_solid_limit_no_inside():
    # a -> 0 (solid magnetic cylinder, no inside) -> S = 4 mu_r / (mu_r + 1)**2
    mu_r = 100.0
    S = shielding_factor_cylinder(1e-9, 1.0, mu_r)
    expected = 4.0 * mu_r / (mu_r + 1.0) ** 2
    assert S == pytest.approx(expected, rel=1e-6)


def test_sphere_solid_limit_no_inside():
    mu_r = 100.0
    S = shielding_factor_sphere(1e-9, 1.0, mu_r)
    expected = 9.0 * mu_r / ((mu_r + 2.0) * (2.0 * mu_r + 1.0))
    assert S == pytest.approx(expected, rel=1e-6)


# ---------------------------------------------------------------------------
# Monotonicity in mu_r and in t/a
# ---------------------------------------------------------------------------


def test_cylinder_S_monotone_decreasing_in_mu_r():
    a, b = 0.9, 1.0
    Ss = [shielding_factor_cylinder(a, b, mu_r) for mu_r in [1.0, 10.0, 100.0, 1000.0]]
    assert all(Ss[i] > Ss[i + 1] for i in range(3))


def test_sphere_S_monotone_decreasing_in_mu_r():
    a, b = 0.9, 1.0
    Ss = [shielding_factor_sphere(a, b, mu_r) for mu_r in [1.0, 10.0, 100.0, 1000.0]]
    assert all(Ss[i] > Ss[i + 1] for i in range(3))


def test_cylinder_S_monotone_decreasing_in_thickness():
    # Thicker wall (smaller a/b) gives stronger shielding (smaller S).
    mu_r, b = 100.0, 1.0
    Ss = [shielding_factor_cylinder(a, b, mu_r) for a in [0.99, 0.95, 0.9, 0.5]]
    assert all(Ss[i] > Ss[i + 1] for i in range(3))


def test_argument_validation_cyl():
    with pytest.raises(ValueError):
        shielding_factor_cylinder(0.0, 1.0, 100.0)
    with pytest.raises(ValueError):
        shielding_factor_cylinder(1.0, 0.5, 100.0)  # a > b
    with pytest.raises(ValueError):
        shielding_factor_cylinder(0.5, 1.0, 0.5)    # mu_r < 1


def test_argument_validation_sphere():
    with pytest.raises(ValueError):
        shielding_factor_sphere(0.0, 1.0, 100.0)
    with pytest.raises(ValueError):
        shielding_factor_sphere(1.0, 0.5, 100.0)
    with pytest.raises(ValueError):
        shielding_factor_sphere(0.5, 1.0, 0.5)
