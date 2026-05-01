"""Tests for radia.analytical_formulas.solenoid_central (Part 4 §4)."""

import math

import numpy as np
import pytest

from radia.analytical_formulas import (
    fabri_F,
    solenoid_axial_field,
    solenoid_central_field,
)
from radia.analytical_formulas.solenoid_central import MU_0


# ---------------------------------------------------------------------------
# Limit checks
# ---------------------------------------------------------------------------


def test_infinite_solenoid_limit():
    """Long solenoid (b -> infinity): B_0 -> mu_0 J (a_2 - a_1)."""
    a1, a2 = 0.05, 0.1
    J = 1.0e7
    B_inf = MU_0 * J * (a2 - a1)
    B_long = solenoid_central_field(a1, a2, b=100.0, J=J)
    assert B_long == pytest.approx(B_inf, rel=1e-3)


def test_thin_radial_layer_limit():
    """Thin radial layer (alpha -> 1) gives F -> 0 and hence B_0 -> 0."""
    a1 = 0.1
    F = fabri_F(alpha=1.0001, beta=1.0)
    assert 0.0 < F < 1e-3
    B = solenoid_central_field(a1, a1 * 1.0001, b=0.1, J=1.0e7)
    assert 0.0 < B < 1e-3 * MU_0 * 1.0e7 * a1


def test_short_pancake_limit():
    """Short pancake (b -> 0): F -> beta * ln(alpha)."""
    alpha = 2.5
    for beta in (1e-2, 1e-3):
        F = fabri_F(alpha, beta)
        F_pancake = beta * math.log(alpha)
        assert F == pytest.approx(F_pancake, rel=1e-2)


def test_axial_at_origin_matches_central_field():
    """axial_field(z=0, ...) == central_field(...)."""
    a1, a2, b, J = 0.05, 0.1, 0.1, 1.0e7
    B_axis = float(solenoid_axial_field(0.0, a1, a2, b, J))
    B_centre = solenoid_central_field(a1, a2, b, J)
    assert B_axis == pytest.approx(B_centre, rel=1e-13)


def test_axial_field_symmetry():
    """B_z(z) is even in z."""
    zs = np.array([-0.05, -0.02, 0.0, 0.02, 0.05])
    B = solenoid_axial_field(zs, a1=0.05, a2=0.1, b=0.1, J=1.0e7)
    assert B[0] == pytest.approx(B[-1], rel=1e-12)
    assert B[1] == pytest.approx(B[-2], rel=1e-12)


def test_axial_field_far_outside_decays_like_dipole():
    """Far from the coil, B_z ~ 1 / z**3 (axial dipole)."""
    a1, a2, b, J = 0.05, 0.1, 0.1, 1.0e7
    z1, z2 = 1.0, 2.0
    B1 = float(solenoid_axial_field(z1, a1, a2, b, J))
    B2 = float(solenoid_axial_field(z2, a1, a2, b, J))
    # B(2 z) / B(z) ~ (z / 2 z)**3 = 1 / 8
    assert B2 / B1 == pytest.approx(1.0 / 8.0, rel=2e-2)


def test_central_field_linear_in_J():
    a1, a2, b = 0.05, 0.1, 0.1
    B1 = solenoid_central_field(a1, a2, b, J=1.0e7)
    B2 = solenoid_central_field(a1, a2, b, J=3.0e7)
    assert B2 / B1 == pytest.approx(3.0, rel=1e-12)


def test_argument_validation():
    with pytest.raises(ValueError):
        fabri_F(alpha=0.5, beta=1.0)        # alpha < 1
    with pytest.raises(ValueError):
        fabri_F(alpha=2.0, beta=-1.0)       # beta < 0
    with pytest.raises(ValueError):
        solenoid_central_field(0.0, 0.1, 0.1, 1.0)
    with pytest.raises(ValueError):
        solenoid_central_field(0.05, 0.04, 0.1, 1.0)   # a2 < a1
