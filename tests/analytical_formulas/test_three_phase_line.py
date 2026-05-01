"""Tests for radia.analytical_formulas.three_phase_line (Part 4 §5, Part 5 §3)."""

import math

import numpy as np
import pytest

from radia.analytical_formulas import (
    ac_locus_axes,
    balanced_six_phase_currents,
    balanced_three_phase_currents,
    field_xy,
    helical_far_field_amplitude,
    helical_near_field_amplitude,
    hexagon_positions,
    planar_far_field_amplitude,
    planar_positions,
    triangle_far_field_amplitude,
    triangle_positions,
    vector_potential_z,
)


# ---------------------------------------------------------------------------
# Direct evaluator: check against single-wire Biot-Savart
# ---------------------------------------------------------------------------


def test_single_wire_field_matches_textbook():
    """B = mu_0 I / (2 pi r) for a single straight wire."""
    I = 100.0
    positions = np.array([[0.0, 0.0]])
    obs = np.array([1.0, 0.0])
    B = field_xy(positions, [I], obs)
    # B should be perpendicular to the radial direction; for I in +z at origin
    # observed at (+1, 0), B is in +y.
    assert B[0] == pytest.approx(0.0, abs=1e-15)
    assert B[1] == pytest.approx(4.0e-7 * math.pi * I / (2.0 * math.pi * 1.0),
                                 rel=1e-14)


def test_balanced_three_phase_sum_is_zero():
    """The sum of balanced three-phase currents is zero."""
    I = balanced_three_phase_currents(100.0)
    assert abs(np.sum(I)) < 1e-10


def test_balanced_six_phase_sum_is_zero():
    I = balanced_six_phase_currents(100.0)
    assert abs(np.sum(I)) < 1e-10


# ---------------------------------------------------------------------------
# Triangle far-field: closed form vs direct evaluation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("r_over_a", [200, 1000, 5000])
def test_triangle_far_field_closed_form(r_over_a):
    """At r >> a the rigorous |B|_max from field_xy + ac_locus_axes
    matches ``triangle_far_field_amplitude`` to leading order, with
    convergence rate ``a / r`` (next multipole is n = 4 with
    coefficient O((a/r)**3) but cross-talk through the lower term gives
    ``a/r`` envelope).
    """
    a = 0.05
    I_peak = 100.0
    r = r_over_a * a
    positions = triangle_positions(a)
    currents = balanced_three_phase_currents(I_peak)
    B = field_xy(positions, currents, [r, 0.0])
    Bmax, Bmin = ac_locus_axes([B[0], B[1], 0.0])
    formula = triangle_far_field_amplitude(r, I_peak, a)
    # Empirical fit: |Bmax/formula - 1| ~ 1.0 / (r/a)
    rel_tol = 2.0 / r_over_a
    assert Bmax == pytest.approx(formula, rel=rel_tol)
    # Far field is circularly polarised: 1 - Bmin / Bmax -> 0 like 2 / (r/a).
    assert 1.0 - Bmin / Bmax < 5.0 / r_over_a


def test_triangle_far_field_decay_one_over_r_squared():
    a = 0.05
    I_peak = 100.0
    Bs = []
    # Use r/a = 100, 200, 400 so that the 1/r**2 leading term clearly dominates.
    for r in (5.0, 10.0, 20.0):
        positions = triangle_positions(a)
        currents = balanced_three_phase_currents(I_peak)
        B = field_xy(positions, currents, [r, 0.0])
        Bmax, _ = ac_locus_axes([B[0], B[1], 0.0])
        Bs.append(Bmax)
    assert Bs[1] / Bs[0] == pytest.approx(0.25, rel=1e-2)
    assert Bs[2] / Bs[1] == pytest.approx(0.25, rel=1e-2)


# ---------------------------------------------------------------------------
# Planar far-field: closed form vs direct evaluation, sqrt(3) coefficient
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("r_over_d", [50, 200, 1000])
def test_planar_far_field_closed_form(r_over_d):
    d = 0.05
    I_peak = 100.0
    r = r_over_d * d
    positions = planar_positions(d)
    currents = balanced_three_phase_currents(I_peak)
    B = field_xy(positions, currents, [r, 0.0])
    Bmax, Bmin = ac_locus_axes([B[0], B[1], 0.0])
    formula = planar_far_field_amplitude(r, I_peak, d)
    rel_tol = 5.0 / r_over_d
    assert Bmax == pytest.approx(formula, rel=rel_tol)
    # Planar far-field is linearly polarised (Bmin -> 0).
    assert Bmin / Bmax < 5.0 / r_over_d


def test_planar_decay_one_over_r_squared():
    d = 0.05
    I_peak = 100.0
    Bs = [
        ac_locus_axes(
            [*field_xy(planar_positions(d), balanced_three_phase_currents(I_peak),
                       [r, 0.0]), 0.0]
        )[0] for r in (5.0, 10.0, 20.0)
    ]
    assert Bs[1] / Bs[0] == pytest.approx(0.25, rel=1e-2)
    assert Bs[2] / Bs[1] == pytest.approx(0.25, rel=1e-2)


# ---------------------------------------------------------------------------
# Hexagon naive arrangement: regression-test that it is dipole-dominant
# ---------------------------------------------------------------------------


def test_hexagon_naive_arrangement_is_dipolar():
    """Sanity: ``hexagon_positions + balanced_six_phase_currents`` is
    a sum of three anti-phase 3-phase pairs; the dipoles do *not*
    cancel and the far-field decays as 1 / r**2.
    """
    a = 0.5
    I_peak = 100.0
    positions = hexagon_positions(a)
    currents = balanced_six_phase_currents(I_peak)
    Bmags = []
    for r in (10.0, 100.0, 1000.0):
        B = field_xy(positions, currents, [r, 0.0])
        Bmax, _ = ac_locus_axes([B[0], B[1], 0.0])
        Bmags.append(Bmax)
    # B(10 r) / B(r) ~ 1/100 (dipole 1/r**2), not 1/1000 (1/r**3).
    assert Bmags[1] / Bmags[0] == pytest.approx(0.01, rel=1e-3)
    assert Bmags[2] / Bmags[1] == pytest.approx(0.01, rel=1e-3)


# ---------------------------------------------------------------------------
# Helical asymptotes
# ---------------------------------------------------------------------------


def test_helical_near_field_matches_triangle():
    """In the regime r << p the helical field is approximated by the
    static triangle field.
    """
    r, p, a, I = 0.01, 1.0, 0.001, 100.0
    near = helical_near_field_amplitude(r, I, a, p)
    tri = triangle_far_field_amplitude(r, I, a)
    assert near == pytest.approx(tri, rel=1e-12)


def test_helical_far_field_decays_exponentially():
    """exp(-2 pi r / p) is the dominant factor at r >> p."""
    p, a, I = 0.5, 0.05, 100.0
    r1, r2 = 1.0, 1.5
    B1 = helical_far_field_amplitude(r1, I, a, p)
    B2 = helical_far_field_amplitude(r2, I, a, p)
    expected_ratio = (
        math.sqrt(r1 / r2) * math.exp(-2.0 * math.pi * (r2 - r1) / p)
    )
    assert B2 / B1 == pytest.approx(expected_ratio, rel=1e-12)


def test_helical_far_field_validates_pitch():
    with pytest.raises(ValueError):
        helical_far_field_amplitude(1.0, 100.0, 0.05, p=0.0)


# ---------------------------------------------------------------------------
# Vector-potential consistency: B = curl A
# ---------------------------------------------------------------------------


def test_vector_potential_z_curl_equals_field_xy():
    a = 0.05
    I_peak = 100.0
    positions = triangle_positions(a)
    currents = balanced_three_phase_currents(I_peak)
    P = np.array([0.5, 0.3])
    h = 1e-5
    A = vector_potential_z(positions, currents, P)
    Apx = vector_potential_z(positions, currents, P + [h, 0])
    Amx = vector_potential_z(positions, currents, P + [-h, 0])
    Apy = vector_potential_z(positions, currents, P + [0, h])
    Amy = vector_potential_z(positions, currents, P + [0, -h])
    Bx_fd = (Apy - Amy) / (2 * h)
    By_fd = -(Apx - Amx) / (2 * h)
    B = field_xy(positions, currents, P)
    assert abs(B[0] - Bx_fd) / abs(B[0]) < 1e-5
    assert abs(B[1] - By_fd) / abs(B[1]) < 1e-5
