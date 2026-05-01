"""Tests for radia.analytical_formulas.ac_locus (Part 5, eq 29-37)."""

import math

import numpy as np
import pytest

from radia.analytical_formulas import ac_locus_axes, ac_locus_axes_batch


def test_zero_phasor():
    Bmax, Bmin = ac_locus_axes([0.0, 0.0, 0.0])
    assert Bmax == 0.0
    assert Bmin == 0.0


def test_linear_polarisation_along_x():
    Bmax, Bmin = ac_locus_axes([1.0, 0.0, 0.0])
    assert Bmax == pytest.approx(1.0, abs=1e-12)
    assert Bmin == pytest.approx(0.0, abs=1e-12)


def test_linear_polarisation_arbitrary_phase():
    # Same phase on each component -- still linear, length |B|.
    phase = np.exp(1j * 0.7)
    B = phase * np.array([0.6, -0.5, 0.4])
    Bmax, Bmin = ac_locus_axes(B)
    expected = np.linalg.norm([0.6, -0.5, 0.4])
    assert Bmax == pytest.approx(expected, abs=1e-12)
    assert Bmin == pytest.approx(0.0, abs=1e-12)


def test_circular_polarisation_in_xy_plane():
    Bmax, Bmin = ac_locus_axes([1.0, 1j, 0.0])
    assert Bmax == pytest.approx(1.0, abs=1e-12)
    assert Bmin == pytest.approx(1.0, abs=1e-12)


def test_elliptical_polarisation_in_xy_plane():
    Bmax, Bmin = ac_locus_axes([2.0, 1j, 0.0])
    assert Bmax == pytest.approx(2.0, abs=1e-12)
    assert Bmin == pytest.approx(1.0, abs=1e-12)


def test_axis_invariant_under_unitary_rotation():
    # Rotating B in physical space leaves the locus (an ellipse) congruent;
    # major / minor axes are unchanged.
    rng = np.random.default_rng(seed=1234)
    B = rng.standard_normal(3) + 1j * rng.standard_normal(3)
    Bmax_ref, Bmin_ref = ac_locus_axes(B)

    # Rotate by a random rotation matrix.
    theta, phi, psi = rng.uniform(0, 2 * np.pi, size=3)
    Rx = np.array([[1, 0, 0], [0, math.cos(theta), -math.sin(theta)],
                   [0, math.sin(theta), math.cos(theta)]])
    Ry = np.array([[math.cos(phi), 0, math.sin(phi)], [0, 1, 0],
                   [-math.sin(phi), 0, math.cos(phi)]])
    Rz = np.array([[math.cos(psi), -math.sin(psi), 0],
                   [math.sin(psi), math.cos(psi), 0], [0, 0, 1]])
    R = Rz @ Ry @ Rx
    Bmax_rot, Bmin_rot = ac_locus_axes(R @ B)
    assert Bmax_rot == pytest.approx(Bmax_ref, rel=1e-12)
    assert Bmin_rot == pytest.approx(Bmin_ref, rel=1e-12)


def test_axes_match_brute_force_time_sweep():
    # Compare against |b(t)| evaluated on a dense time grid.
    rng = np.random.default_rng(seed=42)
    B = rng.standard_normal(3) + 1j * rng.standard_normal(3)
    Bmax_an, Bmin_an = ac_locus_axes(B)

    omega_t = np.linspace(0.0, 2 * np.pi, 4001, endpoint=False)
    bt = np.real(B[None, :] * np.exp(1j * omega_t)[:, None])  # (N, 3)
    mag = np.linalg.norm(bt, axis=1)
    assert mag.max() == pytest.approx(Bmax_an, rel=1e-4)
    assert mag.min() == pytest.approx(Bmin_an, rel=1e-3, abs=1e-4)


def test_batch_matches_loop():
    rng = np.random.default_rng(seed=7)
    B = rng.standard_normal((20, 3)) + 1j * rng.standard_normal((20, 3))
    Bmax_b, Bmin_b = ac_locus_axes_batch(B)
    for k, (mx, mn) in enumerate(zip(Bmax_b, Bmin_b)):
        a, b = ac_locus_axes(B[k])
        assert mx == pytest.approx(a, rel=1e-12)
        assert mn == pytest.approx(b, rel=1e-12)


def test_batch_higher_rank():
    rng = np.random.default_rng(seed=2)
    B = rng.standard_normal((4, 5, 3)) + 1j * rng.standard_normal((4, 5, 3))
    Bmax, Bmin = ac_locus_axes_batch(B)
    assert Bmax.shape == (4, 5)
    assert Bmin.shape == (4, 5)
    # Spot-check one element.
    mx, mn = ac_locus_axes(B[2, 3])
    assert Bmax[2, 3] == pytest.approx(mx, rel=1e-12)


def test_batch_validates_shape():
    with pytest.raises(ValueError):
        ac_locus_axes_batch(np.zeros((4, 2), dtype=complex))


def test_invalid_shape():
    with pytest.raises(ValueError):
        ac_locus_axes([1.0, 0.0])
