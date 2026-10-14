"""Tests for radia.analytical_formulas.plate_eddy (Part 1, eq 26-27)."""

import numpy as np
import pytest

from radia.analytical_formulas import plate_eddy_J, plate_eddy_T


# ---------------------------------------------------------------------------
# Boundary conditions on T_z and J . n
# ---------------------------------------------------------------------------


def test_T_zero_on_x_edges_exact():
    # cos((2n+1) pi a / (2 a)) = cos((2n+1) pi / 2) = 0 exactly, and
    # x**2 - a**2 = 0 at x = +- a, so T_z(+- a, y) = 0 exactly.
    a, b = 0.1, 0.05
    sigma, Bdot = 5e7, 1.0
    y_grid = np.linspace(-b, b, 11)
    for y in y_grid:
        for x in (+a, -a):
            T = plate_eddy_T([x, y], a, b, sigma, Bdot, n_terms=50)
            assert abs(T) < 1e-10 * sigma * Bdot * a * a


def test_T_zero_on_y_edges_with_enough_terms():
    a, b = 0.1, 0.05
    sigma, Bdot = 5e7, 1.0
    x_grid = np.linspace(-a, a, 7)
    scale = sigma * abs(Bdot) * a * a
    for x in x_grid:
        for y in (+b, -b):
            T = plate_eddy_T([x, y], a, b, sigma, Bdot, n_terms=500)
            assert abs(T) < 1e-3 * scale


def test_Jx_vanishes_on_x_edges():
    a, b = 0.1, 0.05
    sigma, Bdot = 5e7, 1.0
    y_grid = np.linspace(-0.8 * b, 0.8 * b, 5)
    for y in y_grid:
        Jp = plate_eddy_J([+a, y], a, b, sigma, Bdot, n_terms=50)
        Jm = plate_eddy_J([-a, y], a, b, sigma, Bdot, n_terms=50)
        # At x = +-a, cos((2n+1) pi / 2) = 0 -> series part of J_x is exactly zero.
        assert abs(Jp[0]) < 1e-6
        assert abs(Jm[0]) < 1e-6


def test_Jy_vanishes_on_y_edges_in_series_limit():
    a, b = 0.1, 0.05
    sigma, Bdot = 5e7, 1.0
    x_grid = np.linspace(-0.8 * a, 0.8 * a, 5)
    scale = sigma * abs(Bdot) * a
    for x in x_grid:
        Jp = plate_eddy_J([x, +b], a, b, sigma, Bdot, n_terms=2000)
        Jm = plate_eddy_J([x, -b], a, b, sigma, Bdot, n_terms=2000)
        assert abs(Jp[1]) < 1e-3 * scale
        assert abs(Jm[1]) < 1e-3 * scale


# ---------------------------------------------------------------------------
# Symmetry
# ---------------------------------------------------------------------------


def test_T_even_in_x_and_y():
    a, b = 0.1, 0.07
    sigma, Bdot = 5e7, 1.0
    Tpp = plate_eddy_T([0.04, 0.03], a, b, sigma, Bdot, n_terms=80)
    Tmp = plate_eddy_T([-0.04, 0.03], a, b, sigma, Bdot, n_terms=80)
    Tpm = plate_eddy_T([0.04, -0.03], a, b, sigma, Bdot, n_terms=80)
    Tmm = plate_eddy_T([-0.04, -0.03], a, b, sigma, Bdot, n_terms=80)
    assert Tpp == pytest.approx(Tmp, rel=1e-12, abs=1e-12)
    assert Tpp == pytest.approx(Tpm, rel=1e-12, abs=1e-12)
    assert Tpp == pytest.approx(Tmm, rel=1e-12, abs=1e-12)


def test_Jx_odd_in_y_even_in_x():
    a, b = 0.1, 0.07
    sigma, Bdot = 5e7, 1.0
    Jpp = plate_eddy_J([0.04, 0.03], a, b, sigma, Bdot, n_terms=80)
    Jmp = plate_eddy_J([-0.04, 0.03], a, b, sigma, Bdot, n_terms=80)
    Jpm = plate_eddy_J([0.04, -0.03], a, b, sigma, Bdot, n_terms=80)
    # J_x even in x, odd in y
    assert Jpp[0] == pytest.approx(Jmp[0], rel=1e-12)
    assert Jpp[0] == pytest.approx(-Jpm[0], rel=1e-12)
    # J_y odd in x, even in y
    assert Jpp[1] == pytest.approx(-Jmp[1], rel=1e-12)
    assert Jpp[1] == pytest.approx(Jpm[1], rel=1e-12)


# ---------------------------------------------------------------------------
# Linearity in (sigma, Bdot)
# ---------------------------------------------------------------------------


def test_linearity_in_sigma_and_Bdot():
    a, b = 0.1, 0.05
    P = [0.04, 0.03]
    T1 = plate_eddy_T(P, a, b, 5.0e7, 1.0, n_terms=100)
    T2 = plate_eddy_T(P, a, b, 1.5e7, 2.0, n_terms=100)
    # Doubling sigma * Bdot doubles T_z.
    T_double = plate_eddy_T(P, a, b, 1.0e8, 1.0, n_terms=100)
    assert T_double == pytest.approx(2.0 * T1, rel=1e-12)
    # T_2 / T_1 should be (1.5e7 * 2.0) / (5.0e7 * 1.0) = 0.6
    assert T2 / T1 == pytest.approx(0.6, rel=1e-12)


# ---------------------------------------------------------------------------
# J = curl(T e_z) verification by finite differences
# ---------------------------------------------------------------------------


def test_J_consistent_with_curl_T():
    a, b = 0.1, 0.05
    sigma, Bdot = 5e7, 1.0
    P = np.array([0.04, 0.02])
    h = 1e-4
    Tpx = plate_eddy_T(P + [h, 0.0], a, b, sigma, Bdot, n_terms=200)
    Tmx = plate_eddy_T(P + [-h, 0.0], a, b, sigma, Bdot, n_terms=200)
    Tpy = plate_eddy_T(P + [0.0, h], a, b, sigma, Bdot, n_terms=200)
    Tmy = plate_eddy_T(P + [0.0, -h], a, b, sigma, Bdot, n_terms=200)
    Jx_fd = (Tpy - Tmy) / (2 * h)
    Jy_fd = -(Tpx - Tmx) / (2 * h)
    J = plate_eddy_J(P, a, b, sigma, Bdot, n_terms=200)
    # FD has O(h**2) truncation; tolerate ~1e-3 rel.
    assert J[0] == pytest.approx(Jx_fd, rel=2e-3)
    assert J[1] == pytest.approx(Jy_fd, rel=2e-3)


def test_div_J_zero_at_interior():
    a, b = 0.1, 0.05
    sigma, Bdot = 5e7, 1.0
    P = np.array([0.04, 0.02])
    # FD step small enough that O(h**2 max|J'''|) truncation is well below
    # |J| / Lscale**2; a typical interior |J''' / J| ~ 1 / Lscale**2 with
    # Lscale ~ a / 5 for first-mode terms, so 1e-4 keeps the truncation
    # error around |J| * 1e-6.
    h = 1e-5
    Jpx = plate_eddy_J(P + [h, 0.0], a, b, sigma, Bdot, n_terms=200)
    Jmx = plate_eddy_J(P + [-h, 0.0], a, b, sigma, Bdot, n_terms=200)
    Jpy = plate_eddy_J(P + [0.0, h], a, b, sigma, Bdot, n_terms=200)
    Jmy = plate_eddy_J(P + [0.0, -h], a, b, sigma, Bdot, n_terms=200)
    div = (Jpx[0] - Jmx[0]) / (2 * h) + (Jpy[1] - Jmy[1]) / (2 * h)
    Jmag = np.linalg.norm(plate_eddy_J(P, a, b, sigma, Bdot, n_terms=200))
    # Mathematically zero; FD truncation gives O(h**2 * max|J'''|).
    assert abs(div) / Jmag < 1e-5


# ---------------------------------------------------------------------------
# Long-strip limit
# ---------------------------------------------------------------------------


def test_long_strip_limit_recovers_simple_solution():
    # b -> infinity: cosh ratio -> 0, so T_z -> (sigma Bdot / 8) (x**2 - a**2)
    # and J_y -> -(sigma Bdot / 4) x, J_x -> 0.
    a, b = 0.05, 5.0       # very long plate
    sigma, Bdot = 5e7, 1.0
    T = plate_eddy_T([0.02, 0.0], a, b, sigma, Bdot, n_terms=20)
    expected = (sigma * Bdot / 8.0) * (0.02 ** 2 - a ** 2)
    assert T == pytest.approx(expected, rel=1e-6)
    J = plate_eddy_J([0.02, 0.0], a, b, sigma, Bdot, n_terms=20)
    assert abs(J[0]) < 1e-6 * abs(J[1])
    assert J[1] == pytest.approx(-(sigma * Bdot / 4.0) * 0.02, rel=1e-6)


# ---------------------------------------------------------------------------
# Vectorisation
# ---------------------------------------------------------------------------


def test_batch_evaluation():
    a, b = 0.1, 0.05
    sigma, Bdot = 5e7, 1.0
    points = np.array([[0.0, 0.0], [0.05, 0.02], [-0.03, -0.01]])
    T_batch = plate_eddy_T(points, a, b, sigma, Bdot, n_terms=80)
    J_batch = plate_eddy_J(points, a, b, sigma, Bdot, n_terms=80)
    for k in range(len(points)):
        T_one = plate_eddy_T(points[k], a, b, sigma, Bdot, n_terms=80)
        J_one = plate_eddy_J(points[k], a, b, sigma, Bdot, n_terms=80)
        assert T_batch[k] == pytest.approx(T_one, rel=1e-12)
        np.testing.assert_allclose(J_batch[k], J_one, rtol=1e-12)
