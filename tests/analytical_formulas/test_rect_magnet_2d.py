"""Tests for radia.analytical_formulas.rect_magnet_2d (Part 2, eq 2-3)."""

import math

import numpy as np
import pytest

from radia.analytical_formulas import rect_magnet_2d_A, rect_magnet_2d_B
from radia.analytical_formulas.rect_magnet_2d import MU_0


# ---------------------------------------------------------------------------
# Slab limits (verifies M_x convention)
# ---------------------------------------------------------------------------


def test_long_bar_M_along_long_axis_gives_full_B():
    # Bar with a << b (long along y) and M along the long axis y:
    # Inside, demag along long axis -> 0, B_y -> mu_0 * M_y.
    a, b = 0.001, 1.0
    My = 1.0e6
    B = rect_magnet_2d_B([0.0, 0.0], a, b, M=(0.0, My))
    assert B[1] == pytest.approx(MU_0 * My, rel=1e-3)
    assert B[0] == pytest.approx(0.0, abs=1e-3 * MU_0 * My)


def test_long_bar_M_perpendicular_to_long_axis_gives_zero_B_inside():
    # Bar a << b, M along x (perpendicular to long axis): N_x -> 1, B_x -> 0.
    a, b = 0.001, 1.0
    Mx = 1.0e6
    B = rect_magnet_2d_B([0.0, 0.0], a, b, M=(Mx, 0.0))
    assert abs(B[0]) < 1e-3 * MU_0 * Mx
    assert abs(B[1]) < 1e-9


def test_square_bar_demag_one_half():
    # Square (a == b), M_x. By symmetry of cross-section the in-plane
    # demagnetization factor along x is 1/2.
    a = b = 0.05
    Mx = 1.0e6
    B = rect_magnet_2d_B([0.0, 0.0], a, b, M=(Mx, 0.0))
    expected = MU_0 * Mx * 0.5
    assert B[0] == pytest.approx(expected, rel=1e-12)
    assert B[1] == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Symmetry
# ---------------------------------------------------------------------------


def test_M_x_bar_B_x_odd_in_y():
    a, b = 0.04, 0.06
    Mx = 1.0e6
    B_pos = rect_magnet_2d_B([0.5, +0.3], a, b, M=(Mx, 0.0))
    B_neg = rect_magnet_2d_B([0.5, -0.3], a, b, M=(Mx, 0.0))
    # B_x is even in y; B_y is odd in y for M_x-bar (curl symmetry).
    # Wait -- A_z is odd in y (top/bottom currents flip sign), so:
    # B_x = dA/dy is even in y, B_y = -dA/dx is odd in y.
    assert B_pos[0] == pytest.approx(B_neg[0], rel=1e-12)
    assert B_pos[1] == pytest.approx(-B_neg[1], rel=1e-12)


def test_M_y_bar_B_x_odd_in_x():
    a, b = 0.04, 0.06
    My = 1.0e6
    B_pos = rect_magnet_2d_B([+0.3, 0.5], a, b, M=(0.0, My))
    B_neg = rect_magnet_2d_B([-0.3, 0.5], a, b, M=(0.0, My))
    # By the mirror symmetry argument: A_z(x, y; My) = -A_z(-x, y; My).
    # So B_x = dA_z/dy: A_z(x,y) = -A_z(-x, y) -> dA/dy at (x, y) = -dA/dy at (-x, y).
    # Hence B_x is odd in x, B_y is even in x.
    assert B_pos[0] == pytest.approx(-B_neg[0], rel=1e-12)
    assert B_pos[1] == pytest.approx(B_neg[1], rel=1e-12)


# ---------------------------------------------------------------------------
# Far-field 2D dipole behaviour
# ---------------------------------------------------------------------------


def test_far_field_falls_off_like_inverse_r_squared():
    # 2D point dipole field magnitude ~ 1 / r**2.
    a, b = 0.01, 0.005
    Mx = 1.0e6
    P1 = np.array([1.0, 0.5])
    P2 = 2.0 * P1
    B1 = rect_magnet_2d_B(P1, a, b, M=(Mx, 0.0))
    B2 = rect_magnet_2d_B(P2, a, b, M=(Mx, 0.0))
    ratio = np.linalg.norm(B2) / np.linalg.norm(B1)
    # |B| at 2r should be ~ 1/4 of that at r.
    assert ratio == pytest.approx(0.25, rel=1e-2)


def test_far_field_dipole_amplitude():
    # 2D magnetostatic dipole far-field magnitude:
    # |B| = mu_0 * m / (2 pi r**2), with m = M_x * 4 a b (per unit length).
    a, b = 0.005, 0.0025
    Mx = 1.0e6
    m = Mx * (4.0 * a * b)
    r = 5.0
    P = np.array([r, 0.0])
    B = rect_magnet_2d_B(P, a, b, M=(Mx, 0.0))
    # On axis y = 0, M along x: dipole pattern -> B_x = +mu_0 m / (2 pi r**2)
    expected_Bx = MU_0 * m / (2.0 * math.pi * r * r)
    assert B[0] == pytest.approx(expected_Bx, rel=2e-3)


# ---------------------------------------------------------------------------
# Curl consistency: B_x = dA/dy, B_y = -dA/dx
# ---------------------------------------------------------------------------


def test_curl_consistency_at_external_point():
    a, b = 0.04, 0.07
    Mx, My = 0.7e6, -0.3e6
    P = np.array([0.10, 0.05])
    h = 1e-5
    A0 = rect_magnet_2d_A(P, a, b, M=(Mx, My))
    Apx = rect_magnet_2d_A(P + [h, 0.0], a, b, M=(Mx, My))
    Amx = rect_magnet_2d_A(P + [-h, 0.0], a, b, M=(Mx, My))
    Apy = rect_magnet_2d_A(P + [0.0, h], a, b, M=(Mx, My))
    Amy = rect_magnet_2d_A(P + [0.0, -h], a, b, M=(Mx, My))
    Bx_fd = (Apy - Amy) / (2 * h)
    By_fd = -(Apx - Amx) / (2 * h)
    B = rect_magnet_2d_B(P, a, b, M=(Mx, My))
    assert B[0] == pytest.approx(Bx_fd, rel=1e-5, abs=1e-9)
    assert B[1] == pytest.approx(By_fd, rel=1e-5, abs=1e-9)


# ---------------------------------------------------------------------------
# Linearity in M
# ---------------------------------------------------------------------------


def test_linearity_in_M():
    a, b = 0.02, 0.03
    P = np.array([0.05, 0.04])
    Bx_only = rect_magnet_2d_B(P, a, b, M=(1.0e6, 0.0))
    By_only = rect_magnet_2d_B(P, a, b, M=(0.0, -2.0e6))
    Bxy = rect_magnet_2d_B(P, a, b, M=(1.0e6, -2.0e6))
    np.testing.assert_allclose(Bxy, Bx_only + By_only, rtol=1e-12)


# ---------------------------------------------------------------------------
# Vectorisation
# ---------------------------------------------------------------------------


def test_batch_evaluation():
    a, b = 0.02, 0.03
    points = np.array([[0.05, 0.04], [0.1, 0.0], [-0.05, -0.05]])
    M = (1.0e6, 0.0)
    B_batch = rect_magnet_2d_B(points, a, b, M=M)
    for k in range(len(points)):
        B_single = rect_magnet_2d_B(points[k], a, b, M=M)
        np.testing.assert_allclose(B_batch[k], B_single, rtol=1e-12)
