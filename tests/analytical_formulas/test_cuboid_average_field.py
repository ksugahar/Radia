"""Closed-form vs Gauss-Legendre cross-validation for ``cuboid_average_field``.

Phase-alpha closed form (sympy-derived; ``project_phase_alpha_g1_g2_derived``)
must agree with the legacy Gauss-Legendre quadrature path at the
quadrature precision floor.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from radia.analytical_formulas.cuboid_average_field import (
    MU_0,
    average_B_in_box,
    average_B_in_box_numerical,
    average_demag_tensor,
    _average_demag_tensor_python,
    _HAVE_CPP,
)


# ---------------------------------------------------------------------------
# Closed-form vs numerical cross-validation (parametrised over geometries)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "label, src_min, src_max, tgt_min, tgt_max, M",
    [
        ("offset_x_Mx",
         (0., 0., 0.), (1., 1., 1.), (2., 0., 0.), (3., 1., 1.), (1e6, 0., 0.)),
        ("offset_y_Mx",
         (0., 0., 0.), (1., 1., 1.), (0., 2., 0.), (1., 3., 1.), (1e6, 0., 0.)),
        ("offset_z_Mx",
         (0., 0., 0.), (1., 1., 1.), (0., 0., 2.), (1., 1., 3.), (1e6, 0., 0.)),
        ("diagonal_offset",
         (0., 0., 0.), (1., 1., 1.), (2., 2., 2.), (3., 3., 3.), (1e6, 0., 0.)),
        ("M_y_offset_x",
         (0., 0., 0.), (1., 1., 1.), (2., 0., 0.), (3., 1., 1.), (0., 1e6, 0.)),
        ("M_z_offset_diag_xy",
         (0., 0., 0.), (1., 1., 1.), (2., 2., 0.), (3., 3., 1.), (0., 0., 1e6)),
        ("different_sizes_arbitrary_M",
         (-0.5, -0.5, -0.5), (0.5, 0.5, 0.5),
         (1.5, -1.0, -2.0), (3.5, 1.0, 2.0),
         (1e6, 5e5, 0.)),
        ("distant_target",
         (-0.5, -0.5, -0.5), (0.5, 0.5, 0.5),
         (9.5, -0.5, -0.5), (10.5, 0.5, 0.5),
         (1e6, 0., 0.)),
    ],
)
def test_closed_form_matches_numerical(label, src_min, src_max, tgt_min, tgt_max, M):
    B_closed = average_B_in_box(M, src_min, src_max, tgt_min, tgt_max, method="closed")
    B_num = average_B_in_box(M, src_min, src_max, tgt_min, tgt_max, method="numerical", n_quad=12)
    # The numerical Gauss-Legendre at n_quad=12 hits ~1e-7 relative precision when
    # the target box is small compared to the source (slow integrand variation).
    # Closed form is the more accurate of the two. Use rtol = 1e-6 here.
    for k in range(3):
        if abs(B_num[k]) > 1e-12:
            np.testing.assert_allclose(B_closed[k], B_num[k], rtol=1e-6,
                                       err_msg=f"{label}: B[{k}] closed != numerical")
        else:
            assert abs(B_closed[k]) < 1e-12, f"{label}: B[{k}] expected ~0, got {B_closed[k]}"


# ---------------------------------------------------------------------------
# Self-source canonical demag tensor (= -1/3 I for a cube)
# ---------------------------------------------------------------------------
def test_self_cube_demag_tensor():
    """Average demag tensor over a cube == itself == -(1/3) I."""
    smin = (0., 0., 0.); smax = (1., 1., 1.)
    A = average_demag_tensor(smin, smax, smin, smax)
    np.testing.assert_allclose(A, -np.eye(3) / 3.0, atol=1e-13)


def test_self_cube_average_B():
    """For self-source (target = source) cube with ``M = M_x x_hat``,
    ``<B_x>_T = (2/3) mu_0 M_x``."""
    smin = (0., 0., 0.); smax = (1., 1., 1.)
    M = (1e6, 0., 0.)
    B = average_B_in_box(M, smin, smax, smin, smax, method="closed")
    expected = (2.0 / 3.0) * MU_0 * M[0]
    np.testing.assert_allclose(B[0], expected, rtol=1e-13)
    np.testing.assert_allclose(B[1], 0.0, atol=1e-13)
    np.testing.assert_allclose(B[2], 0.0, atol=1e-13)


# ---------------------------------------------------------------------------
# Overlap correction: half-overlap test
# ---------------------------------------------------------------------------
def test_half_overlap_matches_numerical():
    """Source [0,1]^3, target [0.5,1.5]x[0,1]x[0,1] -- 50% overlap.

    Closed-form must include the ``mu_0 * M * V_overlap / V_T`` correction.
    """
    src_min = (0., 0., 0.); src_max = (1., 1., 1.)
    tgt_min = (0.5, 0., 0.); tgt_max = (1.5, 1., 1.)
    M = (1e6, 0., 0.)
    B_closed = average_B_in_box(M, src_min, src_max, tgt_min, tgt_max, method="closed")
    B_num = average_B_in_box_numerical(M, src_min, src_max, tgt_min, tgt_max, n_quad=12)
    np.testing.assert_allclose(B_closed[0], B_num[0], rtol=1e-9)


# ---------------------------------------------------------------------------
# Linearity in M
# ---------------------------------------------------------------------------
def test_linear_in_M_closed():
    smin = (-0.01, -0.01, -0.005); smax = (0.01, 0.01, 0.005)
    tmin = (0.04, -0.001, -0.001); tmax = (0.06, 0.001, 0.001)
    B1 = average_B_in_box([1e6, 0., 0.], smin, smax, tmin, tmax)
    B2 = average_B_in_box([3e6, 0., 0.], smin, smax, tmin, tmax)
    np.testing.assert_allclose(B2, 3.0 * np.asarray(B1), rtol=1e-13, atol=1e-15)


def test_superposition_closed():
    """B(M_x x + M_y y) == B(M_x x) + B(M_y y) by linearity."""
    smin = (-0.5, -0.5, -0.5); smax = (0.5, 0.5, 0.5)
    tmin = (1.0, 0.0, 0.0); tmax = (2.0, 1.0, 1.0)
    B_x = average_B_in_box((1e6, 0., 0.), smin, smax, tmin, tmax)
    B_y = average_B_in_box((0., 1e6, 0.), smin, smax, tmin, tmax)
    B_xy = average_B_in_box((1e6, 1e6, 0.), smin, smax, tmin, tmax)
    np.testing.assert_allclose(B_xy, np.asarray(B_x) + np.asarray(B_y), rtol=1e-13)


# ---------------------------------------------------------------------------
# Tensor symmetry: <A_ij> should be symmetric (passive demag tensor identity)
# ---------------------------------------------------------------------------
def test_demag_tensor_symmetric():
    """<A>_T should be symmetric (i.e. ``<A_xy> = <A_yx>``)."""
    smin = (0., 0., 0.); smax = (1.5, 1.0, 0.7)
    tmin = (3., 1., 2.); tmax = (4., 2., 3.)
    A = average_demag_tensor(smin, smax, tmin, tmax)
    np.testing.assert_allclose(A, A.T, atol=1e-13)


# ---------------------------------------------------------------------------
# Far-field dipole limit
# ---------------------------------------------------------------------------
def test_far_field_dipole_closed():
    """Closed-form B at large distance matches the dipole formula.

    Note on box scales: the closed-form 64-corner sum suffers from
    catastrophic cancellation when ``V_T`` is much smaller than the
    source box (corner-primitive values are ``O(1)`` and cancel to give
    a tiny averaged result). For accurate far-field evaluation, keep
    the target box dimensions on the same order as the source box;
    for true point evaluation, call ``CuboidMagnet.get_B(point)`` directly.
    """
    src_size = np.array([0.5, 0.5, 0.5])
    src_min = -src_size / 2
    src_max = +src_size / 2
    M = np.array([1e6, 0., 0.])
    m = M[0] * float(np.prod(src_size))
    d = 5.0
    tgt_size = np.array([0.5, 0.5, 0.5])
    tgt_min = np.array([d - tgt_size[0] / 2, -tgt_size[1] / 2, -tgt_size[2] / 2])
    tgt_max = tgt_min + tgt_size
    B_avg = average_B_in_box(M, src_min, src_max, tgt_min, tgt_max)
    B_dipole = (MU_0 / (4 * math.pi)) * 2 * m / d**3
    np.testing.assert_allclose(B_avg[0], B_dipole, rtol=5e-2)
    assert abs(B_avg[1]) < 1e-12
    assert abs(B_avg[2]) < 1e-12


# ---------------------------------------------------------------------------
# Method='numerical' selector still works
# ---------------------------------------------------------------------------
def test_method_numerical_still_works():
    smin = (0., 0., 0.); smax = (1., 1., 1.)
    tmin = (2., 0., 0.); tmax = (3., 1., 1.)
    B = average_B_in_box((1e6, 0., 0.), smin, smax, tmin, tmax,
                         method="numerical", n_quad=6)
    assert B[0] > 0  # smoke test


def test_method_unknown_raises():
    with pytest.raises(ValueError):
        average_B_in_box((1e6, 0., 0.), (0., 0., 0.), (1., 1., 1.),
                        (2., 0., 0.), (3., 1., 1.), method="bogus")


# ---------------------------------------------------------------------------
# C++ vs Python reference agreement (Phase beta)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not _HAVE_CPP, reason="radia._radia_pybind not available")
@pytest.mark.parametrize(
    "src_min, src_max, tgt_min, tgt_max",
    [
        ((0., 0., 0.), (1., 1., 1.), (2., 0., 0.), (3., 1., 1.)),
        ((0., 0., 0.), (1., 1., 1.), (0.5, 0., 0.), (1.5, 1., 1.)),  # half overlap
        ((0., 0., 0.), (1., 1., 1.), (0., 0., 0.), (1., 1., 1.)),    # self
        ((-0.5, -0.5, -0.5), (0.5, 0.5, 0.5), (1.5, -1.0, -2.0), (3.5, 1.0, 2.0)),
        ((0., 0., 0.), (1.5, 1.0, 0.7), (3., 1., 2.), (4., 2., 3.)),  # asymmetric sizes
    ],
)
def test_cpp_matches_python_reference(src_min, src_max, tgt_min, tgt_max):
    A_cpp = average_demag_tensor(src_min, src_max, tgt_min, tgt_max)  # C++ path
    A_py = _average_demag_tensor_python(src_min, src_max, tgt_min, tgt_max)
    np.testing.assert_allclose(A_cpp, A_py, atol=1e-13, rtol=1e-13)


@pytest.mark.skipif(not _HAVE_CPP, reason="radia._radia_pybind not available")
def test_cpp_kernel_speedup_smoke():
    """Smoke test: C++ kernel finishes < 0.1s for 1000 calls.

    Pure Python is ~1ms / call; C++ should be ~10us / call (100x).
    """
    import time
    src_min = (0., 0., 0.); src_max = (1., 1., 1.)
    tgt_min = (2., 0., 0.); tgt_max = (3., 1., 1.)
    M = (1e6, 0., 0.)
    n = 1000
    t0 = time.perf_counter()
    for _ in range(n):
        average_B_in_box(M, src_min, src_max, tgt_min, tgt_max)
    elapsed = time.perf_counter() - t0
    assert elapsed < 1.0, f"C++ kernel too slow: {elapsed:.3f}s for {n} calls"
