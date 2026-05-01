"""Tests for radia.analytical_formulas.gauss_legendre (Part 3 §4)."""

import math

import numpy as np
import pytest

from radia.analytical_formulas import (
    gauss_legendre_integrate,
    gauss_legendre_integrate_2d,
    gauss_legendre_nodes_weights,
)
from radia.analytical_formulas.gauss_legendre import N_MAX


# ---------------------------------------------------------------------------
# Nodes and weights
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [1, 2, 3, 5, 10, 16, 24])
def test_nodes_weights_basic_properties(n):
    x, w = gauss_legendre_nodes_weights(n)
    assert x.shape == (n,)
    assert w.shape == (n,)
    # Nodes are in [-1, +1] sorted ascending.
    assert np.all((x >= -1.0) & (x <= 1.0))
    assert np.all(np.diff(x) > 0)
    # Weights are positive and sum to 2.
    assert np.all(w > 0)
    assert np.sum(w) == pytest.approx(2.0, rel=1e-13)


def test_n_max_ceiling_enforced():
    with pytest.raises(ValueError):
        gauss_legendre_nodes_weights(N_MAX + 1)
    with pytest.raises(ValueError):
        gauss_legendre_nodes_weights(0)


def test_symmetric_nodes_around_zero():
    """Gauss-Legendre nodes are symmetric around 0."""
    for n in (3, 5, 8, 12):
        x, _ = gauss_legendre_nodes_weights(n)
        np.testing.assert_allclose(x, -x[::-1], atol=1e-14)


# ---------------------------------------------------------------------------
# Polynomial exactness: 2 n - 1 degree
# ---------------------------------------------------------------------------


def test_exact_for_polynomials_up_to_degree_2n_minus_1():
    """An n-point rule integrates polynomials of degree 2 n - 1 exactly."""
    for n in (1, 3, 8, 12):
        deg = 2 * n - 1
        # Use x**deg whose integral on [-1, 1] is 0 if deg odd, else 2/(deg+1).
        if deg % 2 == 0:
            expected = 2.0 / (deg + 1)
        else:
            expected = 0.0
        val = gauss_legendre_integrate(lambda x, d=deg: x ** d, -1.0, 1.0, n=n)
        assert val == pytest.approx(expected, abs=1e-12)


# ---------------------------------------------------------------------------
# Application: known transcendental integrals
# ---------------------------------------------------------------------------


def test_integrate_sin_on_zero_pi():
    val = gauss_legendre_integrate(np.sin, 0.0, math.pi, n=10)
    assert val == pytest.approx(2.0, abs=1e-12)


def test_integrate_exp_on_minus_one_one():
    val = gauss_legendre_integrate(np.exp, -1.0, 1.0, n=12)
    expected = math.e - 1.0 / math.e
    assert val == pytest.approx(expected, rel=1e-12)


def test_integrate_smooth_function_general_interval():
    """Affine map for a non-symmetric interval."""
    val = gauss_legendre_integrate(np.cos, 0.0, 1.0, n=12)
    assert val == pytest.approx(math.sin(1.0), rel=1e-13)


def test_integrate_2d_polynomial_exact():
    """An n x n tensor product is exact for a polynomial that has degree
    <= 2 n - 1 in each variable separately.
    """
    val = gauss_legendre_integrate_2d(
        lambda X, Y: X ** 2 + Y ** 2, 0.0, 1.0, 0.0, 1.0, nx=4, ny=4
    )
    assert val == pytest.approx(2.0 / 3.0, abs=1e-13)


def test_integrate_2d_smooth_function():
    val = gauss_legendre_integrate_2d(
        lambda X, Y: np.sin(X) * np.cos(Y),
        0.0, math.pi, 0.0, math.pi / 2.0,
        nx=10, ny=10,
    )
    # int_0^pi sin x dx * int_0^{pi/2} cos y dy = 2 * 1 = 2.
    assert val == pytest.approx(2.0, rel=1e-13)
