"""Tests for radia.analytical_formulas.elliptic_integrals (Part 3 §3)."""

import math

import numpy as np
import pytest

from radia.analytical_formulas import (
    E_hastings_2,
    E_hastings_4,
    K_hastings_2,
    K_hastings_4,
)


def _scipy_K(k):
    from scipy.special import ellipk
    return ellipk(k * k)   # scipy uses parameter m = k**2


def _scipy_E(k):
    from scipy.special import ellipe
    return ellipe(k * k)


# ---------------------------------------------------------------------------
# Limit values
# ---------------------------------------------------------------------------


def test_K_at_k_zero_equals_pi_over_2():
    assert K_hastings_4(0.0) == pytest.approx(math.pi / 2.0, abs=1e-10)
    assert K_hastings_2(0.0) == pytest.approx(math.pi / 2.0, abs=1e-7)


def test_E_at_k_zero_equals_pi_over_2():
    assert E_hastings_4(0.0) == pytest.approx(math.pi / 2.0, abs=1e-10)
    assert E_hastings_2(0.0) == pytest.approx(math.pi / 2.0, abs=1e-7)


# ---------------------------------------------------------------------------
# Accuracy tables vs scipy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k", [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95])
def test_K_hastings_4_accuracy(k):
    """Degree-4 Hastings: PDF claims ``< 2e-8`` absolute error."""
    err = abs(K_hastings_4(k) - _scipy_K(k))
    assert err < 1e-7


@pytest.mark.parametrize("k", [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95])
def test_E_hastings_4_accuracy(k):
    err = abs(E_hastings_4(k) - _scipy_E(k))
    assert err < 1e-7


@pytest.mark.parametrize("k", [0.0, 0.3, 0.5, 0.7, 0.9])
def test_K_hastings_2_accuracy(k):
    """Degree-2 Hastings: PDF claims ``< 4e-5`` absolute error."""
    err = abs(K_hastings_2(k) - _scipy_K(k))
    assert err < 5e-5


@pytest.mark.parametrize("k", [0.0, 0.3, 0.5, 0.7, 0.9])
def test_E_hastings_2_accuracy(k):
    err = abs(E_hastings_2(k) - _scipy_E(k))
    assert err < 5e-5


# ---------------------------------------------------------------------------
# Vectorisation
# ---------------------------------------------------------------------------


def test_K_vectorised():
    ks = np.array([0.1, 0.5, 0.9])
    K = K_hastings_4(ks)
    for i, k in enumerate(ks):
        assert K[i] == pytest.approx(K_hastings_4(float(k)), rel=1e-13)


def test_E_vectorised():
    ks = np.array([0.1, 0.5, 0.9])
    E = E_hastings_4(ks)
    for i, k in enumerate(ks):
        assert E[i] == pytest.approx(E_hastings_4(float(k)), rel=1e-13)
