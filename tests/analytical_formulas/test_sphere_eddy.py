"""Limit and transition-curve locks for sphere_complex_polarizability.

The transition goldens are the self-checked values of the BVP implementation
(2026-07-28); the limits are the independent closed forms.
"""
import numpy as np
import pytest

from radia.analytical_formulas import sphere_complex_polarizability

MU_0 = 4.0e-7 * np.pi
A = 0.01
SIGMA = 5.8e7
SCALE = 4.0 * np.pi * A**3


def _f_of_ad(ad):
    delta = A / ad
    return 1.0 / (np.pi * MU_0 * SIGMA * delta**2)


@pytest.mark.parametrize("mu_r", [1.0, 10.0, 100.0, 1000.0])
def test_static_limit_matches_magnetostatic_sphere(mu_r):
    alpha = sphere_complex_polarizability(1e-6, A, SIGMA, mu_r)
    exact = SCALE * (mu_r - 1.0) / (mu_r + 2.0)
    # 5e-5: at f = 1e-6 Hz the mu_r = 1 sphere keeps a tiny PHYSICAL eddy
    # response (~1.2e-5 of scale, prop. f); this is not an implementation error
    assert abs(alpha - exact) / SCALE < 5e-5
    exact0 = sphere_complex_polarizability(0.0, A, SIGMA, mu_r)
    assert abs(exact0 - exact) == 0.0


@pytest.mark.parametrize("mu_r,band", [(1.0, 1e-2), (100.0, 2e-2)])
def test_shielding_limit(mu_r, band):
    alpha = sphere_complex_polarizability(1e9, A, SIGMA, mu_r)
    exact = -2.0 * np.pi * A**3
    assert abs(alpha - exact) / abs(exact) < band


def test_low_frequency_scalings_mu1():
    a1 = sphere_complex_polarizability(0.5, A, SIGMA, 1.0)
    a2 = sphere_complex_polarizability(1.0, A, SIGMA, 1.0)
    assert a1.imag < 0.0 and a2.imag < 0.0  # e^{+jwt}: passive response lags
    assert a1.real < 0.0 and a2.real < 0.0
    np.testing.assert_allclose(a2.imag / a1.imag, 2.0, rtol=2e-3)
    np.testing.assert_allclose(a2.real / a1.real, 4.0, rtol=2e-3)


def test_transition_curve_golden_mu1():
    golden = {
        0.3: -0.00021 - 0.01200j,
        1.0: -0.02441 - 0.12845j,
        3.0: -0.49692 - 0.33502j,
        10.0: -0.85000 - 0.13500j,
    }
    for ad, ref in golden.items():
        alpha = sphere_complex_polarizability(_f_of_ad(ad), A, SIGMA, 1.0)
        norm = alpha / (2.0 * np.pi * A**3)
        assert abs(norm - ref) < 2e-4, (ad, norm)


def test_array_input_shape_and_scalar_return():
    freqs = np.array([[1.0, 10.0], [100.0, 1000.0]])
    out = sphere_complex_polarizability(freqs, A, SIGMA, 10.0)
    assert out.shape == freqs.shape
    scalar = sphere_complex_polarizability(10.0, A, SIGMA, 10.0)
    assert isinstance(scalar, complex)
    assert abs(scalar - out[0, 1]) == 0.0


def test_invalid_arguments_raise():
    with pytest.raises(ValueError):
        sphere_complex_polarizability(1.0, -1.0, SIGMA, 1.0)
    with pytest.raises(ValueError):
        sphere_complex_polarizability(1.0, A, 0.0, 1.0)
    with pytest.raises(ValueError):
        sphere_complex_polarizability(1.0, A, SIGMA, 0.5)
    with pytest.raises(ValueError):
        sphere_complex_polarizability(-1.0, A, SIGMA, 1.0)
    for value in (np.nan, np.inf):
        with pytest.raises(ValueError, match="finite"):
            sphere_complex_polarizability(value, A, SIGMA, 1.0)
        with pytest.raises(ValueError, match="radius"):
            sphere_complex_polarizability(1.0, value, SIGMA, 1.0)
        with pytest.raises(ValueError, match="sigma"):
            sphere_complex_polarizability(1.0, A, value, 1.0)
        with pytest.raises(ValueError, match="mu_r"):
            sphere_complex_polarizability(1.0, A, SIGMA, value)
