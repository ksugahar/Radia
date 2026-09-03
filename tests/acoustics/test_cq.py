"""Fast CQ formula and input-contract regressions.

The NGSolve BEM comparison lives in
``validation_test/acoustics/test_cq_bem.py``.
"""
import numpy as np

from radia.acoustics import cq
from radia.acoustics import soft_sphere_scattering


def test_complex_k_series_matches_real_soft_sphere():
    # at a REAL wavenumber the complex-argument series equals the real soft-sphere series
    obs = np.array([[0.0, 0.0, -3.0], [3.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
    ref = soft_sphere_scattering(2.0, 1.0, obs)["scattered"]
    chk = cq.soft_sphere_scattering_complex_k(2.0, 1.0, obs)
    np.testing.assert_allclose(chk, ref, atol=1e-12)


def test_bdf_delta_values():
    z = 0.3 + 0.1j
    np.testing.assert_allclose(cq.bdf_delta(z, "BDF1"), 1 - z)
    np.testing.assert_allclose(cq.bdf_delta(z, "BDF2"), 1.5 - 2 * z + 0.5 * z**2)


def test_complex_k_series_finite_at_complex_k():
    # the CQ Laplace nodes are complex; the reference series must stay finite there
    obs = np.array([[0.0, 0.0, -3.0], [3.0, 0.0, 0.0]])
    v = cq.soft_sphere_scattering_complex_k(2.0 + 1.0j, 1.0, obs)
    assert np.all(np.isfinite(v)) and np.max(np.abs(v)) > 0
