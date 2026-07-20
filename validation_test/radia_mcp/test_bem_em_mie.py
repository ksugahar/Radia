# -*- coding: utf-8 -*-
"""Analytic gate for the electromagnetic (vector Maxwell) Mie PEC-sphere RCS.

These are CLOSED-FORM checks (no BEM solve): the EM Mie cross-sections must hit
their two textbook limits.  This is the analytic gate a future MaxwellSingleLayer
EFIE/CFIE PEC-sphere RCS solve in ngsolve.bem will be validated against.
"""
import pytest

pytest.importorskip("scipy")
pytest.importorskip("ngsolve")

from radia_mcp.radia_ngsolve.bem_integral import (
    mie_pec_coefficients,
    mie_pec_sphere_rcs,
    maxwell_efie_pec_sphere_rcs,
)


def test_em_mie_rayleigh_limit():
    """ka -> 0 : sigma_back / (pi a^2) -> 9 (ka)^4 (electric + magnetic dipole)."""
    for ka in (0.05, 0.1):
        r = mie_pec_sphere_rcs(ka)
        ratio = r["sigma_back_over_pa2"] / (9.0 * ka ** 4)
        assert abs(ratio - 1.0) < 0.01, (ka, ratio)


def test_em_mie_geometric_limit():
    """ka -> inf : sigma_back/(pi a^2) -> 1 (specular), sigma_sca/(pi a^2) -> 2 (extinction paradox)."""
    r = mie_pec_sphere_rcs(100.0)
    assert abs(r["sigma_back_over_pa2"] - 1.0) < 0.01
    assert abs(r["sigma_sca_over_pa2"] - 2.0) < 0.02


def test_em_mie_coefficients_shapes_and_magnitude():
    """a_n,b_n are bounded by 1 (lossless PEC partial-wave reflection) and decay past n~ka."""
    import numpy as np
    n, a, b = mie_pec_coefficients(5.0)
    assert n[0] == 1 and len(n) == len(a) == len(b)
    assert np.all(np.abs(a) <= 1.0 + 1e-9) and np.all(np.abs(b) <= 1.0 + 1e-9)
    assert abs(a[-1]) < 1e-3 and abs(b[-1]) < 1e-3   # high-order terms negligible


@pytest.mark.parametrize("ka", [0.5, 1.0, 2.0])
def test_em_mie_efie_solve_matches_rcs(ka):
    """ACTUAL Maxwell EFIE BEM solve (MaxwellSingleLayerPotentialOperator) reproduces the analytic EM
    Mie monostatic RCS across Rayleigh -> resonance, with ONE global normalization (not a per-ka fit)."""
    r = maxwell_efie_pec_sphere_rcs(ka, maxh=0.35, order=2)
    assert r["rel_err"] < 5e-3, (ka, r["sigma_back_over_pa2"], r["mie"], r["rel_err"])
