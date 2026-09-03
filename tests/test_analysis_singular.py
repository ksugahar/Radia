"""No-Fallback tests for radia.analysis MultiPortResult S / Y parameters.

A singular Z at one frequency must give NaN at that frequency (clearly invalid)
plus a warning -- NOT a silent pseudo-inverse (get_Y) or a plausible zero
(get_S) that the caller would trust.
"""
import warnings

import numpy as np

from radia.analysis import AnalysisType, MultiPortResult, SolverType


def _result(freqs, Z):
    return MultiPortResult(
        analysis_type=list(AnalysisType)[0], solver_type=list(SolverType)[0],
        success=True, frequencies=np.asarray(freqs, dtype=float), n_ports=2,
        Z_matrix=np.asarray(Z))


def test_get_Y_singular_is_nan_not_pinv():
    """A singular Z -> Y = NaN at that frequency (no defined admittance) plus a
    warning, NOT a silent pseudo-inverse (a different, un-auditable min-norm)."""
    Z = np.zeros((2, 2, 2), dtype=complex)
    Z[0] = 50.0 * np.eye(2)                        # invertible
    Z[1] = np.ones((2, 2))                         # rank-1 -> singular
    res = _result([1.0e6, 2.0e6], Z)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        Y = res.get_Y()
    assert np.isfinite(Y[0]).all()
    assert np.isnan(Y[1]).all()                    # singular -> NaN, not min-norm
    assert any("Singular" in str(x.message) for x in w)


def test_get_S_singular_is_nan_not_zero():
    """A singular (Z + Z0 I) -> S = NaN at that frequency, NOT a plausible zero
    S-matrix (which reads as a perfectly-matched / fully-absorbing port)."""
    Z0 = 50.0
    Z = np.zeros((2, 2, 2), dtype=complex)
    Z[0] = 50.0 * np.eye(2)                         # Z + Z0 I = 100 I (invertible)
    Z[1] = np.ones((2, 2)) - Z0 * np.eye(2)         # Z + Z0 I = ones -> singular
    res = _result([1.0e6, 2.0e6], Z)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        S = res.get_S(Z0=Z0)
    assert np.isfinite(S[0]).all()
    assert np.isnan(S[1]).all()                     # singular -> NaN, not zero
    assert any("Singular" in str(x.message) for x in w)


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v"])
