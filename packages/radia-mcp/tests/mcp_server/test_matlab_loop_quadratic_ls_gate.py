"""Public-safe MATLAB loop optimization gate.

The MATLAB Gypsilab lane keeps a readable quadratic least-squares example:
objective value, analytic gradient, central finite-difference gradient, and
initial point are all explicit.  This test mirrors that solver-independent
contract so the public MCP side can check the same teaching identity without
depending on MATLAB or private paths.
"""

import numpy as np


def _objective(A, b, x):
    r = A @ x - b
    return 0.5 * float(r @ r)


def _gradient(A, b, x):
    return A.T @ (A @ x - b)


def _finite_difference_gradient(A, b, x, h=1.0e-6):
    out = np.zeros_like(x, dtype=float)
    for i in range(x.size):
        xp = x.copy()
        xm = x.copy()
        xp[i] += h
        xm[i] -= h
        out[i] = (_objective(A, b, xp) - _objective(A, b, xm)) / (2.0 * h)
    return out


def test_quadratic_ls_gradient_contract_matches_finite_difference():
    A = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, -1.0],
        ]
    )
    x_true = np.array([2.0, -1.0])
    b = A @ x_true
    x0 = np.array([0.25, -0.5])

    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    analytic = _gradient(A, b, x0)
    finite_difference = _finite_difference_gradient(A, b, x0)

    assert np.allclose(x, x_true, atol=1.0e-12)
    assert np.linalg.norm(_gradient(A, b, x)) < 1.0e-12
    assert _objective(A, b, x0) > _objective(A, b, x)
    assert np.max(np.abs(analytic - finite_difference)) < 1.0e-8
