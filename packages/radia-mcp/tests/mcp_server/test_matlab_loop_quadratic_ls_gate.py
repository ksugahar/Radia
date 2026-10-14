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


def _trace_objective(T, M, g, alpha, u):
    r = T @ u - g
    return 0.5 * float(r @ M @ r) + 0.5 * float(alpha) * float(u @ u)


def _trace_gradient(T, M, g, alpha, u):
    return T.T @ M @ (T @ u - g) + float(alpha) * u


def _trace_finite_difference_gradient(T, M, g, alpha, u, h=1.0e-6):
    out = np.zeros_like(u, dtype=float)
    for i in range(u.size):
        up = u.copy()
        um = u.copy()
        up[i] += h
        um[i] -= h
        out[i] = (
            _trace_objective(T, M, g, alpha, up)
            - _trace_objective(T, M, g, alpha, um)
        ) / (2.0 * h)
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
    assert np.max(np.abs(analytic - finite_difference)) < 1.0e-6


def test_fem_bem_trace_least_squares_contract_keeps_interior_unforced():
    # Four boundary trace rows and one interior H1 unknown. The last column is
    # invisible to the trace and should be driven to zero by Tikhonov
    # regularization rather than by a hidden constraint.
    T = np.hstack([np.eye(4), np.zeros((4, 1))])
    M = np.array(
        [
            [2.0, 0.2, 0.1, 0.0],
            [0.2, 1.5, 0.0, 0.1],
            [0.1, 0.0, 1.2, 0.2],
            [0.0, 0.1, 0.2, 1.8],
        ]
    )
    g = np.array([10.0, 20.0, 30.0, 40.0])
    alpha = 1.0e-3
    normal = T.T @ M @ T + alpha * np.eye(5)
    rhs = T.T @ M @ g
    u = np.linalg.solve(normal, rhs)
    u0 = np.ones(5)

    analytic = _trace_gradient(T, M, g, alpha, u0)
    finite_difference = _trace_finite_difference_gradient(T, M, g, alpha, u0)

    assert abs(u[-1]) < 1.0e-12
    assert np.linalg.norm(_trace_gradient(T, M, g, alpha, u)) < 1.0e-10
    assert _trace_objective(T, M, g, alpha, u0) > _trace_objective(T, M, g, alpha, u)
    assert np.max(np.abs(analytic - finite_difference)) < 1.0e-6


def test_fem_bem_tikhonov_path_exposes_regularization_tradeoff():
    T = np.hstack([np.eye(4), np.zeros((4, 1))])
    M = np.array(
        [
            [2.0, 0.2, 0.1, 0.0],
            [0.2, 1.5, 0.0, 0.1],
            [0.1, 0.0, 1.2, 0.2],
            [0.0, 0.1, 0.2, 1.8],
        ]
    )
    g = np.array([10.0, 20.0, 30.0, 40.0])
    alphas = np.array([0.0, 1.0e-3, 1.0e-1, 1.0])
    solution_norms = []
    trace_residuals = []

    for alpha in alphas:
        normal = T.T @ M @ T + alpha * np.eye(5)
        rhs = T.T @ M @ g
        if alpha == 0.0:
            u = np.linalg.lstsq(normal, rhs, rcond=None)[0]
        else:
            u = np.linalg.solve(normal, rhs)
        solution_norms.append(float(np.linalg.norm(u)))
        trace_residuals.append(float(np.linalg.norm(T @ u - g)))
        assert np.linalg.norm(_trace_gradient(T, M, g, alpha, u)) < 1.0e-10

    assert all(a >= b - 1.0e-10 for a, b in zip(solution_norms, solution_norms[1:]))
    assert all(a <= b + 1.0e-10 for a, b in zip(trace_residuals, trace_residuals[1:]))
    assert solution_norms[-1] < solution_norms[0]
    assert trace_residuals[-1] > trace_residuals[0]
