"""Linear conjugate gradient -- matrix-free SPD solve for optimization inner loops."""
import os
import sys

import pytest

np = pytest.importorskip("numpy")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.topology_optimization.krylov import linear_conjugate_gradient


def _spd_system(n=8, seed=4):
    rng = np.random.default_rng(seed)
    M = rng.normal(size=(n, n))
    A = M.T @ M + n * np.eye(n)
    b = rng.normal(size=n)
    return A, b


def test_linear_cg_matches_dense_solve_and_finishes_within_dimension():
    A, b = _spd_system()
    out = linear_conjugate_gradient(A, b, tol=1e-13)
    ref = np.linalg.solve(A, b)

    assert out["converged"]
    assert out["n_iter"] <= b.size
    assert out["relative_residual"] < 1e-13
    assert np.linalg.norm(out["x"] - ref) < 1e-13
    assert out["residual_history"][0] == pytest.approx(np.linalg.norm(b))
    assert out["residual_history"][-1] == pytest.approx(out["residual_norm"])


def test_linear_cg_accepts_matrix_free_operator_and_initial_guess():
    A, b = _spd_system(n=6, seed=9)
    x0 = np.ones_like(b)
    out = linear_conjugate_gradient(lambda v: A @ v, b, x0=x0, tol=1e-12)

    assert out["converged"]
    assert np.linalg.norm(A @ out["x"] - b) / max(np.linalg.norm(b), 1.0) < 1e-12
    assert out["n_iter"] <= b.size


def test_linear_cg_preconditioner_and_zero_residual_case():
    A = np.diag([2.0, 4.0, 8.0])
    b = np.array([1.0, -2.0, 3.0])
    out = linear_conjugate_gradient(A, b, preconditioner=lambda r: r / np.diag(A), tol=1e-14)
    assert out["converged"]
    assert out["n_iter"] == 1
    assert out["x"] == pytest.approx(np.linalg.solve(A, b))

    exact = linear_conjugate_gradient(A, b, x0=np.linalg.solve(A, b))
    assert exact["converged"]
    assert exact["n_iter"] == 0


def test_linear_cg_validation_errors():
    A, b = _spd_system(n=3)
    with pytest.raises(ValueError, match="square"):
        linear_conjugate_gradient(np.ones((2, 3)), b)
    with pytest.raises(ValueError, match="vector"):
        linear_conjugate_gradient(A, np.ones((3, 1)))
    with pytest.raises(ValueError, match="tol"):
        linear_conjugate_gradient(A, b, tol=0.0)
    with pytest.raises(ValueError, match="same shape"):
        linear_conjugate_gradient(A, b, x0=np.ones(2))
    with pytest.raises(ValueError, match="symmetric positive definite"):
        linear_conjugate_gradient(np.diag([1.0, -1.0, 2.0]), b)
    with pytest.raises(ValueError, match="preconditioner"):
        linear_conjugate_gradient(A, b, preconditioner=lambda r: -r)
