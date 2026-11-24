"""Surrogate tests for PeriodicStepReuseSolver.

Guards: one NGSolve sparsecholesky construction reused and exact pitch-step eddy-current state
permutation against an independent dense per-step solve; signed permutations;
per-step true residual; loud rejection of non-integer, non-finite or complex
inputs without counting a failed step.
Assumes: linear SPD A = K + M/dt on a periodic DOF ring; non-SPD A raises
without any factorization attempt or fallback.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csc_matrix, diags

import radia.periodic_step_reuse as periodic_step_reuse
from radia._residual_gate import RELATIVE_LIMIT
from radia.periodic_step_reuse import PeriodicStepReuseSolver


def _periodic_ring_matrices(n: int = 12):
    # Symmetric positive-definite 1D periodic ring K and conductivity mass M
    main = np.full(n, 4.0)
    off = np.full(n - 1, -1.0)
    K = diags([off, main, off], [-1, 0, 1], shape=(n, n), format="lil")
    K[0, n - 1] = -1.0
    K[n - 1, 0] = -1.0
    M = diags([np.full(n, 0.5)], [0], shape=(n, n), format="csc")
    return csc_matrix(K), M


@pytest.mark.parametrize("signed", [False, True])
def test_periodic_step_reuse_matches_dense_reference_with_single_factorization(signed, monkeypatch):
    n = 12
    dt = 0.01
    K, M = _periodic_ring_matrices(n)
    perm = np.roll(np.arange(n), -1)  # 1-pitch rotational shift
    sign = np.where(np.arange(n) % 2 == 0, 1.0, -1.0) if signed else None
    constructions = []
    real_factor = periodic_step_reuse._sparsecholesky

    def counting_factor(A):
        constructions.append(A.shape)
        return real_factor(A)

    monkeypatch.setattr(periodic_step_reuse, "_sparsecholesky", counting_factor)
    solver = PeriodicStepReuseSolver(K, M, dt=dt, perm=perm, sign=sign)
    rng = np.random.default_rng(42)

    # Independent reference: dense matrices, explicit permutation matrix P with
    # P[i, perm[i]] = sign[i], and numpy.linalg.solve per step.
    A_dense = K.toarray() + M.toarray() / dt
    P = np.zeros((n, n))
    P[np.arange(n), perm] = 1.0 if sign is None else sign
    x_reuse = np.zeros(n)
    x_ref = np.zeros(n)
    for step in range(8):
        rhs = rng.standard_normal(n)
        x_reuse = solver.step(rhs, x_prev=x_reuse)
        x_ref = np.linalg.solve(A_dense, rhs + M.toarray() @ (P @ x_ref) / dt)
        rel_err = np.linalg.norm(x_reuse - x_ref) / np.linalg.norm(x_ref)
        assert rel_err <= 1e-12, f"step {step}: actual={rel_err}, expected<=1e-12"
        assert solver.last_relative_residual <= RELATIVE_LIMIT, (
            f"step {step}: actual={solver.last_relative_residual}, expected<={RELATIVE_LIMIT}")

    assert constructions == [(n, n)], f"actual sparsecholesky constructions={constructions}, expected one"
    assert solver.factorization_count == 1
    assert solver.steps_solved == 8


def test_periodic_step_reuse_rejects_invalid_construction():
    K, M = _periodic_ring_matrices(6)
    bad = [
        dict(dt=0.0, perm=np.arange(6)),
        dict(dt=float("nan"), perm=np.arange(6)),
        dict(dt=0.01, perm=[0, 0, 2, 3, 4, 5]),
        dict(dt=0.01, perm=np.arange(6), sign=[1, 1, 0, 1, 1, 1]),
        dict(dt=0.01, perm=np.arange(6), sign=[1, 1, np.nan, 1, 1, 1]),
    ]
    for kwargs in bad:
        with pytest.raises(ValueError):
            PeriodicStepReuseSolver(K, M, **kwargs)
    for perm in ([0.0, 1.5, 2.0, 3.0, 4.0, 5.0], [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [True] * 6):
        with pytest.raises(TypeError):
            PeriodicStepReuseSolver(K, M, dt=0.01, perm=perm)
    with pytest.raises(TypeError):
        PeriodicStepReuseSolver(K.astype(complex), M, dt=0.01, perm=np.arange(6))
    K_nan = K.copy()
    K_nan[0, 0] = np.nan
    with pytest.raises(ValueError):
        PeriodicStepReuseSolver(K_nan, M, dt=0.01, perm=np.arange(6))


def test_periodic_step_reuse_rejects_non_spd_without_fallback(monkeypatch):
    K, M = _periodic_ring_matrices(6)
    calls = []
    monkeypatch.setattr(periodic_step_reuse, "_sparsecholesky", lambda A: calls.append(A) or None)
    K_asym = K.tolil()
    K_asym[0, 1] = -2.0
    with pytest.raises(ValueError, match="symmetric positive definite"):
        PeriodicStepReuseSolver(csc_matrix(K_asym), M, dt=0.01, perm=np.arange(6))
    with pytest.raises(ValueError, match="positive definite"):
        PeriodicStepReuseSolver(-K, M * 0.0, dt=0.01, perm=np.arange(6))
    assert calls == [], f"actual factorization attempts={len(calls)}, expected=0"


def test_periodic_step_reuse_rejects_bad_step_inputs_without_counting():
    K, M = _periodic_ring_matrices(6)
    solver = PeriodicStepReuseSolver(K, M, dt=0.01, perm=np.roll(np.arange(6), 1))
    good = np.ones(6)
    for rhs, x_prev, error in (
        (np.array([np.nan, 1, 1, 1, 1, 1]), None, ValueError),
        (np.array([np.inf, 1, 1, 1, 1, 1]), None, ValueError),
        (np.array([1j, 1, 1, 1, 1, 1]), None, TypeError),
        (good, np.array([np.nan, 1, 1, 1, 1, 1]), ValueError),
        (good, np.ones(5), ValueError),
    ):
        with pytest.raises(error):
            solver.step(rhs, x_prev=x_prev)
    assert solver.steps_solved == 0, f"actual={solver.steps_solved}, expected=0"
    solver.step(good)
    assert solver.steps_solved == 1, f"actual={solver.steps_solved}, expected=1"


def test_periodic_step_reuse_raises_on_true_residual_without_counting():
    K, M = _periodic_ring_matrices(6)
    solver = PeriodicStepReuseSolver(K, M, dt=0.01, perm=np.arange(6))
    exact = solver._solve
    solver._solve = lambda b: exact(b) * (1.0 + 1e-3)
    with pytest.raises(RuntimeError, match="true relative residual"):
        solver.step(np.ones(6))
    assert solver.steps_solved == 0, f"actual={solver.steps_solved}, expected=0"
