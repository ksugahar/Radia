"""Surrogate tests for PeriodicStepReuseSolver.

Guards: single-factorization reuse and exact pitch-step eddy-current state
permutation against naive per-step re-assembly and re-factorization.
Assumes: linear SPD/nonsingular sparse K and M on a periodic DOF ring.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import spsolve

_MOD_PATH = Path(__file__).resolve().parents[1] / "src" / "radia" / "periodic_step_reuse.py"
_SPEC = importlib.util.spec_from_file_location("radia_periodic_step_reuse", _MOD_PATH)
_MOD = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MOD)
PeriodicStepReuseSolver = _MOD.PeriodicStepReuseSolver


def _periodic_ring_matrices(n: int = 12):
    # Symmetric positive-definite 1D periodic ring K and conductivity mass M
    main = np.full(n, 4.0)
    off = np.full(n - 1, -1.0)
    K = diags([off, main, off], [-1, 0, 1], shape=(n, n), format="lil")
    K[0, n - 1] = -1.0
    K[n - 1, 0] = -1.0
    M = diags([np.full(n, 0.5)], [0], shape=(n, n), format="csc")
    return csc_matrix(K), M


def test_periodic_step_reuse_matches_naive_steps_with_single_factorization():
    n = 12
    dt = 0.01
    K, M = _periodic_ring_matrices(n)
    perm = np.roll(np.arange(n), -1)  # 1-pitch rotational shift
    sign = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)  # signed edge orientations

    solver = PeriodicStepReuseSolver(K, M, dt=dt, perm=perm, sign=sign)
    rng = np.random.default_rng(42)

    x_reuse = np.zeros(n)
    x_naive = np.zeros(n)
    for _ in range(8):
        rhs = rng.standard_normal(n)
        x_reuse = solver.step(rhs, x_prev=x_reuse)
        # Naive per-step re-assembly and solve
        A_naive = csc_matrix(K + (1.0 / dt) * M)
        x_naive = spsolve(A_naive, rhs + (1.0 / dt) * (M @ (sign * x_naive[perm])))
        rel_err = np.linalg.norm(x_reuse - x_naive) / max(np.linalg.norm(x_naive), 1e-15)
        assert rel_err <= 1e-12, f"actual={rel_err}, expected<=1e-12"

    assert solver.factorization_count == 1
    assert solver.steps_solved == 8


def test_periodic_step_reuse_rejects_invalid_perm_and_dt():
    K, M = _periodic_ring_matrices(6)
    with pytest.raises(ValueError):
        PeriodicStepReuseSolver(K, M, dt=0.0, perm=np.arange(6))
    with pytest.raises(ValueError):
        PeriodicStepReuseSolver(K, M, dt=0.01, perm=[0, 0, 2, 3, 4, 5])
    with pytest.raises(ValueError):
        PeriodicStepReuseSolver(K, M, dt=0.01, perm=np.arange(6), sign=[1, 1, 0, 1, 1, 1])
