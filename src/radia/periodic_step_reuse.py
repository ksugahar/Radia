"""Periodic pitch-matched time-stepping with matrix, factorization, and state reuse.

For rotating machines or periodic linear guides where each time step advances by
an exact mesh symmetry pitch, the rotated/shifted mesh coincides with the original
mesh up to a DOF permutation (with optional +/-1 edge-orientation sign).
This allows:
  1. Assembling the system matrix A = K + (1/dt) * M once.
  2. Factoring (or building the AMS/preconditioner hierarchy for) A once and
     reusing it across all time steps.
  3. Carrying eddy-current / flux state from the previous step by permuting
     the previous solution vector along the pitch mapping without re-meshing.

Inputs are real and finite; anything else raises instead of being cast. Every
step checks its true relative residual against the shared limit and counts
only steps that passed.
"""

from __future__ import annotations

import math
import numbers

import numpy as np
from scipy.sparse import csc_matrix, issparse
from scipy.sparse.linalg import splu

from ._residual_gate import RELATIVE_LIMIT, residual_scale


def _real_finite_vector(value, n: int, name: str) -> np.ndarray:
    array = np.asarray(value)
    if array.dtype == np.bool_ or not (np.issubdtype(array.dtype, np.integer)
                                       or np.issubdtype(array.dtype, np.floating)):
        raise TypeError(f"{name} must be a real numeric array, got dtype {array.dtype}")
    if array.shape != (n,):
        raise ValueError(f"{name} must have shape ({n},), got {array.shape}")
    array = array.astype(np.float64)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _real_finite_matrix(matrix, name: str) -> csc_matrix:
    if not issparse(matrix):
        raise TypeError(f"{name} must be a scipy sparse matrix")
    if not (np.issubdtype(matrix.dtype, np.integer) or np.issubdtype(matrix.dtype, np.floating)):
        raise TypeError(f"{name} must be real, got dtype {matrix.dtype}")
    result = csc_matrix(matrix, dtype=np.float64)
    if not np.all(np.isfinite(result.data)):
        raise ValueError(f"{name} must be finite")
    return result


def _permutation(perm, n: int) -> np.ndarray:
    array = np.asarray(perm)
    if array.dtype == np.bool_ or not np.issubdtype(array.dtype, np.integer):
        raise TypeError(f"perm must be an integer array, got dtype {array.dtype}")
    if array.shape != (n,):
        raise ValueError(f"perm must have shape ({n},), got {array.shape}")
    array = array.astype(np.int64)
    if np.any(array < 0) or np.any(array >= n) or len(np.unique(array)) != n:
        raise ValueError("perm must be a bijective permutation of 0..n-1")
    return array


class PeriodicStepReuseSolver:
    """Reusable linear time-stepping solver for pitch-matched periodic motion."""

    def __init__(self, K, M, dt: float, perm, sign=None) -> None:
        if isinstance(dt, bool) or not isinstance(dt, numbers.Real) or not math.isfinite(dt) or dt <= 0.0:
            raise ValueError(f"dt must be a finite positive real number, got {dt!r}")
        K = _real_finite_matrix(K, "K")
        M = _real_finite_matrix(M, "M")
        if K.shape != M.shape or K.shape[0] != K.shape[1]:
            raise ValueError(f"K and M must be square with matching shape, got {K.shape} and {M.shape}")

        n = K.shape[0]
        perm_arr = _permutation(perm, n)
        if sign is None:
            sign_arr = np.ones(n, dtype=np.float64)
        else:
            sign_arr = _real_finite_vector(sign, n, "sign")
            if not np.all(np.isin(sign_arr, (-1.0, 1.0))):
                raise ValueError("sign entries must be -1 or +1")

        self.n = n
        self.dt = float(dt)
        self.perm = perm_arr
        self.sign = sign_arr
        self.M = M
        self.A = csc_matrix(K + (1.0 / self.dt) * self.M)
        self._lu = splu(self.A)
        self.factorization_count = 1
        self.steps_solved = 0
        self.last_relative_residual = math.nan

    def shift_state(self, x_prev) -> np.ndarray:
        """Map previous-step DOF state onto the current pitch-shifted mesh."""
        x = _real_finite_vector(x_prev, self.n, "x_prev")
        return self.sign * x[self.perm]

    def step(self, rhs, x_prev=None) -> np.ndarray:
        """Advance one periodic pitch step reusing the single matrix factorization.

        Raises RuntimeError, without counting the step, when the solution is not
        finite or its true relative residual exceeds RELATIVE_LIMIT.
        """
        b = _real_finite_vector(rhs, self.n, "rhs")
        if x_prev is not None:
            b = b + (1.0 / self.dt) * (self.M @ self.shift_state(x_prev))
        x_next = self._lu.solve(b)
        residual = self.A @ x_next - b
        relative = float(np.linalg.norm(residual)) / residual_scale(float(np.linalg.norm(b)))
        if not np.all(np.isfinite(x_next)) or not math.isfinite(relative) or relative > RELATIVE_LIMIT:
            raise RuntimeError(
                f"periodic step true relative residual {relative:.3e} exceeds {RELATIVE_LIMIT:g}")
        self.last_relative_residual = relative
        self.steps_solved += 1
        return x_next
