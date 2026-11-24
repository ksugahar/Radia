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
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csc_matrix, issparse
from scipy.sparse.linalg import splu


class PeriodicStepReuseSolver:
    """Reusable linear time-stepping solver for pitch-matched periodic motion."""

    def __init__(self, K, M, dt: float, perm, sign=None) -> None:
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError(f"dt must be finite and positive, got {dt}")
        if not issparse(K) or not issparse(M):
            raise TypeError("K and M must be scipy sparse matrices")
        if K.shape != M.shape or K.shape[0] != K.shape[1]:
            raise ValueError(f"K and M must be square with matching shape, got {K.shape} and {M.shape}")

        n = K.shape[0]
        perm_arr = np.asarray(perm, dtype=np.int64)
        if perm_arr.shape != (n,) or np.any(perm_arr < 0) or np.any(perm_arr >= n) or len(np.unique(perm_arr)) != n:
            raise ValueError("perm must be a bijective permutation of 0..n-1")

        if sign is None:
            sign_arr = np.ones(n, dtype=np.float64)
        else:
            sign_arr = np.asarray(sign, dtype=np.float64)
            if sign_arr.shape != (n,) or not np.all(np.isin(sign_arr, (-1.0, 1.0))):
                raise ValueError("sign must have shape (n,) with entries in {-1, +1}")

        self.n = n
        self.dt = float(dt)
        self.perm = perm_arr
        self.sign = sign_arr
        self.M = csc_matrix(M)
        self.A = csc_matrix(K + (1.0 / self.dt) * self.M)
        self._lu = splu(self.A)
        self.factorization_count = 1
        self.steps_solved = 0

    def shift_state(self, x_prev: np.ndarray) -> np.ndarray:
        """Map previous-step DOF state onto the current pitch-shifted mesh."""
        x = np.asarray(x_prev, dtype=np.float64)
        if x.shape != (self.n,):
            raise ValueError(f"x_prev must have shape ({self.n},), got {x.shape}")
        return self.sign * x[self.perm]

    def step(self, rhs: np.ndarray, x_prev: np.ndarray | None = None) -> np.ndarray:
        """Advance one periodic pitch step reusing the single matrix factorization."""
        b = np.asarray(rhs, dtype=np.float64)
        if b.shape != (self.n,):
            raise ValueError(f"rhs must have shape ({self.n},), got {b.shape}")
        if x_prev is not None:
            b = b + (1.0 / self.dt) * (self.M @ self.shift_state(x_prev))
        x_next = self._lu.solve(b)
        self.steps_solved += 1
        return x_next
