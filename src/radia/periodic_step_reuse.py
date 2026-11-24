"""Periodic pitch-matched time-stepping with matrix, preconditioner, and state reuse.

For rotating machines or periodic linear guides where each time step advances by
an exact mesh symmetry pitch, the rotated/shifted mesh coincides with the original
mesh up to a DOF permutation (with optional +/-1 edge-orientation sign).
This allows:
  1. Assembling the system matrix A = K + (1/dt) * M once.
  2. Building one Radia SparseSolv preconditioner for A and reusing it with
     native PCG across all time steps: ``"amg"`` (Compact AMG, scalar SPD
     operators) or ``"ams"`` (Compact AMS, lowest-order HCurl operators; the
     caller supplies the discrete gradient and vertex coordinates). The choice
     is explicit; there is no default, no direct factorization and no fallback.
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

from ._residual_gate import RELATIVE_LIMIT, residual_scale
from .ngsolve_sparse import real_sparse_matrix_type

# Relative asymmetry allowed in A = K + M/dt, as a multiple of max|A|.
_SYMMETRY_TOLERANCE = 1e-12
# PCG stops on its true relative residual below the shared limit, leaving
# margin for the independent recheck in step().
_PCG_TOLERANCE = 0.5 * RELATIVE_LIMIT
_PRECONDITIONERS = ("amg", "ams")


def _native_matrix(matrix: csc_matrix):
    data = matrix.tocoo()
    return real_sparse_matrix_type().CreateFromCOO(
        data.row.tolist(), data.col.tolist(), data.data.tolist(), *data.shape)


def _build_preconditioner(kind: str, native, gradient, coordinates):
    """Construct the single Radia SparseSolv preconditioner for the solver."""
    from .sparsesolv_ngsolve import CompactAMGPreconditioner, CompactAMSPreconditioner

    if kind == "amg":
        return CompactAMGPreconditioner(native)
    return CompactAMSPreconditioner(
        native, _native_matrix(gradient), freedofs=None,
        coord_x=coordinates[:, 0].tolist(), coord_y=coordinates[:, 1].tolist(),
        coord_z=coordinates[:, 2].tolist())


def _pcg_solver(native, preconditioner, max_iterations: int):
    """Return x = solve(b) by native PCG reusing the one preconditioner."""
    from .sparsesolv_ngsolve import NativePCG

    pcg = NativePCG(native, preconditioner, None)
    source = native.CreateColVector()
    answer = native.CreateColVector()

    def solve(b: np.ndarray) -> np.ndarray:
        source.FV().NumPy()[:] = b
        iterations, _, converged = pcg.Solve(source, answer, _PCG_TOLERANCE, max_iterations)
        if not converged:
            raise RuntimeError(
                f"periodic step PCG did not converge in {iterations} iterations "
                f"(maximum {max_iterations})")
        return answer.FV().NumPy().copy()

    return solve


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

    def __init__(self, K, M, dt: float, perm, sign=None, *, preconditioner: str,
                 gradient=None, coordinates=None, max_iterations: int = 1000) -> None:
        if not isinstance(preconditioner, str) or preconditioner not in _PRECONDITIONERS:
            raise ValueError(f"preconditioner must be one of {_PRECONDITIONERS}, got {preconditioner!r}")
        if (isinstance(max_iterations, bool) or not isinstance(max_iterations, numbers.Integral)
                or max_iterations < 1):
            raise ValueError(f"max_iterations must be a positive integer, got {max_iterations!r}")
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

        if preconditioner == "amg":
            if gradient is not None or coordinates is not None:
                raise ValueError("gradient and coordinates belong to the 'ams' preconditioner")
        else:
            if gradient is None or coordinates is None:
                raise ValueError("'ams' requires the discrete gradient and vertex coordinates")
            gradient = _real_finite_matrix(gradient, "gradient")
            vertices = gradient.shape[1]
            if gradient.shape[0] != n:
                raise ValueError(f"gradient must have {n} rows, got shape {gradient.shape}")
            coordinates = np.asarray(coordinates)
            if coordinates.dtype == np.bool_ or not (np.issubdtype(coordinates.dtype, np.integer)
                                                    or np.issubdtype(coordinates.dtype, np.floating)):
                raise TypeError(f"coordinates must be real, got dtype {coordinates.dtype}")
            if coordinates.shape != (vertices, 3):
                raise ValueError(f"coordinates must have shape ({vertices}, 3), got {coordinates.shape}")
            coordinates = coordinates.astype(np.float64)
            if not np.all(np.isfinite(coordinates)):
                raise ValueError("coordinates must be finite")

        self.n = n
        self.dt = float(dt)
        self.perm = perm_arr
        self.sign = sign_arr
        self.M = M
        self.A = csc_matrix(K + (1.0 / self.dt) * self.M)
        scale = float(np.max(np.abs(self.A.data))) if self.A.nnz else 0.0
        asymmetry = float(np.max(np.abs((self.A - self.A.T).data), initial=0.0))
        if scale == 0.0 or asymmetry > _SYMMETRY_TOLERANCE * scale:
            raise ValueError(
                f"A = K + M/dt must be symmetric positive definite for PCG; "
                f"max|A - A^T| = {asymmetry:.3e}, max|A| = {scale:.3e}")
        if np.any(self.A.diagonal() <= 0.0):
            raise ValueError("A = K + M/dt must be positive definite: nonpositive diagonal entry")
        self.preconditioner = preconditioner
        native = _native_matrix(self.A)
        self._solve = _pcg_solver(
            native, _build_preconditioner(preconditioner, native, gradient, coordinates),
            int(max_iterations))
        self.preconditioner_setups = 1
        self.steps_solved = 0
        self.last_relative_residual = math.nan

    def shift_state(self, x_prev) -> np.ndarray:
        """Map previous-step DOF state onto the current pitch-shifted mesh."""
        x = _real_finite_vector(x_prev, self.n, "x_prev")
        return self.sign * x[self.perm]

    def step(self, rhs, x_prev=None) -> np.ndarray:
        """Advance one periodic pitch step reusing the single preconditioner.

        Raises RuntimeError, without counting the step, when PCG does not
        converge (an indefinite A breaks it loudly), the solution is not finite
        or its true relative residual exceeds RELATIVE_LIMIT.
        """
        b = _real_finite_vector(rhs, self.n, "rhs")
        if x_prev is not None:
            b = b + (1.0 / self.dt) * (self.M @ self.shift_state(x_prev))
        x_next = self._solve(b)
        residual = self.A @ x_next - b
        relative = float(np.linalg.norm(residual)) / residual_scale(float(np.linalg.norm(b)))
        if not np.all(np.isfinite(x_next)) or not math.isfinite(relative) or relative > RELATIVE_LIMIT:
            raise RuntimeError(
                f"periodic step true relative residual {relative:.3e} exceeds {RELATIVE_LIMIT:g}")
        self.last_relative_residual = relative
        self.steps_solved += 1
        return x_next
