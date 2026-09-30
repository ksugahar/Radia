"""Strict true-residual acceptance of FE solves; backward error is diagnostic only."""
from __future__ import annotations

import math

import numpy as np

RELATIVE_LIMIT = 1e-6


def frobenius_norm_on(matrix, free) -> float:
    """||A||_F restricted to the free rows and columns (``CSR()`` holds both triangles)."""
    values, columns, pointers = (np.asarray(item) for item in matrix.CSR())
    pointers = pointers.astype(np.int64)
    rows = np.repeat(np.arange(len(pointers) - 1), np.diff(pointers))
    keep = free[rows] & free[columns.astype(np.int64)]
    return float(np.sqrt(np.sum(np.abs(values[keep]) ** 2)))


def residual_scale(rhs_norm: float, reference_norm: float = 0.0) -> float:
    """Use an explicit fixed physical load scale, never the computed solution.

    A caller may supply its original effective-load norm for Newton correction
    equations. The reference must be fixed before solving, finite and nonnegative.
    """
    if not math.isfinite(reference_norm) or reference_norm < 0:
        raise ValueError("reference_norm must be finite and nonnegative")
    return max(float(rhs_norm), float(reference_norm), 1e-300)


def check_true_residual(matrix, residual, solution, rhs, free, what,
                        *, reference_norm: float = 0.0) -> float:
    """Require the shared limit; return ||r||/max(||b||, reference_norm).

    ``residual``, ``solution`` and ``rhs`` are full-length arrays; ``free`` is
    the boolean mask of the rows that were solved.
    """
    x = solution[free]
    r = residual[free]
    b = rhs[free]
    r_norm = float(np.linalg.norm(r))
    relative = r_norm / residual_scale(float(np.linalg.norm(b)), reference_norm)
    backward = math.inf
    finite = bool(np.all(np.isfinite(x))) and math.isfinite(relative)
    if finite and relative > RELATIVE_LIMIT:
        backward = r_norm / max(frobenius_norm_on(matrix, free) * float(np.linalg.norm(x))
                                + float(np.linalg.norm(b)), 1e-300)
    if not finite or relative > RELATIVE_LIMIT:
        raise RuntimeError(
            f"{what} true relative residual {relative:.3e} exceeds {RELATIVE_LIMIT:g}; diagnostic "
            f"normwise backward error is {backward:.3e}")
    return relative
