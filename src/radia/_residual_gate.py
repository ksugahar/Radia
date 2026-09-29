"""True-residual acceptance of FE linear solves (no solver fallback).

A solve passes when the relative residual ||r|| / ||b|| on the free rows is at
most 1e-8, or when the normwise backward error ||r|| / (||A||_F ||x|| + ||b||)
is at most 1e-12.  The second clause keeps exact solves of ill-conditioned but
legitimate systems: a gauged HCurl curl-curl with a load that is not
discretely divergence free reaches ||r||/||b|| ~ 1e-8 at a backward error of
1e-17.  A factorization that ignored a nonsymmetric triangle or met a singular
pivot fails both clauses (backward error ~1e-1, or non-finite values).
"""
from __future__ import annotations

import math

import numpy as np

RELATIVE_LIMIT = 1e-8
BACKWARD_LIMIT = 1e-12


def frobenius_norm_on(matrix, free) -> float:
    """||A||_F restricted to the free rows and columns (``CSR()`` holds both triangles)."""
    values, columns, pointers = (np.asarray(item) for item in matrix.CSR())
    pointers = pointers.astype(np.int64)
    rows = np.repeat(np.arange(len(pointers) - 1), np.diff(pointers))
    keep = free[rows] & free[columns.astype(np.int64)]
    return float(np.sqrt(np.sum(np.abs(values[keep]) ** 2)))


def check_true_residual(matrix, residual, solution, rhs, free, what) -> float:
    """Raise unless the solve meets one of the two clauses; return ||r||/||b||.

    ``residual``, ``solution`` and ``rhs`` are full-length arrays; ``free`` is
    the boolean mask of the rows that were solved.
    """
    x = solution[free]
    r = residual[free]
    b = rhs[free]
    r_norm = float(np.linalg.norm(r))
    relative = r_norm / max(float(np.linalg.norm(b)), 1e-300)
    backward = math.inf
    finite = bool(np.all(np.isfinite(x))) and math.isfinite(relative)
    if finite and relative > RELATIVE_LIMIT:
        backward = r_norm / max(frobenius_norm_on(matrix, free) * float(np.linalg.norm(x))
                                + float(np.linalg.norm(b)), 1e-300)
    if not finite or (relative > RELATIVE_LIMIT and backward > BACKWARD_LIMIT):
        raise RuntimeError(
            f"{what} true relative residual {relative:.3e} exceeds {RELATIVE_LIMIT:g} and its "
            f"normwise backward error {backward:.3e} exceeds {BACKWARD_LIMIT:g}")
    return relative
