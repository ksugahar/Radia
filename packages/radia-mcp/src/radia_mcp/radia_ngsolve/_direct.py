"""Checked SparseCholesky solve for symmetric (real or complex-symmetric) FE systems.

SparseCholesky returns a wrong answer without an error for a singular or
nonsymmetric matrix, so every solve is accepted by its true residual on the free
rows: ||r||/||b|| <= 1e-6. Backward error cannot relax that contract.
The computational backend supplies the shared acceptance limit.
"""
from __future__ import annotations

import numpy as np

DIRECT_INVERSE = "sparsecholesky"


def check_residual(matrix, residual, solution, rhs, free, what, *, reference_norm=0.0):
    from ._vendor.residual_gate import check_true_residual

    return check_true_residual(matrix, residual, solution, rhs, free, what,
                               reference_norm=reference_norm)


def require_direct_inverse(inverse: str) -> None:
    """Reject any direct backend other than SparseCholesky (PARDISO is not used)."""
    if inverse != DIRECT_INVERSE:
        raise ValueError(
            f"direct solves use inverse={DIRECT_INVERSE!r}; got {inverse!r}")


def factor(matrix, freedofs, what):
    """SparseCholesky factor with an actionable message on the index-overflow ceiling.

    The factorisation of a 3D HCurl (curl-curl) system overflows around ~84k
    order-2 elements and surfaces as "bad array new length": integer overflow in
    the fill-in, not out-of-RAM. Residual checks cannot prevent it.
    """
    try:
        return matrix.Inverse(freedofs, inverse=DIRECT_INVERSE)
    except Exception as e:
        if "array new length" in str(e).lower():
            raise RuntimeError(
                f"Direct solver '{DIRECT_INVERSE}' overflowed the factorization index "
                f"space in {what} -- the ~84k order-2 3D HCurl (curl-curl) ceiling: "
                "integer overflow in the sparse fill-in, NOT out-of-RAM. Switch to an "
                "iterative solver (CG with a BDDC/AMG preconditioner registered before "
                "a.Assemble()) or coarsen the mesh / lower the element order."
            ) from e
        raise


class CheckedInverse:
    """Factor once; every ``inverse * rhs`` is accepted by its true residual."""

    def __init__(self, matrix, freedofs, what):
        self.matrix = matrix
        self.freedofs = freedofs
        self.what = what
        self.inverse = factor(matrix, freedofs, what)
        self._free = None

    def __mul__(self, rhs):
        solution = rhs.CreateVector()
        solution.data = self.inverse * rhs
        residual = rhs.CreateVector()
        residual.data = rhs - self.matrix * solution
        if self._free is None:
            self._free = np.fromiter((bool(bit) for bit in self.freedofs),
                                     dtype=bool, count=len(solution))
        check_residual(self.matrix, residual.FV().NumPy(), solution.FV().NumPy(),
                       rhs.FV().NumPy(), self._free, self.what)
        return solution


def solve_symmetric(matrix, freedofs, rhs, what):
    """Return ``A^-1 rhs`` on the free rows (zero elsewhere), residual-checked."""
    return CheckedInverse(matrix, freedofs, what) * rhs
