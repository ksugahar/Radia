"""Checked SparseCholesky solve for symmetric (real or complex-symmetric) FE systems.

SparseCholesky returns a wrong answer without an error for a singular or
nonsymmetric matrix, so every solve is accepted by its true residual on the free
rows: ||r||/||b|| <= 1e-6. Backward error cannot relax that contract.
The computational backend supplies the shared acceptance limit.
"""
from __future__ import annotations

import numpy as np


def check_residual(matrix, residual, solution, rhs, free, what, *, reference_norm=0.0):
    from radia._residual_gate import check_true_residual

    return check_true_residual(matrix, residual, solution, rhs, free, what,
                               reference_norm=reference_norm)


def solve_symmetric(matrix, freedofs, rhs, what):
    """Return ``A^-1 rhs`` on the free rows (zero elsewhere), residual-checked."""
    solution = rhs.CreateVector()
    solution.data = matrix.Inverse(freedofs, inverse="sparsecholesky") * rhs
    residual = rhs.CreateVector()
    residual.data = rhs - matrix * solution
    free = np.fromiter((bool(bit) for bit in freedofs), dtype=bool, count=len(solution))
    check_residual(matrix, residual.FV().NumPy(), solution.FV().NumPy(),
                   rhs.FV().NumPy(), free, what)
    return solution
