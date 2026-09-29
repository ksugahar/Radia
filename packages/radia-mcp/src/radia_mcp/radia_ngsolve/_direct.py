"""Checked SparseCholesky solve for symmetric (real or complex-symmetric) FE systems.

SparseCholesky returns a wrong answer without an error for a singular or
nonsymmetric matrix, so every solve is accepted by its true residual on the free
rows: ||r||/||b|| <= 1e-8. Backward error cannot relax that contract.
Radia-mcp is a separate distribution, so it carries its own check.
"""
from __future__ import annotations

import math

import numpy as np


def check_residual(matrix, residual, solution, rhs, free, what):
    x, r, b = solution[free], residual[free], rhs[free]
    r_norm = float(np.linalg.norm(r))
    relative = r_norm / max(float(np.linalg.norm(b)), 1e-300)
    if np.all(np.isfinite(x)) and math.isfinite(relative) and relative <= 1e-8:
        return relative
    raise RuntimeError(f"{what}: nonfinite solution or true relative residual "
                       f"{relative:.3e} exceeds 1e-8")


def solve_nonsymmetric(matrix, freedofs, rhs, what):
    """UMFPACK LU for the axisymmetric A-V net-current system, residual-checked.

    The planar A-V systems are symmetrized by the 1/s-scaled potential
    W = Vc / (j omega) (Kameari) and use :func:`solve_symmetric`.  In
    ``solve_axi_eddy`` the constraint weights A by r but the A row does not, so
    the same scaling leaves it nonsymmetric; its current definition is open, and
    until it is settled this is the one LU use."""
    solution = rhs.CreateVector()
    solution.data = matrix.Inverse(freedofs, inverse="umfpack") * rhs
    residual = rhs.CreateVector()
    residual.data = rhs - matrix * solution
    free = np.fromiter((bool(bit) for bit in freedofs), dtype=bool, count=len(solution))
    check_residual(matrix, residual.FV().NumPy(), solution.FV().NumPy(),
                   rhs.FV().NumPy(), free, what)
    return solution


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
