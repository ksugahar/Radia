"""Checked SparseCholesky solve for symmetric (real or complex-symmetric) FE systems.

SparseCholesky returns a wrong answer without an error for a singular or
nonsymmetric matrix, so every solve is accepted by its true residual on the free
rows: ||r||/||b|| <= 1e-8, or ||r||/||b|| <= 1e-6 with normwise backward error
||r|| / (||A||_F ||x|| + ||b||) <= 1e-12 (the rule of ``radia._residual_gate``;
radia-mcp is a separate distribution, so it carries its own copy).
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
    if not relative <= 1e-6:
        raise RuntimeError(f"{what}: true relative residual {relative:.3e} exceeds 1e-8")
    values, columns, pointers = (np.asarray(item) for item in matrix.CSR())
    rows = np.repeat(np.arange(len(pointers) - 1), np.diff(pointers.astype(np.int64)))
    keep = free[rows] & free[columns.astype(np.int64)]
    scale = float(np.sqrt(np.sum(np.abs(values[keep]) ** 2))) * float(np.linalg.norm(x))
    backward = r_norm / max(scale + float(np.linalg.norm(b)), 1e-300)
    if not (np.all(np.isfinite(x)) and math.isfinite(backward) and backward <= 1e-12):
        raise RuntimeError(f"{what}: true relative residual {relative:.3e} exceeds 1e-8 and "
                           f"backward error {backward:.3e} exceeds 1e-12")
    return relative


def solve_nonsymmetric(matrix, freedofs, rhs, what):
    """UMFPACK LU for the A-V net-current systems, residual-checked.

    Their constraint row carries -j omega sigma where the column carries -sigma,
    so the matrix is not symmetric and SparseCholesky (one triangle) would return
    a wrong answer silently.  A symmetric scaling of the constraint row is the
    route to SparseCholesky; until it is adopted this is the one LU use."""
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
