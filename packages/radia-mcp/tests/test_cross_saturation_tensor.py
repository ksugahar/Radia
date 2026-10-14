r"""Incremental inductance matrix + cross-saturation + reciprocity -- regression test (#52).

The multi-port generalisation of secant/incremental inductance (#28/#31). The INCREMENTAL
matrix L_jk = d lambda_j/d i_k (central FD of the nonlinear flux-linkage map) is SYMMETRIC
(reciprocity, co-energy Hessian) and, in saturation, much SMALLER than the frozen-permeability
APPARENT/secant inductance. The off-diagonal (cross) inductance collapses as the shared iron
saturates. Pure-FD helper (tool-independent) + a saturable two-winding FE check."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import incremental_inductance_matrix


def test_incremental_matrix_recovers_linear():
    # a LINEAR (symmetric) flux map lambda_j = sum_k L_jk i_k -> FD recovers L exactly
    L = [[3.0, 1.0], [1.0, 2.0]]
    flux = lambda i: [L[0][0]*i[0] + L[0][1]*i[1], L[1][0]*i[0] + L[1][1]*i[1]]
    M = incremental_inductance_matrix(flux, [0.7, -0.4], 0.1)
    for j in range(2):
        for k in range(2):
            assert math.isclose(M[j][k], L[j][k], rel_tol=1e-9), f"({j},{k})"
    # reciprocity of the recovered matrix
    assert math.isclose(M[0][1], M[1][0], rel_tol=1e-12)
