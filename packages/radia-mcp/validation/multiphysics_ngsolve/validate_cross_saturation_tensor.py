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

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import incremental_inductance_matrix


def validate_incremental_matrix_recovers_linear():
    # a LINEAR (symmetric) flux map lambda_j = sum_k L_jk i_k -> FD recovers L exactly
    L = [[3.0, 1.0], [1.0, 2.0]]
    flux = lambda i: [L[0][0]*i[0] + L[0][1]*i[1], L[1][0]*i[0] + L[1][1]*i[1]]
    M = incremental_inductance_matrix(flux, [0.7, -0.4], 0.1)
    for j in range(2):
        for k in range(2):
            assert math.isclose(M[j][k], L[j][k], rel_tol=1e-9), f"({j},{k})"
    # reciprocity of the recovered matrix
    assert math.isclose(M[0][1], M[1][0], rel_tol=1e-12)


def validate_cross_saturation_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction, sqrt, InnerProduct, TaskManager
    from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue
    from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                               solve_planar_magnetostatic_nonlinear,
                                               coil_flux_linkage_2d, frozen_reluctivity, NU0)
    MUR0, BK, NEXP = 1000.0, 1.4, 6.0
    OUT, WIN, BOX, rw = 0.06, 0.03, 0.12, 0.003
    N, DEPTH, xin, xout = 200, 0.05, 0.0075, 0.045

    frame = MoveTo(-OUT/2, -OUT/2).Rectangle(OUT, OUT).Face()
    window = MoveTo(-WIN/2, -WIN/2).Rectangle(WIN, WIN).Face()
    cAp = WorkPlane().Circle(-xin, 0, rw).Face(); cAp.faces.name = "cApos"
    cAn = WorkPlane().Circle(-xout, 0, rw).Face(); cAn.faces.name = "cAneg"
    cBp = WorkPlane().Circle(xin, 0, rw).Face(); cBp.faces.name = "cBpos"
    cBn = WorkPlane().Circle(xout, 0, rw).Face(); cBn.faces.name = "cBneg"
    iron = frame - window; iron.faces.name = "iron"
    win_air = window - cAp - cBp; win_air.faces.name = "air"
    box = MoveTo(-BOX/2, -BOX/2).Rectangle(BOX, BOX).Face(); box.edges.name = "outer"
    out_air = box - frame - cAn - cBn; out_air.faces.name = "air"
    iron.faces.maxh = 0.004; win_air.faces.maxh = 0.004; out_air.faces.maxh = 0.012
    for c in (cAp, cAn, cBp, cBn):
        c.faces.maxh = 0.0012
    mesh = Mesh(OCCGeometry(Glue([iron, win_air, out_air, cAp, cAn, cBp, cBn]), dim=2).GenerateMesh(maxh=0.012))
    MATS = mesh.GetMaterials()

    def nu_of_B(B):
        Bm = sqrt(InnerProduct(B, B) + 1e-20)
        sat = Bm**NEXP / (Bm**NEXP + BK**NEXP)
        inv = 1.0/MUR0 + (1.0 - 1.0/MUR0)*sat
        ir = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in MATS])
        return ir*NU0*inv + (1.0 - ir)*NU0

    def jz(pos, neg, cur):
        j0 = N*cur/(math.pi*rw*rw)
        return CoefficientFunction([{pos: j0, neg: -j0}.get(m, 0.0) for m in MATS])

    def flux(curr):
        A = solve_planar_magnetostatic_nonlinear(mesh, nu_of_B,
                                                 Jz=jz("cApos", "cAneg", curr[0]) + jz("cBpos", "cBneg", curr[1]),
                                                 order=2, relax=0.4, max_iter=200, tol=1e-8)
        return [coil_flux_linkage_2d(A, mesh, "cApos", "cAneg", N, DEPTH),
                coil_flux_linkage_2d(A, mesh, "cBpos", "cBneg", N, DEPTH)]

    with TaskManager():
        L = incremental_inductance_matrix(flux, [1.5, 1.0], 0.02)
        # reciprocity (exact identity, to FD accuracy)
        assert abs(L[0][1] - L[1][0]) / abs(L[0][1]) < 1e-2
        # frozen-perm SECANT apparent inductance is much larger than the incremental
        A0 = solve_planar_magnetostatic_nonlinear(mesh, nu_of_B,
                                                  Jz=jz("cApos", "cAneg", 1.5) + jz("cBpos", "cBneg", 1.0),
                                                  order=2, relax=0.4, max_iter=200, tol=1e-8)
        Aa = solve_planar_magnetostatic(mesh, frozen_reluctivity(A0, nu_of_B), Jz=jz("cApos", "cAneg", 1.0), order=2)
        Lsec_AA = coil_flux_linkage_2d(Aa, mesh, "cApos", "cAneg", N, DEPTH)
        assert Lsec_AA > 2.0 * L[0][0]            # secant >> incremental in saturation (#31)
        # cross inductance collapsed vs the unsaturated linear mutual
        nu_lin = CoefficientFunction([(NU0/MUR0 if m == "iron" else NU0) for m in MATS])
        Au = solve_planar_magnetostatic(mesh, nu_lin, Jz=jz("cBpos", "cBneg", 1.0), order=2)
        Lunsat = coil_flux_linkage_2d(Au, mesh, "cApos", "cAneg", N, DEPTH)
        assert L[0][1] < 0.3 * Lunsat             # mutual collapses under saturation


if __name__ == "__main__":
    validate_incremental_matrix_recovers_linear()
    validate_cross_saturation_fe()
    print("[OK] incremental inductance matrix: reciprocity, secant>>incremental, cross-sat collapse.")
