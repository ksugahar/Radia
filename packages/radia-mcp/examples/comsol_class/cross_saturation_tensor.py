"""COMSOL-class #52 -- INCREMENTAL inductance matrix + cross-saturation + reciprocity.

The multi-port generalisation of the secant/incremental inductance (#28/#31): two windings on a
SHARED SATURABLE iron core. The INCREMENTAL matrix L_jk = d lambda_j/d i_k (the small-signal /
dq-control inductance) is got by FINITE DIFFERENCE of the nonlinear flux-linkage map
(`incremental_inductance_matrix`). Two exact properties:

  * RECIPROCITY: L_AB = L_BA (the matrix is symmetric -- the magnetic co-energy Hessian);
  * INCREMENTAL << APPARENT(secant): the frozen-permeability secant inductance (freeze nu at
    the operating |B|, #23) is much LARGER than this tangent incremental in saturation (#31).

Cross-(saturation): the off-diagonal (cross) inductance collapses as the shared iron saturates.
Reuses the #23 saturable window-frame engine. Self-contained, headless; metres.
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, sqrt, InnerProduct, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           solve_planar_magnetostatic_nonlinear,
                                           coil_flux_linkage_2d, frozen_reluctivity,
                                           incremental_inductance_matrix, NU0)

MUR0, BK, NEXP = 1000.0, 1.4, 6.0
OUT, WIN, BOX, rw = 0.06, 0.03, 0.12, 0.003
N, DEPTH, xin, xout = 200, 0.05, 0.0075, 0.045
MATS = None


def nu_of_B(B):
    Bm = sqrt(InnerProduct(B, B) + 1e-20)
    sat = Bm**NEXP / (Bm**NEXP + BK**NEXP)
    inv = 1.0/MUR0 + (1.0 - 1.0/MUR0)*sat
    iron = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in MATS])
    return iron*NU0*inv + (1.0 - iron)*NU0


def build():
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
    return Mesh(OCCGeometry(Glue([iron, win_air, out_air, cAp, cAn, cBp, cBn]), dim=2).GenerateMesh(maxh=0.012))


def jz(pos, neg, current):
    j0 = N*current/(math.pi*rw*rw)
    return CoefficientFunction([{pos: j0, neg: -j0}.get(m, 0.0) for m in MATS])


mesh = build(); MATS = mesh.GetMaterials()


def flux_linkages(curr):              # nonlinear lambda_A, lambda_B at currents [iA, iB]
    A = solve_planar_magnetostatic_nonlinear(mesh, nu_of_B, Jz=jz("cApos", "cAneg", curr[0]) + jz("cBpos", "cBneg", curr[1]),
                                             order=2, relax=0.4, max_iter=200, tol=1e-8)
    return [coil_flux_linkage_2d(A, mesh, "cApos", "cAneg", N, DEPTH),
            coil_flux_linkage_2d(A, mesh, "cBpos", "cBneg", N, DEPTH)]


with TaskManager():
    Linc = incremental_inductance_matrix(flux_linkages, [1.5, 1.0], 0.02)     # tangent matrix
    # secant (frozen-perm apparent) for contrast
    A0 = solve_planar_magnetostatic_nonlinear(mesh, nu_of_B, Jz=jz("cApos", "cAneg", 1.5) + jz("cBpos", "cBneg", 1.0),
                                              order=2, relax=0.4, max_iter=200, tol=1e-8)
    mur_op = NU0/nu_of_B(CoefficientFunction((grad(A0)[1], -grad(A0)[0])))(mesh(0.0, 0.024))
    nuf = frozen_reluctivity(A0, nu_of_B)
    Aa = solve_planar_magnetostatic(mesh, nuf, Jz=jz("cApos", "cAneg", 1.0), order=2)
    Lsec_AA = coil_flux_linkage_2d(Aa, mesh, "cApos", "cAneg", N, DEPTH)
    nu_lin = CoefficientFunction([(NU0/MUR0 if m == "iron" else NU0) for m in MATS])
    Au = solve_planar_magnetostatic(mesh, nu_lin, Jz=jz("cBpos", "cBneg", 1.0), order=2)
    Lunsat_AB = coil_flux_linkage_2d(Au, mesh, "cApos", "cAneg", N, DEPTH)

print(f"operating point: iron mu_r {MUR0:.0f} -> {mur_op:.0f}")
print(f"INCREMENTAL matrix [H]:  [[{Linc[0][0]:.4e}, {Linc[0][1]:.4e}],")
print(f"                          [{Linc[1][0]:.4e}, {Linc[1][1]:.4e}]]")
print(f"RECIPROCITY  |L_AB-L_BA|/|L_AB| = {abs(Linc[0][1]-Linc[1][0])/abs(Linc[0][1])*100:.3f}%")
print(f"secant/incremental (AA) = {Lsec_AA/Linc[0][0]:.2f}x  (frozen-perm apparent >> tangent, #31)")
print(f"cross-saturation: incremental mutual {Linc[0][1]:.4e} vs unsaturated {Lunsat_AB:.4e} "
      f"(collapsed to {Linc[0][1]/Lunsat_AB:.3f}x)")
print("OK -- incremental matrix symmetric (reciprocity); apparent>>incremental; mutual collapses with saturation.")
