"""COMSOL-class #23 -- FROZEN-PERMEABILITY superposition (saturated-machine Ld/Lq engine).

In a SATURATED magnetic circuit the flux linkage is NOT a linear sum of its sources
(superposition broken by nu(|B|)), so you cannot separate PM / d-axis / q-axis flux by
differencing.  The frozen-permeability method restores superposition: take a CONVERGED
nonlinear solve, FREEZE nu at its operating-point field, then the LINEAR frozen-nu
problem superposes exactly -- the engine behind saturated Ld(id,iq)/Lq maps and MTPA.

Self-validating, NO external data: with nu frozen at the operating point,

    lambda_A[nonlinear, both windings]  ==  lambda_A[frozen, A only] + lambda_A[frozen, B only]

(exact in the continuum; FE gives a small residual).  Closed window-frame iron core
(flux circulates) with two windings each threading a limb; meshed in METRES.

Two lessons baked in (see frozen_reluctivity docstring): freeze the DIRECT CF (L2-projecting
a clamped nu smears it); use a SMOOTH BOUNDED nu(B) so Picard converges (the freeze only
reproduces a converged solve -- which is the built-in correctness check).
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, sqrt, InnerProduct, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           solve_planar_magnetostatic_nonlinear,
                                           coil_flux_linkage_2d, frozen_reluctivity, NU0)

MUR0, BK, NEXP = 1000.0, 1.4, 6.0       # smooth bounded iron: mu_r0 (soft) -> 1 (saturated)
OUT, WIN, BOX, rw = 0.06, 0.03, 0.12, 0.003
N, DEPTH, xin, xout = 200, 0.05, 0.0075, 0.045
MATS = None


def nu_of_B(B):
    Bm = sqrt(InnerProduct(B, B) + 1e-20)
    sat = Bm**NEXP / (Bm**NEXP + BK**NEXP)
    inv_mur = 1.0 / MUR0 + (1.0 - 1.0 / MUR0) * sat
    iron = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in MATS])
    return iron * NU0 * inv_mur + (1.0 - iron) * NU0


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
    return Mesh(OCCGeometry(Glue([iron, win_air, out_air, cAp, cAn, cBp, cBn]),
                            dim=2).GenerateMesh(maxh=0.012))


def jz(pos, neg, current):
    j0 = N * current / (math.pi * rw * rw)
    return CoefficientFunction([{pos: j0, neg: -j0}.get(m, 0.0) for m in MATS])


mesh = build()
MATS = mesh.GetMaterials()
JzA, JzB = jz("cApos", "cAneg", 1.5), jz("cBpos", "cBneg", 1.0)

with TaskManager():
    A0 = solve_planar_magnetostatic_nonlinear(mesh, nu_of_B, Jz=JzA + JzB, order=2,
                                              relax=0.4, max_iter=200, tol=1e-8)
    B0 = CoefficientFunction((grad(A0)[1], -grad(A0)[0]))
    mur_op = NU0 / nu_of_B(B0)(mesh(0.0, 0.024))            # iron mu_r at the operating point
    lam_total = coil_flux_linkage_2d(A0, mesh, "cApos", "cAneg", N, DEPTH)

    nu_frozen = frozen_reluctivity(A0, nu_of_B)             # FROZEN at the operating point
    A_a = solve_planar_magnetostatic(mesh, nu_frozen, Jz=JzA, order=2)
    A_b = solve_planar_magnetostatic(mesh, nu_frozen, Jz=JzB, order=2)
    lam_a = coil_flux_linkage_2d(A_a, mesh, "cApos", "cAneg", N, DEPTH)
    lam_b = coil_flux_linkage_2d(A_b, mesh, "cApos", "cAneg", N, DEPTH)

recomb = lam_a + lam_b
err = abs(recomb - lam_total) / abs(lam_total) * 100
print(f"operating point: iron mu_r {MUR0:.0f} -> {mur_op:.0f} (nonlinear)")
print(f"nonlinear lambda_A         = {lam_total:.6e}")
print(f"frozen A-only + B-only     = {recomb:.6e}  ({lam_a:.4e} + {lam_b:.4e})")
print(f"recombination error        = {err:.2f} %")
print("OK" if err < 2.0 else f"OFF ({err:.2f}%)")
