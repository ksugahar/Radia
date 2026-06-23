r"""Frozen-permeability superposition on a SATURABLE 2-winding iron core -- self-validating.

In a saturated magnetic circuit superposition is broken by nu(|B|).  The frozen-permeability
method restores it: freeze nu at a converged nonlinear operating point, then the linear
frozen-nu problem superposes exactly:

    lambda_A[nonlinear, both windings]  ==  lambda_A[frozen, A only] + lambda_A[frozen, B only]

(exact in the continuum; FE residual ~0.5 %).  NO external data -- the recombination is the
ground truth.  This is the engine behind saturated-machine Ld(id,iq)/Lq maps.
See examples/comsol_class/frozen_permeability.py.
"""
import math
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, grad, sqrt, InnerProduct, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           solve_planar_magnetostatic_nonlinear,
                                           coil_flux_linkage_2d, frozen_reluctivity, NU0)

MUR0, BK, NEXP = 1000.0, 1.4, 6.0
OUT, WIN, BOX, rw = 0.06, 0.03, 0.12, 0.003
N, DEPTH, xin, xout = 200, 0.05, 0.0075, 0.045
_MATS = []


def _nu_of_B(B):
    Bm = sqrt(InnerProduct(B, B) + 1e-20)
    sat = Bm**NEXP / (Bm**NEXP + BK**NEXP)
    inv_mur = 1.0 / MUR0 + (1.0 - 1.0 / MUR0) * sat
    iron = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in _MATS])
    return iron * NU0 * inv_mur + (1.0 - iron) * NU0


def _build():
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


def _jz(pos, neg, current):
    j0 = N * current / (math.pi * rw * rw)
    return CoefficientFunction([{pos: j0, neg: -j0}.get(m, 0.0) for m in _MATS])


def validate_frozen_permeability_superposition():
    global _MATS
    mesh = _build()
    _MATS = mesh.GetMaterials()
    JzA, JzB = _jz("cApos", "cAneg", 1.5), _jz("cBpos", "cBneg", 1.0)
    with TaskManager():
        A0 = solve_planar_magnetostatic_nonlinear(mesh, _nu_of_B, Jz=JzA + JzB, order=2,
                                                  relax=0.4, max_iter=200, tol=1e-8)
        B0 = CoefficientFunction((grad(A0)[1], -grad(A0)[0]))
        mur_op = NU0 / _nu_of_B(B0)(mesh(0.0, 0.024))
        lam_total = coil_flux_linkage_2d(A0, mesh, "cApos", "cAneg", N, DEPTH)
        nu_frozen = frozen_reluctivity(A0, _nu_of_B)
        A_a = solve_planar_magnetostatic(mesh, nu_frozen, Jz=JzA, order=2)
        A_b = solve_planar_magnetostatic(mesh, nu_frozen, Jz=JzB, order=2)
        lam_a = coil_flux_linkage_2d(A_a, mesh, "cApos", "cAneg", N, DEPTH)
        lam_b = coil_flux_linkage_2d(A_b, mesh, "cApos", "cAneg", N, DEPTH)
    err = abs(lam_a + lam_b - lam_total) / abs(lam_total)
    print(f"[frozen-perm] mu_r {MUR0:.0f}->{mur_op:.0f}  lambda_total={lam_total:.5e}  "
          f"A+B={lam_a + lam_b:.5e}  recomb err={100*err:.2f} %")
    assert mur_op < 0.3 * MUR0, f"operating point not nonlinear (mu_r {mur_op:.0f})"
    assert err < 0.02, f"frozen-nu superposition recombination off by {100*err:.2f} %"


if __name__ == "__main__":
    validate_frozen_permeability_superposition()
    print("[OK] frozen-permeability superposition validated (self-consistent).")
