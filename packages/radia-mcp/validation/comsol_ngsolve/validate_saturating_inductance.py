r"""Saturating winding inductance L(I) -- the Ld saturation knee (saturated-dq foundation).

Closed window-frame iron core (smooth bounded nu(B)), one winding; the SECANT inductance
L_sec(I)=lambda_nl(I)/I via solve.saturated_secant_inductance. Tool-independent self-consistency:
  (1) at LOW current (unsaturated) L_sec -> the linear inductance L0 (nu=nu0 everywhere);
  (2) L_sec DECREASES monotonically with current (the iron saturates);
  (3) deep saturation pulls L_sec well below L0 (the knee).
Frozen-perm superposition for the dq split is separately locked in validate_frozen_permeability.py.
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

from ngsolve import Mesh, CoefficientFunction, sqrt, InnerProduct, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           saturated_secant_inductance,
                                           coil_flux_linkage_2d, NU0)

MUR0, BK, NEXP = 1000.0, 1.4, 6.0
OUT, WIN, BOX, rw = 0.06, 0.03, 0.12, 0.003
N, DEPTH, xin, xout = 200, 0.05, 0.0075, 0.045
_MATS = None


def _nu_of_B(B):
    Bm = sqrt(InnerProduct(B, B) + 1e-20)
    sat = Bm**NEXP / (Bm**NEXP + BK**NEXP)
    iron = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in _MATS])
    return iron * NU0 * (1.0/MUR0 + (1.0 - 1.0/MUR0)*sat) + (1.0 - iron) * NU0


def _nu_lin():
    iron = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in _MATS])
    return iron * NU0/MUR0 + (1.0 - iron) * NU0


def _build():
    frame = MoveTo(-OUT/2, -OUT/2).Rectangle(OUT, OUT).Face()
    window = MoveTo(-WIN/2, -WIN/2).Rectangle(WIN, WIN).Face()
    cAp = WorkPlane().Circle(-xin, 0, rw).Face(); cAp.faces.name = "cApos"
    cAn = WorkPlane().Circle(-xout, 0, rw).Face(); cAn.faces.name = "cAneg"
    iron = frame - window; iron.faces.name = "iron"
    win_air = window - cAp; win_air.faces.name = "air"
    box = MoveTo(-BOX/2, -BOX/2).Rectangle(BOX, BOX).Face(); box.edges.name = "outer"
    out_air = box - frame - cAn; out_air.faces.name = "air"
    iron.faces.maxh = 0.004; win_air.faces.maxh = 0.004; out_air.faces.maxh = 0.012
    cAp.faces.maxh = 0.0012; cAn.faces.maxh = 0.0012
    return Mesh(OCCGeometry(Glue([iron, win_air, out_air, cAp, cAn]), dim=2).GenerateMesh(maxh=0.012))


def _jz(current):
    j0 = N * current / (math.pi * rw * rw)
    return CoefficientFunction([{"cApos": j0, "cAneg": -j0}.get(m, 0.0) for m in _MATS])


def validate_saturating_inductance_knee():
    global _MATS
    mesh = _build(); _MATS = mesh.GetMaterials()
    with TaskManager():
        A_lin = solve_planar_magnetostatic(mesh, _nu_lin(), Jz=_jz(1.0), order=2)
        L0 = coil_flux_linkage_2d(A_lin, mesh, "cApos", "cAneg", N, DEPTH)
        Ls = {I: saturated_secant_inductance(mesh, _nu_of_B, _jz(I), "cApos", "cAneg", N, I,
                                             DEPTH, order=2, relax=0.4, max_iter=200, tol=1e-8)
              for I in (0.1, 1.0, 2.5)}
    print(f"[saturating L] L0={L0:.4e}  L_sec={{{', '.join(f'{i}:{Ls[i]:.4e}' for i in Ls)}}}")
    assert Ls[0.1] / L0 > 0.97, f"low-I L_sec {Ls[0.1]:.3e} should approach linear L0 {L0:.3e}"
    assert Ls[0.1] > Ls[1.0] > Ls[2.5], "L_sec must decrease monotonically (saturation)"
    assert Ls[2.5] / Ls[0.1] < 0.5, f"deep saturation should pull L_sec well below L0 ({Ls[2.5]/Ls[0.1]:.2f})"


if __name__ == "__main__":
    validate_saturating_inductance_knee()
    print("[OK] saturating inductance knee validated.")
