r"""Incremental inductance L_inc=dlambda/dI + magnetic co-energy -- regression test.

On the saturating lambda(I) of a closed iron core (smooth nu(B), one winding): the SECANT
L_sec=lambda/I and the INCREMENTAL L_inc=dlambda/dI (central FD, incremental_inductance), with
the co-energy W'=int lambda dI (magnetic_coenergy) and energy W=lambda*I-W'. Robust self-
consistent physics of a CONCAVE lambda(I) (tool-independent): L_inc<L_sec everywhere, both
L_inc and L_sec fall monotonically (saturation), co-energy W'>W once saturating, and the exact
identity W+W'=lambda*I. The differential dq / control inductance + energy-method co-energy.
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
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic_nonlinear,
                                           coil_flux_linkage_2d, incremental_inductance,
                                           magnetic_coenergy, NU0)

MUR0, BK, NEXP = 1000.0, 1.4, 2.0
OUT, WIN, BOX, rw = 0.06, 0.03, 0.12, 0.003
N, DEPTH, xin, xout = 200, 0.05, 0.0075, 0.045
_MATS = None


def _nu_of_B(B):
    Bm = sqrt(InnerProduct(B, B) + 1e-20)
    sat = Bm**NEXP / (Bm**NEXP + BK**NEXP)
    iron = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in _MATS])
    return iron * NU0 * (1.0/MUR0 + (1.0 - 1.0/MUR0)*sat) + (1.0 - iron) * NU0


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


def _jz(I):
    j0 = N * I / (math.pi * rw * rw)
    return CoefficientFunction([{"cApos": j0, "cAneg": -j0}.get(m, 0.0) for m in _MATS])


def validate_incremental_inductance_and_coenergy():
    global _MATS
    mesh = _build(); _MATS = mesh.GetMaterials()
    Is = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2]
    with TaskManager():
        lam = [coil_flux_linkage_2d(solve_planar_magnetostatic_nonlinear(
            mesh, _nu_of_B, Jz=_jz(I), order=2, relax=0.4, max_iter=200, tol=1e-8),
            mesh, "cApos", "cAneg", N, DEPTH) for I in Is]
    Wco = magnetic_coenergy(Is, lam)
    W = [lam[k]*Is[k] - Wco[k] for k in range(len(Is))]
    Lsec = [lam[k]/Is[k] for k in range(len(Is))]
    Linc = {k: incremental_inductance(lam[k+1], lam[k-1], (Is[k+1]-Is[k-1])/2)
            for k in range(1, len(Is)-1)}
    print(f"[incremental L] L_sec={[f'{v:.3e}' for v in Lsec]}")
    print(f"               L_inc={[f'{Linc[k]:.3e}' for k in sorted(Linc)]}")
    # concave lambda(I): tangent below secant everywhere
    assert all(Linc[k] < Lsec[k] for k in Linc), "L_inc must be < L_sec (concave lambda)"
    # both inductances fall monotonically (saturation)
    assert all(Lsec[k] > Lsec[k+1] for k in range(len(Is)-1)), "L_sec must fall (saturation)"
    inc = [Linc[k] for k in sorted(Linc)]
    assert all(inc[i] > inc[i+1] for i in range(len(inc)-1)), "L_inc must fall (saturation)"
    # co-energy exceeds energy once saturating; exact identity W+W'=lambda*I
    assert Wco[-1] > W[-1], "co-energy must exceed energy in saturation"
    assert all(abs((W[k]+Wco[k]) - lam[k]*Is[k]) < 1e-9*max(1e-12, lam[k]*Is[k]) for k in range(len(Is)))


if __name__ == "__main__":
    validate_incremental_inductance_and_coenergy()
    print("[OK] incremental inductance + co-energy validated.")
