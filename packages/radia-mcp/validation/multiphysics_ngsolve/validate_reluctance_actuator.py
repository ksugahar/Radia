r"""Reluctance-actuator (electromagnet) holding force -- regression test.

Iron window-frame magnetic circuit (coil NI, mu_r, one air gap g); the pole-face
HOLDING FORCE via Maxwell stress  F = B_gap^2 A_face/(2 mu0).  Checks (tool-independent):
  (1) the FE force matches the reluctance model F = mu0(NI)^2 A/[2(g+l_fe/mu_r)^2]
      (magnetic_circuit_gap_force) to a few % (residual = fringing/leakage);
  (2) F*(g+l_fe/mu_r)^2 is constant across gaps -- the 1/g^2 reluctance force law.
Current-driven A_z in SI metres (not scale-invariant). Classic FEMM solenoid/relay task.
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

from ngsolve import Mesh, CoefficientFunction, grad, Integrate, dx, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic,
                                           magnetic_circuit_gap_force, NU0, MU0)

W, ww, b = 0.24, 0.12, 0.06
mu_r, NI, Rbox = 2000.0, 5000.0, 0.60
A_FACE = b * 1.0


def _rect(x0, y0, w, h):
    return MoveTo(x0, y0).Rectangle(w, h).Face()


def _solve_gap(g):
    l_fe = 2.0 * ((W + ww) / 2.0) * 2 - g
    big = _rect(-W/2, -W/2, W, W); window = _rect(-ww/2, -ww/2, ww, ww)
    gap = _rect(-g/2, ww/2, g, b)
    cin = _rect(-ww/2 + 0.005, -0.04, 0.02, 0.08); cout = _rect(-W/2 - 0.025, -0.04, 0.02, 0.08)
    iron = big - window - gap
    gap.faces.name = "gap"; cin.faces.name = "cin"; cout.faces.name = "cout"; iron.faces.name = "iron"
    box = _rect(-Rbox, -Rbox, 2*Rbox, 2*Rbox)
    air = box - iron - gap - cin - cout; air.faces.name = "air"
    for f in (air,):
        f.edges.Max(X).name = "outer"; f.edges.Min(X).name = "outer"
        f.edges.Max(Y).name = "outer"; f.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([air, iron, gap, cin, cout]), dim=2).GenerateMesh(maxh=0.04))
    Ain = Integrate(mesh.MaterialCF({"cin": 1.0}, default=0.0)*dx, mesh)
    Aout = Integrate(mesh.MaterialCF({"cout": 1.0}, default=0.0)*dx, mesh)
    Jz = mesh.MaterialCF({"cin": NI/Ain, "cout": -NI/Aout}, default=0.0)
    nu = mesh.MaterialCF({"iron": NU0/mu_r}, default=NU0)
    with TaskManager():
        gfu = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    Agap = Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*dx, mesh)
    Bx = abs(Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*B[0]*dx, mesh) / Agap)
    F_fe = Bx*Bx*A_FACE/(2*MU0)
    F_model = magnetic_circuit_gap_force(1.0, NI, g, l_fe, mu_r, A_FACE)
    return F_fe, F_model, g + l_fe/mu_r


def validate_reluctance_actuator_holding_force():
    res = [(_solve_gap(g), g) for g in (0.004, 0.006, 0.009)]
    scaled = []
    for (F_fe, F_model, geff), g in res:
        rel = abs(F_fe - F_model) / F_model
        print(f"[reluctance actuator] g={g*1e3:.1f}mm  F_FE={F_fe:.1f}N  model={F_model:.1f}N  "
              f"({100*rel:.2f}%)  F*geff^2={F_fe*geff**2:.4f}")
        assert rel < 0.05, f"g={g}: FE force {F_fe:.1f} vs model {F_model:.1f} ({100*rel:.1f}%)"
        scaled.append(F_fe * geff ** 2)
    # 1/g^2 reluctance force law: F*(g+l_fe/mu_r)^2 constant across gaps
    spread = (max(scaled) - min(scaled)) / (sum(scaled) / len(scaled))
    assert spread < 0.06, f"F*geff^2 not constant across gaps (spread {100*spread:.1f}%) -> not 1/g^2"


if __name__ == "__main__":
    validate_reluctance_actuator_holding_force()
    print("[OK] reluctance-actuator holding force = B^2 A/2mu0, 1/g^2 law validated.")
