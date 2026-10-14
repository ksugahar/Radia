r"""Thermal bending (through-thickness gradient -> cantilever curvature) -- regression test.

A homogeneous cantilever with a LINEAR through-thickness temperature gradient bends stress-free
into constant curvature kappa = alpha*dT/h; Euler-Bernoulli tip delta = alpha*dT*L^2/(2h)
(thermal_bending_tip_deflection). radia solve_linear_elasticity with a y-linear thermal
pre-strain reproduces it for a slender beam. Tool-independent (beam-theory reference). The
BENDING thermo-mechanical case -- distinct from the AXIAL blocked-bar thermal stress.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, TaskManager, y as Y_cf
from netgen.occ import OCCGeometry, MoveTo, X, Y
from radia_mcp.radia_ngsolve.elasticity import (solve_linear_elasticity,
                                                thermal_bending_tip_deflection)

L, h = 0.10, 0.004
E, nu, alpha, dT = 70e9, 0.33, 23e-6, 50.0


def validate_thermal_bending_cantilever():
    beam = MoveTo(0, 0).Rectangle(L, h).Face(); beam.faces.name = "beam"
    beam.edges.Min(X).name = "clamp"; beam.edges.Max(X).name = "tip"
    mesh = Mesh(OCCGeometry(beam, dim=2).GenerateMesh(maxh=h / 3))
    dT_field = dT * (Y_cf / h - 0.5)
    with TaskManager():
        u = solve_linear_elasticity(mesh, E, nu, dirichlet="clamp",
                                    thermal=(alpha, dT_field), plane="stress", order=2)
    tip_fe = abs(u(mesh(L, h / 2))[1])
    tip_eb = thermal_bending_tip_deflection(alpha, dT, h, L)
    rel = abs(tip_fe - tip_eb) / tip_eb
    print(f"[thermal bending] FE={tip_fe*1e3:.4f} mm  beam-theory={tip_eb*1e3:.4f} mm  ({100*rel:.2f}%)")
    assert rel < 0.04, f"tip {tip_fe:.4e} vs beam theory {tip_eb:.4e} ({100*rel:.1f}%)"
    assert tip_fe > 0, "thermal gradient must produce a deflection"


if __name__ == "__main__":
    validate_thermal_bending_cantilever()
    print("[OK] thermal bending (bimorph curvature) validated.")
