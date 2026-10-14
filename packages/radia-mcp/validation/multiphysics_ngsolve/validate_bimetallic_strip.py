r"""Bimetallic strip (two bonded layers, different alpha -> bending) -- regression test.

Two bonded layers of EQUAL thickness/modulus, different expansion (alpha1, alpha2), uniform dT;
the strip curves with kappa = 3(alpha2-alpha1)dT/(2h) (Timoshenko 1925), cantilever tip
delta = kappa L^2/2 (bimetal_tip_deflection). radia solve_linear_elasticity with a region-wise
alpha (mesh.MaterialCF) reproduces it for a slender strip. Tool-independent (Timoshenko ref).
The thermostat/thermal-switch case -- distinct from the single-material through-thickness gradient.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.elasticity import solve_linear_elasticity, bimetal_tip_deflection

L, h = 0.10, 0.004
E, nu, dT = 70e9, 0.33, 50.0
a1, a2 = 10e-6, 23e-6


def validate_bimetallic_strip_curvature():
    bot = MoveTo(0, 0).Rectangle(L, h/2).Face(); bot.faces.name = "lay1"
    top = MoveTo(0, h/2).Rectangle(L, h/2).Face(); top.faces.name = "lay2"
    strip = Glue([bot, top]); strip.edges.Min(X).name = "clamp"
    mesh = Mesh(OCCGeometry(strip, dim=2).GenerateMesh(maxh=h/4))
    alpha = mesh.MaterialCF({"lay1": a1, "lay2": a2}, default=0.0)
    with TaskManager():
        u = solve_linear_elasticity(mesh, E, nu, dirichlet="clamp",
                                    thermal=(alpha, dT), plane="stress", order=2)
    tip_fe = abs(u(mesh(L, h/2))[1])
    tip_cf = bimetal_tip_deflection(a1, a2, dT, h, L)
    rel = abs(tip_fe - tip_cf) / tip_cf
    print(f"[bimetal] FE={tip_fe*1e3:.4f} mm  Timoshenko={tip_cf*1e3:.4f} mm  ({100*rel:.2f}%)")
    assert rel < 0.04, f"tip {tip_fe:.4e} vs Timoshenko {tip_cf:.4e} ({100*rel:.1f}%)"
    assert tip_fe > 0, "the Delta-alpha must produce a deflection"


if __name__ == "__main__":
    validate_bimetallic_strip_curvature()
    print("[OK] bimetallic strip (Timoshenko curvature) validated.")
