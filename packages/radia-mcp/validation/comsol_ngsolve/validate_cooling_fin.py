r"""Straight cooling fin (extended surface) -- regression test (#47).

Thin fin, base dT_base, convective (h) faces, adiabatic tip: dT(x)=dT_base cosh(m(L-x))/cosh(mL),
m=sqrt(2h/(kt)), efficiency eta=tanh(mL)/(mL). radia's Robin heat solver with a base Dirichlet
reproduces the profile, tip and efficiency; the 1-D fin closed form is the tool-independent gate
(a reference commercial FE code matches to ~0.14% live, recorded internally). The
conduction-convection / heat-sink block; generalises the convective slab (#44)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import (fin_parameter, fin_tip_temperature_rise,
                                                  fin_efficiency)

L, T, K, DTB = 0.05, 0.002, 200.0, 100.0


def validate_fin_helpers():
    m = fin_parameter(500.0, K, T)
    assert math.isclose(m, math.sqrt(2 * 500.0 / (K * T)), rel_tol=1e-12)   # = 50
    assert math.isclose(m, 50.0, rel_tol=1e-12)
    # efficiency limits: short/fat fin (mL->0) -> 1 ; long/thin (mL->inf) -> 0
    assert abs(fin_efficiency(m, 1e-4) - 1.0) < 1e-3
    assert fin_efficiency(m, 1.0) < 0.05            # mL=50 -> ~1/50
    assert fin_efficiency(m, 0.01) > fin_efficiency(m, 0.05)   # monotone decreasing in L
    # tip rise = dT_base/cosh(mL); falls as the fin lengthens
    assert math.isclose(fin_tip_temperature_rise(DTB, m, L), DTB / math.cosh(m * L), rel_tol=1e-12)
    assert fin_tip_temperature_rise(DTB, m, 2 * L) < fin_tip_temperature_rise(DTB, m, L)


def validate_cooling_fin_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction, Integrate, ds, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane, X, Y
    from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady_robin

    fin = WorkPlane().MoveTo(0, -T/2).Rectangle(L, T).Face(); fin.faces.name = "fin"
    fin.edges.Min(X).name = "base"
    fin.edges.Max(Y).name = "conv"; fin.edges.Min(Y).name = "conv"
    fin.edges.Max(X).name = "tip"
    mesh = Mesh(OCCGeometry(fin, dim=2).GenerateMesh(maxh=T / 2))

    h = 500.0
    m = fin_parameter(h, K, T)
    gT = solve_heat_steady_robin(mesh, CoefficientFunction(0.0), K, "fin", "conv", h,
                                 order=3, dirichlet="base", dirichlet_value=DTB)
    # profile dT(x) = dT_base cosh(m(L-x))/cosh(mL)
    for xr in (0.25, 0.5, 0.75):
        x = xr * L
        dT_an = DTB * math.cosh(m * (L - x)) / math.cosh(m * L)
        assert abs(gT(mesh(x, 0.0)) - dT_an) / dT_an < 5e-3, f"x/L={xr}"
    # tip rise
    tip_an = fin_tip_temperature_rise(DTB, m, L)
    assert abs(gT(mesh(L * 0.999, 0.0)) - tip_an) / tip_an < 5e-3
    # fin efficiency from the convected heat == tanh(mL)/(mL)
    Qconv = Integrate(h * gT * ds(definedon=mesh.Boundaries("conv")), mesh)
    eta_fe = Qconv / (h * (2 * L) * DTB)
    assert abs(eta_fe - fin_efficiency(m, L)) / fin_efficiency(m, L) < 5e-3


if __name__ == "__main__":
    validate_fin_helpers()
    validate_cooling_fin_fe()
    print("[OK] cooling fin: cosh profile + tip + efficiency tanh(mL)/(mL) validated.")
