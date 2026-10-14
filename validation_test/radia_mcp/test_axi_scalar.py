"""Axisymmetric scalar-potential FEM (FEMM axi electrostatics / heat), validated
on the concentric-SPHERE capacitor / thermal shell, where (full 3D)
    C  = 4 pi eps ab/(b-a)        (capacitance)
    G_th = 4 pi k  ab/(b-a)        (thermal conductance = Q/dT)
and the potential profile is  u(r) = u_in * (a/r) (b-r)/(b-a) ... here checked
via the integral conductance 2W/V^2 = 4 pi c ab/(b-a).

Exercises radia_mcp.radia_ngsolve.scalar_fem2d.{solve_poisson_axi, capacitance_axi}.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, Integrate, grad, dx, TaskManager, x as RC
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue
from radia_mcp.radia_ngsolve.scalar_fem2d import (
    EPS0, solve_poisson_axi, capacitance_axi)

A, B, U0 = 1e-3, 5e-3, 1.0


def build_mesh(maxh=1.2e-4):
    outer = WorkPlane().Circle(0, 0, B).Face()
    inner = WorkPlane().Circle(0, 0, A).Face()
    half = MoveTo(0, -2 * B).Rectangle(2 * B, 4 * B).Face()   # r >= 0, no tangency
    region = (outer - inner) * half
    region.faces.name = "diel"
    s = 0.7071
    for pt in ((B * s, B * s), (B * s, -B * s)):
        region.edges.Nearest(pt).name = "outer"
    for pt in ((A * s, A * s), (A * s, -A * s)):
        region.edges.Nearest(pt).name = "inner"
    return Mesh(OCCGeometry(region, dim=2).GenerateMesh(maxh=maxh))


def main():
    mesh = build_mesh()
    sphere_factor = 4 * math.pi * A * B / (B - A)
    with TaskManager():
        # electrostatics: C = 4 pi eps ab/(b-a)
        eps = CoefficientFunction(EPS0 * 2.0)
        V = solve_poisson_axi(mesh, eps, {"inner": U0, "outer": 0.0})
        C = capacitance_axi(V, mesh, eps, U0)
        C_ex = EPS0 * 2.0 * sphere_factor
        eC = (C - C_ex) / C_ex
        # heat: thermal conductance G_th = 4 pi k ab/(b-a) (same integral, k)
        k = CoefficientFunction(30.0)
        T = solve_poisson_axi(mesh, k, {"inner": U0, "outer": 0.0})
        Gth = 2.0 * (math.pi * Integrate(k * grad(T) * grad(T) * RC * dx, mesh)) / (U0 ** 2)
        Gth_ex = 30.0 * sphere_factor
        eT = (Gth - Gth_ex) / Gth_ex
    print(f"  sphere C  = {C:.5e}  exact {C_ex:.5e}  err {100*eC:+.2f}%")
    print(f"  sphere Gth= {Gth:.5e}  exact {Gth_ex:.5e}  err {100*eT:+.2f}%")
    for nm, e in [("C", eC), ("Gth", eT)]:
        assert abs(e) < 0.01, f"{nm} off by {100*e:.2f}% (>1%)"
    print("\n[OK] axisymmetric scalar FEM (electrostatic / heat) validated on the "
          "concentric-sphere 4 pi c ab/(b-a) to <1%.")


if __name__ == "__main__":
    main()
