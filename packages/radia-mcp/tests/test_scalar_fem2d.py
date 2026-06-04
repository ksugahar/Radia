"""2D scalar-potential FEM family (FEMM csolv / hsolv / current-flow analogs),
all validated on the coaxial annulus where each "conductance per length" equals
    2 pi c / ln(b/a)        (c = eps, sigma, or k)
and the potential profile is  u(r) = u_in * ln(b/r)/ln(b/a).

Exercises radia_mcp.radia_ngsolve.scalar_fem2d.{solve_electrostatic, capacitance,
solve_current_flow, conductance, solve_thermal}.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, Integrate, grad, dx, TaskManager
from netgen.occ import OCCGeometry, WorkPlane
from radia_mcp.radia_ngsolve.scalar_fem2d import (
    EPS0, solve_electrostatic, capacitance,
    solve_current_flow, conductance, solve_thermal)

A, B, U0 = 1e-3, 5e-3, 1.0
LN = math.log(B / A)


def build_mesh(maxh=2e-4):
    ann = WorkPlane().Circle(0, 0, B).Face() - WorkPlane().Circle(0, 0, A).Face()
    ann.faces.name = "mat"
    ann.edges.Nearest((B, 0)).name = "outer"
    ann.edges.Nearest((A, 0)).name = "inner"
    return Mesh(OCCGeometry(ann, dim=2).GenerateMesh(maxh=maxh))


def main():
    mesh = build_mesh()
    with TaskManager():
        # 1) Electrostatics: C/L = 2 pi eps / ln(b/a)
        epsr = 2.5
        eps = CoefficientFunction(EPS0 * epsr)
        V = solve_electrostatic(mesh, eps, {"inner": U0, "outer": 0.0})
        C = capacitance(V, mesh, eps, U0)
        C_ex = 2 * math.pi * EPS0 * epsr / LN
        eC = (C - C_ex) / C_ex
        print(f"  electrostatic C/L = {C:.5e}  exact {C_ex:.5e}  err {100*eC:+.2f}%")

        # 2) DC current flow: G/L = 2 pi sigma / ln(b/a)
        sig = CoefficientFunction(5.8e7)
        Vc = solve_current_flow(mesh, sig, {"inner": U0, "outer": 0.0})
        G = conductance(Vc, mesh, sig, U0)
        G_ex = 2 * math.pi * 5.8e7 / LN
        eG = (G - G_ex) / G_ex
        print(f"  current-flow G/L  = {G:.5e}  exact {G_ex:.5e}  err {100*eG:+.2f}%")

        # 3) Steady heat: thermal conductance/L = 2 pi k / ln(b/a); check T(mid)
        k = CoefficientFunction(50.0)
        T = solve_thermal(mesh, k, {"inner": U0, "outer": 0.0})
        Gth = 2.0 * (0.5 * Integrate(k * grad(T) * grad(T) * dx, mesh)) / (U0 ** 2)
        Gth_ex = 2 * math.pi * 50.0 / LN
        eT = (Gth - Gth_ex) / Gth_ex
        rmid = math.sqrt(A * B)
        Tmid = T(mesh(rmid, 0.0))
        Tmid_ex = U0 * math.log(B / rmid) / LN
        eTp = (Tmid - Tmid_ex) / Tmid_ex
        print(f"  thermal  G/L = {Gth:.5e}  exact {Gth_ex:.5e}  err {100*eT:+.2f}%; "
              f"T(mid)={Tmid:.4f} exact {Tmid_ex:.4f} ({100*eTp:+.2f}%)")

    for nm, e in [("C", eC), ("G", eG), ("Gth", eT), ("Tprofile", eTp)]:
        assert abs(e) < 0.01, f"{nm} off by {100*e:.2f}% (>1%)"
    print("\n[OK] 2D scalar-FEM family (electrostatic / current-flow / thermal) "
          "validated on coaxial 2 pi c / ln(b/a) to <1%.")


if __name__ == "__main__":
    main()
