"""Nonlinear (saturating-iron) 2D planar A_z magnetostatics, validated by
Ampere's law on a wire wrapped in a nonlinear iron annulus: at radius r in the
iron, H(r) = I/(2 pi r) is INDEPENDENT of the material, so B(r) = BH_curve(H(r))
is an exact reference for any nonlinearity.

Constitutive (Froehlich/Kennelly):  B = H/(alpha + beta H),  nu(B) = alpha/(1-beta B).
Exercises radia_mcp.radia_ngsolve.solve.solve_planar_magnetostatic_nonlinear.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, CoefficientFunction, InnerProduct, sqrt, grad, IfPos,
                     TaskManager)
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic_nonlinear

MU0 = 4e-7 * math.pi
NU0 = 1.0 / MU0
A_W, A_FE, B_FE, R_FAR = 2e-3, 2e-3, 10e-3, 0.06
I = 20.0
MUR0, BSAT = 5000.0, 2.0
ALPHA, BETA = 1.0 / (MU0 * MUR0), 1.0 / BSAT


def B_of_H(H):
    return H / (ALPHA + BETA * H)


def build_mesh(maxh_fe=1.5e-4, maxh_air=4e-3):
    wire = WorkPlane().Circle(0, 0, A_W).Face(); wire.faces.name = "wire"; wire.maxh = maxh_fe
    iron = WorkPlane().Circle(0, 0, B_FE).Face() - WorkPlane().Circle(0, 0, A_FE).Face()
    iron.faces.name = "iron"; iron.maxh = maxh_fe
    box = MoveTo(-R_FAR, -R_FAR).Rectangle(2 * R_FAR, 2 * R_FAR).Face()
    air = box - WorkPlane().Circle(0, 0, B_FE).Face(); air.faces.name = "air"
    for sel in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        sel.name = "outer"
    return Mesh(OCCGeometry(Glue([air, iron, wire]), dim=2).GenerateMesh(maxh=maxh_air))


def main():
    mesh = build_mesh()
    iron_ind = CoefficientFunction([1.0 if m == "iron" else 0.0 for m in mesh.GetMaterials()])
    Jz = mesh.MaterialCF({"wire": I / (math.pi * A_W * A_W)}, default=0.0)

    def nu_of_B(B):
        Bmag = sqrt(InnerProduct(B, B) + 1e-20)
        Bc = IfPos(Bmag - 0.98 * BSAT, 0.98 * BSAT, Bmag)
        return iron_ind * ALPHA / (1.0 - BETA * Bc) + (1.0 - iron_ind) * NU0

    with TaskManager():
        gfu = solve_planar_magnetostatic_nonlinear(mesh, nu_of_B, Jz=Jz, order=2, relax=0.5)
    Bmag = sqrt(InnerProduct(CoefficientFunction((grad(gfu)[1], -grad(gfu)[0])),
                             CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))))
    # sample in the iron, away from the inner interface (r>=4mm) where the
    # tangential-B jump otherwise contaminates point gradient recovery.
    worst = 0.0
    print("  Ampere B(r)=BH(I/2 pi r):")
    for r in (4e-3, 6e-3, 9e-3):
        Bfem = Bmag(mesh(r, 0.0))
        Bref = B_of_H(I / (2 * math.pi * r))
        e = abs(Bfem - Bref) / Bref
        worst = max(worst, e)
        print(f"    r={r*1e3:.0f}mm: B_fem={Bfem:.4f}  B_exact={Bref:.4f}  err={100*e:+.2f}%")
    assert worst < 0.01, f"nonlinear B off by {100*worst:.2f}% (>1%)"
    print(f"\n[OK] nonlinear saturating-iron planar A_z validated by Ampere's law "
          f"(worst {100*worst:.2f}% in the iron).")


if __name__ == "__main__":
    main()
