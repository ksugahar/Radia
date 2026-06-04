"""2D inductance via the energy method L=2W/I^2 (force.inductance_2d), validated
on the classic radius-independent round-wire INTERNAL inductance
    L_int = mu0/(8 pi) = 5.0e-8 H/m
(energy of the field inside the conductor, B(r)=mu0 I r/(2 pi a^2))."""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import reluctivity, solve_planar_magnetostatic
from radia_mcp.radia_ngsolve.force import inductance_2d

MU0 = 4e-7 * math.pi
A_W, R_FAR, I = 1e-3, 30e-3, 100.0
L_INT_EXACT = MU0 / (8 * math.pi)


def build_mesh(maxh_w=4e-5, maxh_air=3e-3):
    wire = WorkPlane().Circle(0, 0, A_W).Face(); wire.faces.name = "wire"; wire.maxh = maxh_w
    box = MoveTo(-R_FAR, -R_FAR).Rectangle(2 * R_FAR, 2 * R_FAR).Face()
    air = box - wire; air.faces.name = "air"
    for sel in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        sel.name = "outer"
    return Mesh(OCCGeometry(Glue([air, wire]), dim=2).GenerateMesh(maxh=maxh_air))


def main():
    mesh = build_mesh()
    nu = reluctivity(mesh, {})
    Jz = mesh.MaterialCF({"wire": I / (math.pi * A_W * A_W)}, default=0.0)
    with TaskManager():
        gfu = solve_planar_magnetostatic(mesh, nu, Jz=Jz, order=2)
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    L_int = inductance_2d(B, mesh, nu, I, region="wire")
    err = (L_int - L_INT_EXACT) / L_INT_EXACT
    print(f"  L_int FEM = {L_int:.6e} H/m  exact mu0/8pi = {L_INT_EXACT:.6e}  err = {100*err:+.2f}%")
    assert abs(err) < 0.01, f"internal inductance off by {100*err:.2f}% (>1%)"
    print("\n[OK] 2D inductance (energy method) validated vs round-wire "
          f"L_int = mu0/8pi ({100*abs(err):.2f}%).")


if __name__ == "__main__":
    main()
