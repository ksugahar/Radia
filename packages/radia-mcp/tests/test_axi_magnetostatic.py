"""Validates solve_axi_magnetostatic (H1Henrotte wrapper) on the uniformly
magnetized linear sphere.

Same geometry and analytical reference as
examples/axifem/research/verification/test_magnetized_sphere.py, but exercises the
reusable API function rather than the raw weak form:

    B_in = 2 mu0 mu_r Hc / (mu_r + 2)     [uniform axial, inside sphere]

For mu_r=2, Hc=3e5 A/m -> B_in_exact = 0.376991 T.  Expected FEM error: -0.05%.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, Integrate, grad, dx, TaskManager
from ngsolve import x as r_cf
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import solve_axi_magnetostatic

MU0 = 4e-7 * math.pi
A_SPHERE = 1.0
MU_R = 2.0
HC = 3.0e5
R_FAR = 50.0
B_IN_EXACT = 2.0 * MU0 * MU_R * HC / (MU_R + 2.0)  # 0.376991 T


def build_mesh(maxh_mag=0.05, maxh_air=0.5):
    disk = WorkPlane().Circle(0, 0, A_SPHERE).Face()
    half = MoveTo(0, -R_FAR).Rectangle(R_FAR, 2 * R_FAR).Face()
    magnet = disk * half
    magnet.faces.name = "magnet"
    magnet.maxh = maxh_mag
    magnet.edges.Min(X).name = "axis"
    air = MoveTo(0, -R_FAR).Rectangle(R_FAR, 2 * R_FAR).Face() - magnet
    air.faces.name = "air"
    air.edges.Min(X).name = "axis"
    air.edges.Max(X).name = "outer"
    air.edges.Max(Y).name = "outer"
    air.edges.Min(Y).name = "outer"
    return Mesh(OCCGeometry(Glue([air, magnet]), dim=2).GenerateMesh(maxh=maxh_air))


def main():
    mesh = build_mesh()
    nu = CoefficientFunction(
        1.0 / (MU0 * mesh.MaterialCF({"magnet": MU_R}, default=1.0)))

    with TaskManager():
        gfu = solve_axi_magnetostatic(
            mesh, nu,
            magnets={"magnet": (HC, 90.0)},  # theta=90 => axial (+z) magnetization
            order=2)

    # B_z = grad(A)[0] + A/r  (Henrotte convention: r = NGSolve x)
    B_z = grad(gfu)[0] + gfu / r_cf

    # volume-averaged B_z over the magnet (robust off-axis estimate)
    vol = Integrate(r_cf, mesh, definedon=mesh.Materials("magnet"))
    bz_avg = Integrate(r_cf * B_z, mesh, definedon=mesh.Materials("magnet")) / vol
    err = (bz_avg - B_IN_EXACT) / B_IN_EXACT

    print(f"  <B_z>_magnet = {bz_avg:.6f} T  exact {B_IN_EXACT:.6f}  err {100*err:+.3f}%")
    assert abs(err) < 5e-3, f"B_z off by {100*err:.3f}% (>0.5%)"
    print(f"\n[OK] solve_axi_magnetostatic validated on magnetized sphere "
          f"({100*abs(err):.3f}% error vs analytical).")


if __name__ == "__main__":
    main()
