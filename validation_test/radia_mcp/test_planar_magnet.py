"""2D PLANAR (Cartesian) permanent-magnet capability -- FEMM prob1big analog
on standard NGSolve H1, validated against the analytical magnetized cylinder.

Exercises radia_mcp.radia_ngsolve.solve.solve_planar_magnetostatic (+
planar_magnet_source). A transverse-magnetized linear cylinder (radius a,
mu_r, Hc) has UNIFORM internal B and a 2D line-dipole outside:
    B_in = mu0 * mu_r * Hc / (mu_r + 1)        (2D transverse demag = 1/2)
    B_equator(r) = -B_in * (a/r)^2             (perp-to-M axis)
For mu_r=2, Hc=3e5 -> B_in = 0.251327 T.
"""
import math
import os
import sys

# Test the DEV-repo source (radia_mcp is installed non-editable in site-packages).
_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, CoefficientFunction, Integrate, dx, TaskManager,
                     grad, x as X_cf)
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import (
    reluctivity, solve_planar_magnetostatic)

MU0 = 4e-7 * math.pi
A_CYL, MU_R, HC, R_FAR = 1.0, 2.0, 3.0e5, 50.0
B_IN_EXACT = MU0 * MU_R * HC / (MU_R + 1.0)


def build_mesh(maxh_mag=0.05, maxh_air=0.5):
    magnet = WorkPlane().Circle(0, 0, A_CYL).Face()
    magnet.faces.name = "magnet"
    magnet.maxh = maxh_mag
    box = MoveTo(-R_FAR, -R_FAR).Rectangle(2 * R_FAR, 2 * R_FAR).Face()
    air = box - magnet
    air.faces.name = "air"
    for sel in (air.edges.Max(X), air.edges.Min(X),
                air.edges.Max(Y), air.edges.Min(Y)):
        sel.name = "outer"
    return Mesh(OCCGeometry(Glue([air, magnet]), dim=2).GenerateMesh(maxh=maxh_air))


def main():
    mesh = build_mesh()
    nu = reluctivity(mesh, {"magnet": MU_R})
    with TaskManager():
        gfu = solve_planar_magnetostatic(
            mesh, nu, magnets={"magnet": (HC, 90.0)}, order=2)  # +y magnetization
    B_x = grad(gfu)[1]
    B_y = -grad(gfu)[0]

    vol = Integrate(CoefficientFunction(1.0), mesh, definedon=mesh.Materials("magnet"))
    by_avg = Integrate(B_y, mesh, definedon=mesh.Materials("magnet")) / vol
    bx_avg = Integrate(B_x, mesh, definedon=mesh.Materials("magnet")) / vol
    err_in = (by_avg - B_IN_EXACT) / B_IN_EXACT
    print(f"  <B_y>_magnet = {by_avg:.6f} T  (exact {B_IN_EXACT:.6f})  err = {100*err_in:+.3f}%")
    print(f"  <B_x>_magnet = {bx_avg:+.2e} T  (expect ~0)")

    ext_ok = True
    print("  external equatorial B_y vs 2D dipole -B_in*(a/r)^2:")
    for rr in (3.0, 5.0):
        by = B_y(mesh(rr, 0.0))
        ref = -B_IN_EXACT * (A_CYL / rr) ** 2
        rel = abs(by - ref) / abs(ref)
        print(f"    r={rr:.1f}: B_y={by:+.6e}  ref={ref:+.6e}  rel={100*rel:.2f}%")
        ext_ok = ext_ok and rel < 0.03

    assert abs(err_in) < 5e-3, f"internal B_y off by {100*err_in:.3f}% (>0.5%)"
    assert abs(bx_avg) < 1e-3, f"interior B_x should be ~0, got {bx_avg:.3e}"
    assert ext_ok, "external 2D dipole off by >3%"
    print("\n[OK] 2D planar permanent-magnet (FEMM prob1big analog) validated on "
          f"standard NGSolve H1 (<B_y> within {100*abs(err_in):.3f}% of analytical).")


if __name__ == "__main__":
    main()
