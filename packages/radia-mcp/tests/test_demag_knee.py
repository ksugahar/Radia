# -*- coding: utf-8 -*-
"""Per-element demag-knee post-processor (radia_ngsolve.demag) validated on the
transverse-magnetized cylinder, where the interior operating point is UNIFORM
and analytic:

    B_op = mu0*mu_r*Hc/(mu_r+1)
    H_m  = B_op/(mu0*mu_r) - Hc = -mu_r*Hc/(mu_r+1)

For mu_r=2, Hc=3e5: B_op = 0.2513 T, H_m = -2.0e5 A/m. The knee check must flag
the magnet as UNSAFE when the knee field is LESS negative than the operating
H_m (knee above the operating point => already past the knee), and SAFE when
the knee is more negative. This is the radia-ngsolve answer to FEMM's manual
"is every point above the knee?" check (cross-learning capability).
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import reluctivity, solve_planar_magnetostatic
from radia_mcp.radia_ngsolve.demag import demag_operating_field, demag_knee_report

MU0 = 4e-7 * math.pi
A_CYL, MU_R, HC, R_FAR = 1.0, 2.0, 3.0e5, 50.0


def build_mesh():
    magnet = WorkPlane().Circle(0, 0, A_CYL).Face()
    magnet.faces.name = "magnet"; magnet.maxh = 0.05
    box = MoveTo(-R_FAR, -R_FAR).Rectangle(2 * R_FAR, 2 * R_FAR).Face()
    air = box - magnet; air.faces.name = "air"
    for s in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        s.name = "outer"
    return Mesh(OCCGeometry(Glue([air, magnet]), dim=2).GenerateMesh(maxh=0.5))


def main():
    mesh = build_mesh()
    with TaskManager():
        nu = reluctivity(mesh, {"magnet": MU_R})
        gfu = solve_planar_magnetostatic(mesh, nu, magnets={"magnet": (HC, 90.0)})
        # knee ABOVE the operating point (-1e5 > -2e5): magnet is past the knee
        rep_unsafe = demag_knee_report(gfu, mesh, {"magnet": (HC, 90.0)}, MU_R, H_knee=-1.0e5)
        # knee BELOW the operating point (-3e5 < -2e5): magnet is safe
        rep_safe = demag_knee_report(gfu, mesh, {"magnet": (HC, 90.0)}, MU_R, H_knee=-3.0e5)

    Hm_exact = -MU_R * HC / (MU_R + 1.0)          # -2.0e5 A/m
    Bop_exact = MU0 * MU_R * HC / (MU_R + 1.0)    # 0.2513 T
    eH = (rep_unsafe["mean_Hm"] - Hm_exact) / abs(Hm_exact)
    eB = (rep_unsafe["mean_Bm"] - Bop_exact) / abs(Bop_exact)

    print(f"  mean_Hm = {rep_unsafe['mean_Hm']:.4e} A/m  (exact {Hm_exact:.4e}, {100*eH:+.2f}%)")
    print(f"  mean_Bm = {rep_unsafe['mean_Bm']:.6f} T    (exact {Bop_exact:.6f}, {100*eB:+.2f}%)")
    print(f"  knee=-1e5: frac_below={rep_unsafe['area_fraction_below_knee']:.3f}  safe={rep_unsafe['safe']}")
    print(f"  knee=-3e5: frac_below={rep_safe['area_fraction_below_knee']:.3f}  safe={rep_safe['safe']}")

    assert abs(eH) < 0.03, f"H_m off {100*eH:.2f}%"
    assert abs(eB) < 0.03, f"B_m off {100*eB:.2f}%"
    # operating H_m = -2e5: knee at -1e5 -> operating point is below the knee -> unsafe
    assert rep_unsafe["area_fraction_below_knee"] > 0.9 and not rep_unsafe["safe"]
    # knee at -3e5 -> operating point (-2e5) is above the knee -> safe
    assert rep_safe["area_fraction_below_knee"] < 0.1 and rep_safe["safe"]

    print(f"\n[OK] demag-knee post-processor: per-element H_m = {Hm_exact:.3e} A/m "
          f"(<3% err), distribution flips safe/unsafe across the knee.")


if __name__ == "__main__":
    main()
