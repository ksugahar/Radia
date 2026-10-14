"""Helicity diagnostic for volume-current stream functions  (Stage B / B1).

Before trying to extract WINDABLE WIRES from a continuous div-free volume
current J = curl T, ask: is a clean (Clebsch / level-set) wire extraction even
POSSIBLE?  The answer is the CURRENT HELICITY:

    H = integral_Omega  T . J  dV  =  integral T . curl T  dV       (gauge-
        invariant when J.n = 0 on dOmega).

VERIFIED (sympy 2026-06-11): a Clebsch current J = grad(lambda) x grad(mu) has
T = lambda*grad(mu) with T.curl T = 0 -> H = 0.  CONVERSELY a NONZERO-helicity
current (linked / knotted current lines) has NO global Clebsch -> its wires are
NOT clean {lambda=const} INT {mu=const} level-set intersections.

Dimensionless diagnostic (Cauchy-Schwarz bounded in [0,1]):

    H_rel = |integral T . J| / (||T||_L2 ||J||_L2)

  * H_rel ~ 0  : T perp J everywhere (Clebsch-representable) -> wire-extractable.
  * H_rel ~ 1  : maximally helical (Beltrami curl T = k T -> T.J = k|T|^2, so
                 H_rel = 1 exactly) -> linked current lines, no clean Clebsch.

Two analytic checks: an AXIAL T = (0,0,f) has T.J = 0 pointwise (curl is in-plane)
-> H_rel = 0; an ABC/Beltrami flow has curl T = T -> H_rel = 1.

Run:  python validation_test/feec/vim_legacy/helicity_diagnostic.py
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))

from ngsolve import (
    Mesh, HCurl, GridFunction, CF, curl, InnerProduct, dx, sin, cos,
    x, y, z, TaskManager, Integrate,
)
from netgen.csg import CSGeometry, OrthoBrick, Pnt

PI = math.pi


def build_box(maxh=0.25):
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(0, 0, 0), Pnt(1, 1, 1)).mat("dom"))
    return Mesh(geo.GenerateMesh(maxh=maxh))


def helicity_relative(mesh, T_cf, order=3):
    """H_rel = |int T . curl T| / (||T|| ||curl T||) in [0,1] for a CF T."""
    gfT = GridFunction(HCurl(mesh, order=order))
    gfT.Set(T_cf)
    J = curl(gfT)
    H = Integrate(InnerProduct(gfT, J) * dx, mesh)
    tn = math.sqrt(Integrate(InnerProduct(gfT, gfT) * dx, mesh))
    jn = math.sqrt(Integrate(InnerProduct(J, J) * dx, mesh))
    return abs(H) / (tn * jn + 1e-30)


def main():
    print("=== helicity diagnostic (Stage B / B1) ===")
    with TaskManager():
        mesh = build_box()

        # (a) AXIAL T = (0,0,f): curl is in-plane -> T . curl T = 0 -> H_rel = 0
        T_axial = CF((0, 0, sin(PI * x) * sin(PI * y) * sin(PI * z)))
        h_axial = helicity_relative(mesh, T_axial)

        # (b) ABC / Beltrami flow (A=B=C=1, k=2pi on [0,1]): curl T = k T
        #     -> T . curl T = k|T|^2 -> H_rel = 1 (exactly, up to projection)
        k = 2 * PI
        T_abc = CF((sin(k * z) + cos(k * y),
                    sin(k * x) + cos(k * z),
                    sin(k * y) + cos(k * x)))
        h_abc = helicity_relative(mesh, T_abc)

    print(f"   axial T   (Clebsch-type)  H_rel = {h_axial:.3e}   (expect ~0)")
    print(f"   ABC/Beltrami flow         H_rel = {h_abc:.3f}      (expect ~1)")
    print("   -> H_rel ~ 0  => clean Clebsch / level-set wire extraction possible")
    print("   -> H_rel ~ 1  => linked current lines, NO global Clebsch (frontier)")

    ok = (h_axial < 5e-2 and h_abc > 0.9)
    print(f"\n=== PASS: {ok} ===")
    return {"h_axial": float(h_axial), "h_abc": float(h_abc)}, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
