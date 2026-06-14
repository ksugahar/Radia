"""2D magnetostatic vector-potential A_z member of the scalar_fem2d family, validated
on the SOLID-CONDUCTOR COAX where the per-length self-inductance has the closed form

    L' = mu0/(2 pi) * ( mu_r * ln(b/a) + 1/4 )

    external flux (a<r<b) -> mu0 mu_r/(2 pi) ln(b/a)   [annulus reluctivity]
    internal distributed-current energy (r<a) -> mu0/(8 pi) = mu0/(2 pi)*(1/4)

This is the CURRENT-carrying (vector-potential) complement to the current-free magnetic
scalar potential (mfnc) already in scalar_fem2d: same elliptic operator
-div(nu grad u) = f, here u = A_z, f = J_z, nu = 1/(mu0 mu_r). The energy-method
inductance L' = 2W/I^2 with a *frozen* nu field is the linear limit of frozen-permeability
machine Ld/Lq (the nonlinear converged-nu case is test_nonlinear_magnetostatic_newton).

Exercises radia_mcp.radia_ngsolve.scalar_fem2d.{solve_magnetostatic_az, inductance_2d}.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, Glue
from radia_mcp.radia_ngsolve.scalar_fem2d import (
    solve_magnetostatic_az, inductance_2d)

MU0 = 4e-7 * math.pi
A_R, B_R, I = 2e-3, 8e-3, 1.0          # inner conductor / shield radii [m]; total current [A]
JZ = I / (math.pi * A_R * A_R)         # uniform conductor current density [A/m^2]
LN = math.log(B_R / A_R)


def build_coax(maxh):
    cond = WorkPlane().Circle(0, 0, A_R).Face()
    cond.faces.name = "cond"
    diel = WorkPlane().Circle(0, 0, B_R).Face() - WorkPlane().Circle(0, 0, A_R).Face()
    diel.faces.name = "diel"
    shape = Glue([cond, diel])
    shape.edges.Nearest((B_R, 0.0)).name = "outer"
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh))


def solve_L(maxh, mur, order=3):
    mesh = build_coax(maxh)
    nu = mesh.MaterialCF({"cond": 1.0 / MU0, "diel": 1.0 / (MU0 * mur)},
                         default=1.0 / MU0)
    A = solve_magnetostatic_az(mesh, nu, {"cond": JZ}, {"outer": 0.0}, order=order)
    return inductance_2d(A, mesh, nu, I)


def main():
    with TaskManager():
        # 1) air coax: convergence toward mu0/(2pi)(ln(b/a)+1/4)
        Lref1 = MU0 / (2 * math.pi) * (LN + 0.25)
        L_coarse = solve_L(0.5e-3, 1.0)
        L_fine = solve_L(0.2e-3, 1.0)
        e_coarse = abs(L_coarse - Lref1) / Lref1
        e_fine = abs(L_fine - Lref1) / Lref1
        print(f"  air coax  L'(coarse)={L_coarse:.5e} err {100*e_coarse:+.2f}%")
        print(f"  air coax  L'(fine)  ={L_fine:.5e} err {100*e_fine:+.2f}%  "
              f"(exact {Lref1:.5e})")

        # 2) magnetic-fill annulus (spatially varying nu = frozen-permeability essence)
        mur = 50.0
        Lref2 = MU0 / (2 * math.pi) * (mur * LN + 0.25)
        L_mu = solve_L(0.2e-3, mur)
        e_mu = abs(L_mu - Lref2) / Lref2
        print(f"  mur={mur:.0f}    L'={L_mu:.5e} err {100*e_mu:+.2f}%  (exact {Lref2:.5e})")

    assert e_fine < 5e-3, f"air-coax L' off by {100*e_fine:.2f}% (>0.5%)"
    assert e_fine < e_coarse, "L' did not converge (fine err >= coarse err)"
    assert e_mu < 5e-3, f"magnetic-fill L' off by {100*e_mu:.2f}% (>0.5%)"
    print("\n[OK] 2D magnetostatic A_z + inductance_2d validated on the solid-conductor "
          "coax L' = mu0/(2pi)(mu_r ln(b/a)+1/4) to <0.5% (air convergent; mu_r-fill nu-varying).")


if __name__ == "__main__":
    main()
