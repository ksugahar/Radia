"""Nonlinear permanent magnet (curved demagnetization / field-dependent recoil
permeability) via solve.solve_planar_magnetostatic_nonlinear with a ``magnets``
source. Two checks on the transverse-magnetized cylinder (uniform interior B):

  1. REDUCTION: a CONSTANT recoil mu_r reproduces the linear analytical interior
     field  B_in = mu0 mu_r Hc/(mu_r+1).
  2. SELF-CONSISTENCY: with a field-dependent recoil mu_r(|B|), the converged
     uniform interior B_op must satisfy the SAME relation evaluated at the LOCAL
     recoil permeability,  B_op = mu0 mu_r(B_op) Hc /(mu_r(B_op)+1)  -- exact for
     the uniform-interior cylinder, so it pins the nonlinear-magnet Picard.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, CoefficientFunction, Integrate, grad, sqrt,
                     InnerProduct, TaskManager)
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import (
    reluctivity, solve_planar_magnetostatic_nonlinear, NU0)

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


def mean_By(gfu, mesh):
    vol = Integrate(CoefficientFunction(1.0), mesh, definedon=mesh.Materials("magnet"))
    return Integrate(-grad(gfu)[0], mesh, definedon=mesh.Materials("magnet")) / vol


def mu_r_curve(b):
    """Field-dependent recoil permeability (curved demag): mu_r ~ 3 at B=0,
    falling toward ~1.5 at high B. Works on floats and CoefficientFunctions."""
    return 1.5 + 1.5 / (1.0 + (b / 0.3) ** 2)


def main():
    mesh = build_mesh()
    iron = mesh.MaterialCF({"magnet": 1.0}, default=0.0)

    with TaskManager():
        # 1. constant recoil -> linear analytic interior field
        nu_const = reluctivity(mesh, {"magnet": MU_R})
        g1 = solve_planar_magnetostatic_nonlinear(
            mesh, lambda B: nu_const, magnets={"magnet": (HC, 90.0)},
            relax=0.6, max_iter=60)
        By1 = mean_By(g1, mesh)
        B_in_exact = MU0 * MU_R * HC / (MU_R + 1.0)
        e1 = (By1 - B_in_exact) / B_in_exact
        print(f"  [reduction] <B_y>={By1:.6f} T  linear exact={B_in_exact:.6f}  "
              f"({100*e1:+.3f}%)")
        assert abs(e1) < 5e-3, f"nonlinear-magnet const-mu off {100*e1:.3f}%"

        # 2. field-dependent recoil -> self-consistent operating point
        def nu_of_B(B):
            Bmag = sqrt(InnerProduct(B, B) + 1e-20)
            return iron * (NU0 / mu_r_curve(Bmag)) + (1.0 - iron) * NU0

        g2 = solve_planar_magnetostatic_nonlinear(
            mesh, nu_of_B, magnets={"magnet": (HC, 90.0)}, relax=0.5, max_iter=120)
        B_op = mean_By(g2, mesh)
        mur_op = mu_r_curve(B_op)
        B_pred = MU0 * mur_op * HC / (mur_op + 1.0)
        e2 = (B_op - B_pred) / B_pred
        print(f"  [self-consistent] B_op={B_op:.6f} T  mu_r(B_op)={mur_op:.4f}  "
              f"-> mu0 mu_r Hc/(mu_r+1)={B_pred:.6f}  ({100*e2:+.3f}%)")
        assert abs(e2) < 5e-3, f"nonlinear-magnet self-consistency off {100*e2:.3f}%"

    print(f"\n[OK] nonlinear permanent magnet validated: const recoil -> linear "
          f"B_in ({100*abs(e1):.2f}%); curved recoil self-consistent ({100*abs(e2):.2f}%).")


if __name__ == "__main__":
    main()
