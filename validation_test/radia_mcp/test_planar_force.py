"""2D PLANAR Maxwell-stress FORCE & TORQUE (eggshell), validated against the
exact two-parallel-wire force  F/L = mu0 * I1 * I2 / (2 pi d).

Exercises radia_mcp.radia_ngsolve.force.{eggshell_force_2d, eggshell_torque_2d}
on a field from solve_planar_magnetostatic. Two same-direction out-of-plane
line currents attract; the right wire feels Fx<0, |Fx|=mu0 I^2/(2 pi d) [N/m].
Torque about an offset pivot (0,-h) must equal h*|F| (lever x force).
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
from netgen.occ import OCCGeometry, MoveTo, WorkPlane, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import reluctivity, solve_planar_magnetostatic
from radia_mcp.radia_ngsolve.force import eggshell_force_2d, eggshell_torque_2d

MU0 = 4e-7 * math.pi
I, D, RW, R_FAR = 1000.0, 0.1, 0.005, 2.0
F_EXACT = MU0 * I * I / (2 * math.pi * D)   # |F|/L, attractive


def build_mesh(maxh_w=0.001, maxh_fine=0.0025, maxh_air=0.05):
    w1 = WorkPlane().Circle(+D / 2, 0, RW).Face(); w1.faces.name = "wire1"; w1.maxh = maxh_w
    w2 = WorkPlane().Circle(-D / 2, 0, RW).Face(); w2.faces.name = "wire2"; w2.maxh = maxh_w
    fine = WorkPlane().Circle(+D / 2, 0, 0.045).Face() - w1
    fine.faces.name = "air"; fine.maxh = maxh_fine
    box = MoveTo(-R_FAR, -R_FAR).Rectangle(2 * R_FAR, 2 * R_FAR).Face()
    air = box - w1 - w2 - fine
    air.faces.name = "air"
    for sel in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        sel.name = "outer"
    return Mesh(OCCGeometry(Glue([air, fine, w1, w2]), dim=2).GenerateMesh(maxh=maxh_air))


def main():
    mesh = build_mesh()
    area = math.pi * RW * RW
    Jz = mesh.MaterialCF({"wire1": I / area, "wire2": I / area}, default=0.0)
    with TaskManager():
        gfu = solve_planar_magnetostatic(mesh, reluctivity(mesh, {}), Jz=Jz, order=2)
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))

    Fx, Fy = eggshell_force_2d(B, mesh, (D / 2, 0.0), 0.012, 0.035)
    err_F = (Fx - (-F_EXACT)) / F_EXACT
    print(f"  force: Fx={Fx:+.6f}  Fy={Fy:+.3e} N/m   exact -{F_EXACT:.6f}   err={100*err_F:+.2f}%")

    # torque about pivot (0,-h): r' x F = h*|Fx|  (lever-arm consistency)
    h = 0.1
    tau = eggshell_torque_2d(B, mesh, (D / 2, 0.0), 0.012, 0.035, pivot=(0.0, -h))
    tau_ref = h * F_EXACT
    err_T = (tau - tau_ref) / tau_ref
    print(f"  torque about (0,{-h}): tau={tau:+.6f}  exact {tau_ref:+.6f}  err={100*err_T:+.2f}%")

    assert abs(err_F) < 0.02, f"two-wire force off by {100*err_F:.2f}% (>2%)"
    assert abs(Fy) < 0.05 * F_EXACT, f"transverse Fy should be ~0, got {Fy:.3e}"
    assert abs(err_T) < 0.03, f"torque off by {100*err_T:.2f}% (>3%)"
    print("\n[OK] 2D planar Maxwell-stress force & torque validated vs two-wire "
          f"mu0 I^2/(2 pi d) ({100*abs(err_F):.2f}% force, {100*abs(err_T):.2f}% torque).")


if __name__ == "__main__":
    main()
