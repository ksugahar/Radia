"""Induction-motor rotor coupling vs slip frequency -- the FEMM induction-motor
method (Meeker, femm.info/wiki/InductionMotorExample): with the rotor at rest the
slip frequency equals the applied frequency, so a time-harmonic solve at the slip
frequency, with a CONDUCTING rotor, develops the induced rotor eddy currents that
do the electromechanical work. This locks in that radia-ngsolve reproduces the IM
rotor branch with EXISTING capability (solve_planar_eddy's complex A_z + j w sigma
eddy term + ohmic_loss_2d) -- no new solver code.

A stranded stator coil drives current; a conducting rotor cylinder develops
induced eddy currents (-j w sigma A_z) that SCREEN flux and dissipate power. Sweep
omega (= slip frequency) and check the induction-coupling signature:
  - Re(lambda) DROPS with slip  (the rotor screens the coil's flux at high slip)
  - Im(lambda) reactive coupling peaks at low/mid slip, decays toward high slip
    (single-pole slip-frequency inductance, Meeker eq. (1))
  - rotor eddy loss P_rotor ~ 0 at slip->0 and rises monotonically with slip
    (the air-gap power that becomes torque x mechanical-speed in a real machine)

This is the parameter-identification quantity Meeker fits (complex inductance vs
slip). The torque-slip curve / circuit-fit (M, Ll, Rr; peak at tau*ws=1) builds on
exactly this solve. See memory reference_femm_motor_fea.md.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, Integrate, dx, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.solve import reluctivity, solve_planar_eddy
from radia_mcp.radia_ngsolve.force import ohmic_loss_2d

SIG = 3.45e7          # aluminium S/m (rotor)
Rr = 0.02             # rotor radius [m]
I = 100.0             # stator coil current [A]


def build_mesh():
    rotor = WorkPlane().Circle(0, 0, Rr).Face(); rotor.faces.name = "rotor"
    coil = MoveTo(0.03, -0.01).Rectangle(0.02, 0.02).Face(); coil.faces.name = "coil"
    box = MoveTo(-0.15, -0.15).Rectangle(0.30, 0.30).Face()
    air = box - rotor - coil; air.faces.name = "air"
    for sel in (air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y), air.edges.Min(Y)):
        sel.name = "outer"
    rotor.faces.maxh = 0.0035; coil.faces.maxh = 0.004
    return Mesh(OCCGeometry(Glue([rotor, coil, air]), dim=2).GenerateMesh(maxh=0.02))


def main():
    mesh = build_mesh()
    nu = reluctivity(mesh, {})                                  # all NU0 (air/rotor)
    sigma = mesh.MaterialCF({"rotor": SIG}, default=0.0)        # only rotor conducts
    Jz = mesh.MaterialCF({"coil": I / (0.02 * 0.02)}, default=0.0)
    omegas = [10.0, 30.0, 100.0, 300.0, 1e3, 3e3, 1e4, 3e4, 1e5]

    reL, imL, P = [], [], []
    with TaskManager():
        for w in omegas:
            Az = solve_planar_eddy(mesh, nu, sigma, w, Jz=Jz)
            lam = Integrate(Az * dx(definedon=mesh.Materials("coil")), mesh)  # complex
            reL.append(lam.real); imL.append(lam.imag)
            P.append(ohmic_loss_2d(-1j * w * Az, mesh, sigma, region="rotor"))
    for w, r, im, p in zip(omegas, reL, imL, P):
        print(f"  w={w:8.0f}  Re(lam)={r:.4e}  Im(lam)={im:.4e}  P_rotor={p:.4e}")

    pmax = max(P)
    print(f"\n  Re(lam): {reL[0]:.3e} -> {reL[-1]:.3e};  P_rotor: {P[0]:.3e} -> {pmax:.3e};"
          f"  Im trough min={min(imL):.3e}")

    # --- induction-coupling signature (robust qualitative checks) ---
    assert reL[-1] < 0.85 * reL[0], "rotor flux-screening (Re lambda drop) not seen"
    assert P[0] < 0.05 * pmax, "rotor loss should be ~0 at low slip"
    assert all(P[i] <= P[i + 1] * 1.05 for i in range(len(P) - 1)), \
        "rotor eddy loss should rise monotonically with slip"
    assert min(imL) < imL[-1], "reactive coupling should peak at low/mid slip (single-pole)"
    print("\n[OK] IM rotor coupling validated: flux screening (Re lam "
          f"{reL[-1]/reL[0]:.2f}x), rotor loss rises from ~0 to {pmax:.2e} W/m, "
          "single-pole reactive coupling -- the FEMM induction-motor slip solve.")


if __name__ == "__main__":
    main()
