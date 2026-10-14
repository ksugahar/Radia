# -*- coding: utf-8 -*-
r"""End-to-end: the promoted build123d EM archetypes (ipm_rotor + slotted_stator) and the AGE
rotating-machine solver, driven by ONE shared parameter set.

Two halves, one machine:
  CAD  -- ipm_rotor / slotted_stator build the manufacturable 3D geometry: labelled multi-region
          solids (iron + buried magnets with radial, N/S-alternating easy-axis labels), Netgen-meshable.
  SOLVE-- an Air-Gap-Element (AGE) magnetostatic solve of the SAME machine (same pole count, Br) on a
          fast 2D concentric-ring cross-section verifies the synchronous torque physics: phase-locked
          stator current gives CONSTANT torque at quadrature, and the torque-angle law is T(delta) =
          T_max*sin(delta) (the non-salient baseline; a flux-barrier mesh would add the IPM reluctance
          sin(2*delta) term -- out of scope here).

The thread is the shared (N_POLES, BR, radii): the archetype is the CAD of the rotor whose torque the
AGE computes.  Mirrors test_age_synchronous_torque's verified AGE pattern, with the geometry now coming
from the formal archetype API.
"""
import math
import os
import sys

import numpy as np

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.build123d.archetypes import ipm_rotor, slotted_stator, parse_magnetization

mm = 1e-3
N_POLES = 2                                                  # 1 pole pair (p=1) -- matches the AGE source
BR = 1.2
R0, Rr, Rm, Rs, Rw, REXT = 3 * mm, 8 * mm, 11 * mm, 11.5 * mm, 12.5 * mm, 22 * mm
MUR_FE, MUR_PM = 5.0e4, 1.05
MU0 = 4e-7 * math.pi
LSTK, J0 = 0.05, 6.0e6
HARM = [1, 3, 5]


# ---- CAD half: the promoted archetypes build a valid, labelled, meshable machine -------------------
def test_ipm_rotor_and_stator_cad():
    rotor = ipm_rotor(R0 * 1e3, Rr * 1e3, N_POLES, 3.0, 1.0, 40.0, 5.0, name="pm")
    assert len(rotor.children) == 1 + 2 * N_POLES and all(s.is_valid for s in rotor.solids())
    mags = [c for c in rotor.children if c.label != "rotor_iron"]
    assert len(mags) == 2 * N_POLES
    assert all(math.isfinite(parse_magnetization(m.label)) for m in mags), "magnets carry easy-axis labels"
    # alternation is only visible for >= 4 poles (a 2-pole IPM is diametric -> one Cartesian angle):
    quad = ipm_rotor(R0 * 1e3, Rr * 1e3, 4, 2.0, 0.8, 35.0, 5.0, name="pm4")
    qmags = [c for c in quad.children if c.label != "rotor_iron"]
    assert len({round(parse_magnetization(m.label)) % 360 for m in qmags}) >= 3, "4-pole N/S alternation"
    stator = slotted_stator(Rm * 1e3, Rs * 1e3 + 6, 12, 4.0, 12.0, 5.0, name="stator")
    assert stator.is_valid and stator.volume > 0

    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    iron = [c for c in rotor.children if c.label == "rotor_iron"][0]
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "iron.step")
        export_step(iron, f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=2.0))
    assert mesh.ne > 50, "the IPM rotor iron tet-meshes through STEP -> Netgen"


# ---- SOLVE half: AGE synchronous torque of the same machine ----------------------------------------
def _geo():
    from netgen.geom2d import SplineGeometry
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=2.0 * mm)
    g.AddCircle((0, 0), Rw, leftdomain=5, rightdomain=1, bc="wind_outer", maxh=0.4 * mm)
    g.AddCircle((0, 0), Rs, leftdomain=0, rightdomain=5, bc="stator_ring", maxh=0.12 * mm)
    g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=0, bc="rotor_ring", maxh=0.12 * mm)
    g.AddCircle((0, 0), Rr, leftdomain=2, rightdomain=3, bc="mag_inner", maxh=0.5 * mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=2, bc="rotor_inner", maxh=0.6 * mm)
    g.SetMaterial(1, "stator"); g.SetMaterial(2, "rotoriron")
    g.SetMaterial(3, "magnet"); g.SetMaterial(5, "winding")
    from ngsolve import Mesh
    return Mesh(g.GenerateMesh(maxh=1.0 * mm))


def test_ipm_age_synchronous_torque():
    from ngsolve import (H1, BilinearForm, LinearForm, CoefficientFunction, grad, dx, x, y,
                         atan2, cos, sin, IfPos, TaskManager)
    from radia_mcp.radia_ngsolve.airgap_machine import (airgap_coupling, airgap_factorize, airgap_solve,
                                                        airgap_torque)
    mesh = _geo()
    nu = mesh.MaterialCF({"stator": 1.0 / MUR_FE, "rotoriron": 1.0 / MUR_FE,
                          "magnet": 1.0 / MUR_PM, "winding": 1.0}, default=1.0)
    fes = H1(mesh, order=3, dirichlet="outer|rotor_inner", complex=True)
    u, w = fes.TnT()
    a = BilinearForm(fes); a += nu * grad(u) * grad(w) * dx

    def source(theta_m, delta):
        th = atan2(y, x)
        chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
        sgn = IfPos(cos(th - theta_m), 1.0, -1.0)             # 2-pole magnet pattern (p=1)
        chiw = mesh.MaterialCF({"winding": 1.0}, default=0.0)
        phi_s = theta_m + delta - math.pi / 2
        L = LinearForm(fes)
        L += chim * (1.0 / MUR_PM) * (BR * sgn * cos(th) * grad(w)[1] - BR * sgn * sin(th) * grad(w)[0]) * dx
        L += chiw * MU0 * (J0 * cos(th - phi_s)) * w * dx
        L.Assemble(); return L

    with TaskManager():
        a.Assemble()
        coup = airgap_coupling(fes, Rm, Rs, "rotor_ring", "stator_ring", HARM)
        fac = airgap_factorize(a.mat, coup, fes.FreeDofs())

        def torque(tm, de):
            g = airgap_solve(fac, fes, dirichlet_cf=None, source_lf=source(tm, de))
            return airgap_torque(coup, g, axial_length=LSTK).real

        thetas = np.linspace(0, 2 * math.pi, 7)[:-1]          # phase-locked at quadrature
        Tlock = np.array([torque(tm, math.pi / 2) for tm in thetas])
        deltas = np.linspace(0, 2 * math.pi, 13)[:-1]         # torque-angle law at fixed rotor
        Tdelta = np.array([torque(0.0, de) for de in deltas])

    ripple = (Tlock.max() - Tlock.min()) / abs(Tlock.mean())
    s = np.sin(deltas)
    Tmax = np.sum(Tdelta * s) / np.sum(s * s)
    resid = np.sqrt(np.mean((Tdelta - Tmax * s) ** 2)) / np.max(np.abs(Tdelta))
    print(f"  AGE 2-pole IPM: |T_lock|={abs(Tlock.mean())*1e3:.2f} mN.m, ripple={ripple:.2e}; "
          f"T_max={abs(Tmax)*1e3:.2f} mN.m, single-sin resid={resid:.2e}")
    assert abs(Tlock.mean()) > 1e-3, "synchronous torque should be physical"
    assert ripple < 1e-2, f"phase-locked torque must be steady (ripple {ripple:.2e})"
    assert resid < 5e-3, f"torque-angle should be a clean single sinusoid (resid {resid:.2e})"


def main():
    test_ipm_rotor_and_stator_cad()
    test_ipm_age_synchronous_torque()
    print("[OK] end-to-end: ipm_rotor / slotted_stator (promoted CAD archetypes, meshable, labelled) + "
          "AGE synchronous torque of the same machine (phase-locked constant torque, T=T_max*sin(delta)).")


if __name__ == "__main__":
    main()
