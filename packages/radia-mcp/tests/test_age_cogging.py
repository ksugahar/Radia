# -*- coding: utf-8 -*-
r"""Open-circuit COGGING torque of a slotted SPM via the AGE, vs the brute meshed-gap FE.

Cogging = the PM-only (zero stator current) reluctance torque from the stator slot permeance. The
slotting is a smooth angular reluctivity modulation nu~(theta) in an annular stator zone (tooth = iron,
slot = air, with a sigmoid tooth/slot transition ~ the physical slot-opening flux fringing / Carter
effect), so the mesh stays smooth and AGE & brute share the identical nu~ field. The rotor PM is rotated
through one cogging period 360/LCM(slots,poles); the AGE keeps the gap analytic (no remesh = no mesh
noise -- the analytic version of the classic "freeze the mesh, move the material" cogging recipe), with
the mesh-free closed-form torque. Validated against the brute meshed-gap Arkkio torque.

Note (honest): cogging needs MANY gap harmonics (the slot permeance has slow ~1/n harmonics), unlike the
smooth-gap loaded torque which needs only n=1. A physical (smoothed) slot opening converges with ~3*Q
harmonics; an idealised sharp slot needs more and shows Gibbs ringing. 2-pole / Q=6 -> period 60 deg.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction,
                     grad, dx, x, y, sqrt, atan2, cos, sin, exp, IfPos, Integrate, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import (
    airgap_coupling, airgap_factorize, airgap_solve, airgap_torque)

mm = 1e-3
R0, Rr, Rm, Rs, Rslot, REXT = 3*mm, 8*mm, 11*mm, 11.5*mm, 13.5*mm, 20*mm
MUR_FE, MUR_PM, BR = 1000.0, 1.05, 1.2
MU0 = 4e-7*math.pi
LSTK = 0.05
NSLOTS, POLES = 6, 2
HARM = list(range(1, 21))
PERIOD = 360.0/np.lcm(NSLOTS, POLES)        # 60 deg


def _geo(meshed_gap):
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=1.5*mm)
    g.AddCircle((0, 0), Rslot, leftdomain=6, rightdomain=1, bc="slot_outer", maxh=0.5*mm)
    if meshed_gap:
        g.AddCircle((0, 0), Rs, leftdomain=4, rightdomain=6, bc="stator_ring", maxh=0.1*mm)
        g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=4, bc="rotor_ring", maxh=0.1*mm)
    else:
        g.AddCircle((0, 0), Rs, leftdomain=0, rightdomain=6, bc="stator_ring", maxh=0.1*mm)
        g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=0, bc="rotor_ring", maxh=0.1*mm)
    g.AddCircle((0, 0), Rr, leftdomain=2, rightdomain=3, bc="mag_inner", maxh=0.4*mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=2, bc="rotor_inner", maxh=0.6*mm)
    g.SetMaterial(1, "stator_back"); g.SetMaterial(2, "rotor_iron")
    g.SetMaterial(3, "magnet"); g.SetMaterial(6, "slotzone")
    if meshed_gap:
        g.SetMaterial(4, "gap")
    return Mesh(g.GenerateMesh(maxh=1.0*mm))


def _nu(mesh):
    th = atan2(y, x)
    base = mesh.MaterialCF({"stator_back": 1/MUR_FE, "rotor_iron": 1/MUR_FE,
                            "magnet": 1/MUR_PM, "gap": 1.0, "slotzone": 0.0}, default=1.0)
    chi = mesh.MaterialCF({"slotzone": 1.0}, default=0.0)
    t = 1.0/(1.0 + exp(-3.0*cos(NSLOTS*th)))            # smooth tooth(1)/slot(0)
    return base + chi*((1.0/MUR_FE)*t + 1.0*(1.0 - t))


def _pm(fes, mesh, theta_r):
    w = fes.TestFunction(); th = atan2(y, x)
    chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
    sgn = IfPos(cos(th - theta_r), 1.0, -1.0)
    L = LinearForm(fes)
    L += chim*(1.0/MUR_PM)*(BR*sgn*cos(th)*grad(w)[1] - BR*sgn*sin(th)*grad(w)[0])*dx
    L.Assemble(); return L


def test_cogging_age_matches_brute():
    ma = _geo(False)
    fa = H1(ma, order=3, dirichlet="outer|rotor_inner", complex=True)
    ua, wa = fa.TnT(); aa = BilinearForm(fa); aa += _nu(ma)*grad(ua)*grad(wa)*dx
    mb = _geo(True)
    fb = H1(mb, order=3, dirichlet="outer|rotor_inner", complex=False)
    ub, wb = fb.TnT(); ab = BilinearForm(fb); ab += _nu(mb)*grad(ub)*grad(wb)*dx
    with TaskManager():
        aa.Assemble(); ab.Assemble()
        coup = airgap_coupling(fa, Rm, Rs, "rotor_ring", "stator_ring", HARM)
        fac = airgap_factorize(aa.mat, coup, fa.FreeDofs())
        Kb = ab.mat.Inverse(fb.FreeDofs(), inverse="umfpack")
        r_cf = sqrt(x*x+y*y); gapdx = dx(definedon=mb.Materials("gap"))
        angs = np.linspace(0, PERIOD, 7)
        Ta, Tb = [], []
        for th in angs:
            tr = math.radians(th)
            ga = airgap_solve(fac, fa, dirichlet_cf=None, source_lf=_pm(fa, ma, tr))
            Ta.append(airgap_torque(coup, ga, axial_length=LSTK).real)
            gb = GridFunction(fb); gb.vec.data = Kb*_pm(fb, mb, tr).vec
            B = CoefficientFunction((grad(gb)[1], -grad(gb)[0]))
            Br_ = (B[0]*x+B[1]*y)/r_cf; Bth = (-B[0]*y+B[1]*x)/r_cf
            Tb.append((LSTK/(MU0*(Rs-Rm)))*Integrate(r_cf*Br_*Bth*gapdx, mb))
    Ta, Tb = np.array(Ta), np.array(Tb)
    pk = np.max(np.abs(Tb))
    rel = np.max(np.abs(Ta-Tb))/pk
    print(f"cogging: brute peak={pk*1e3:.2f} mN.m, AGE peak={np.max(np.abs(Ta))*1e3:.2f}; "
          f"max|AGE-brute|={rel:.2e} of peak; endpoints brute={Tb[0]*1e3:.2f}/{Tb[-1]*1e3:.2f}")
    assert pk > 1e-3, "cogging should be a physical (mN.m) reluctance torque, not ~0"
    assert abs(Tb[0]) < 0.1*pk and abs(Tb[-1]) < 0.1*pk, "cogging ~0 at the magnet-aligned period ends"
    assert rel < 0.10, f"AGE cogging off the brute meshed-gap by {rel:.2e} of peak"


def main():
    test_cogging_age_matches_brute()
    print("[OK] AGE cogging: PM-only slotted-SPM reluctance torque over one 360/LCM period, un-meshed "
          "gap (no remesh), reproduces the brute meshed-gap Arkkio cogging curve.")


if __name__ == "__main__":
    main()
