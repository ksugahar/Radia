# -*- coding: utf-8 -*-
r"""SLIDING-BAND FRONTIER: salient rotor AND slotted stator together, via per-angle AGE vs brute.

This locks the ONE case the salient-rotor test (test_age_salient_ipm) flagged as breaking AGE's
one-factorization: BOTH the rotor (salient, nu~ depends on theta - theta_m, rotating) and the stator
(slotted, nu~ depends on theta, fixed) are angle-structured, so NEITHER reference frame freezes the
other -- K(theta_m) genuinely changes with rotor angle and must be RE-ASSEMBLED + RE-FACTORIZED per
angle (the classical moving-mesh / sliding-band regime).

The thesis, locked here: the AGE still solves it CORRECTLY and MESH-FREE in the gap.  The gap rings
(Rm, Rs) stay circular, so airgap_coupling is built ONCE and reused at every angle; only the FE
interior K re-factorizes.  Validated: AGE (analytic gap, per-angle re-factorize) == brute (meshed
gap, per-angle re-invert) at every rotor angle across one slot pitch, and the physical slot-passing
torque RIPPLE is reproduced.

HONEST cost (documented, not a speed win): this doubly-structured regime is harmonic-hungry -- the
slot x saliency interaction injects high gap harmonics (Q+-1 = 5,7,11,13 ...), so the AGE coupling
needs MANY harmonics (here 1..15; [1,3,5] alone is ~25% off), and it loses the one-factorization
(N re-factorizations).  So unlike the smooth-stator machines (#1-#6, one factorization, few
harmonics), here the AGE's only remaining advantage over a meshed sliding band is the mesh-free gap
itself (no thin-gap layers), not solve speed.  Phase-locked current (smooth stator would give the
constant operating torque of #3) isolates the slot ripple on top of the operating torque.
"""
import math
import os
import sys

import pytest
import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction,
                     grad, dx, x, y, sqrt, atan2, cos, sin, exp, Integrate, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

pytestmark = pytest.mark.xval

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import (
    airgap_coupling, airgap_factorize, airgap_solve, airgap_torque)

mm = 1e-3
R0, Rm, Rs, Rslot, REXT = 3*mm, 11*mm, 11.5*mm, 13.5*mm, 20*mm
MUR_FE = 1000.0
MU0 = 4e-7*math.pi
LSTK = 0.05
J0 = 6.0e6
NSLOTS = 6
HARM = list(range(1, 16))      # slot x saliency -> high gap harmonics; [1,3,5] alone ~25% off
DELTA = math.pi/4


def _geo(meshed_gap):
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=1.5*mm)
    g.AddCircle((0, 0), Rslot, leftdomain=6, rightdomain=1, bc="slot_outer", maxh=0.4*mm)
    if meshed_gap:
        g.AddCircle((0, 0), Rs, leftdomain=4, rightdomain=6, bc="stator_ring", maxh=0.1*mm)
        g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=4, bc="rotor_ring", maxh=0.1*mm)
    else:
        g.AddCircle((0, 0), Rs, leftdomain=0, rightdomain=6, bc="stator_ring", maxh=0.1*mm)
        g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=0, bc="rotor_ring", maxh=0.1*mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=3, bc="rotor_inner", maxh=0.6*mm)
    g.SetMaterial(1, "stator"); g.SetMaterial(3, "rotor"); g.SetMaterial(6, "slotzone")
    if meshed_gap:
        g.SetMaterial(4, "gap")
    return Mesh(g.GenerateMesh(maxh=1.0*mm))


def _nu(mesh, theta_m):
    th = atan2(y, x)
    t_slot = 1.0/(1.0 + exp(-3.0*cos(NSLOTS*th)))
    nu_slot = (1.0/MUR_FE)*t_slot + 1.0*(1.0 - t_slot)
    sal = (1.0/MUR_FE) + (1.0 - 1.0/MUR_FE)*sin(th - theta_m)**2
    chi_slot = mesh.MaterialCF({"slotzone": 1.0}, default=0.0)
    chi_rot = mesh.MaterialCF({"rotor": 1.0}, default=0.0)
    other = mesh.MaterialCF({"stator": 1.0/MUR_FE, "gap": 1.0}, default=1.0)
    return chi_slot*nu_slot + chi_rot*sal + (1.0 - chi_slot - chi_rot)*other


def _cur(fes, mesh, theta_m):
    w = fes.TestFunction(); th = atan2(y, x)
    chiw = mesh.MaterialCF({"slotzone": 1.0}, default=0.0)
    L = LinearForm(fes)
    L += chiw*MU0*(J0*cos(th - (theta_m + DELTA)))*w*dx
    L.Assemble(); return L


def test_slotted_salient_age_matches_brute():
    ma = _geo(False)
    fa = H1(ma, order=3, dirichlet="outer|rotor_inner", complex=True)
    ua, wa = fa.TnT()
    mb = _geo(True)
    fb = H1(mb, order=3, dirichlet="outer|rotor_inner", complex=False)
    ub, wb = fb.TnT()
    r_cf = sqrt(x*x+y*y)
    with TaskManager():
        coup = airgap_coupling(fa, Rm, Rs, "rotor_ring", "stator_ring", HARM)   # ONCE (gap fixed)
        gapdx = dx(definedon=mb.Materials("gap"))
        thetas = np.linspace(0, 2*math.pi/NSLOTS, 7)               # one slot pitch
        Ta, Tb = [], []
        for tm in thetas:
            aa = BilinearForm(fa); aa += _nu(ma, tm)*grad(ua)*grad(wa)*dx; aa.Assemble()
            fac = airgap_factorize(aa.mat, coup, fa.FreeDofs())     # re-factorize per angle
            ga = airgap_solve(fac, fa, source_lf=_cur(fa, ma, tm))
            Ta.append(airgap_torque(coup, ga, axial_length=LSTK).real)
            ab = BilinearForm(fb); ab += _nu(mb, tm)*grad(ub)*grad(wb)*dx; ab.Assemble()
            Kb = ab.mat.Inverse(fb.FreeDofs(), inverse="umfpack")   # re-invert per angle
            gb = GridFunction(fb); gb.vec.data = Kb*_cur(fb, mb, tm).vec
            B = CoefficientFunction((grad(gb)[1], -grad(gb)[0]))
            Br_ = (B[0]*x+B[1]*y)/r_cf; Bth = (-B[0]*y+B[1]*x)/r_cf
            Tb.append((LSTK/(MU0*(Rs-Rm)))*Integrate(r_cf*Br_*Bth*gapdx, mb))
    Ta, Tb = np.array(Ta)*1e3, np.array(Tb)*1e3
    pk = np.max(np.abs(Tb))
    relmax = np.max(np.abs(Ta-Tb))/pk
    dc = np.mean(Tb); ripple = (Tb.max()-Tb.min())/abs(dc)
    period = abs(Tb[0]-Tb[-1])/pk
    print(f"slotted+salient (per-angle, {len(HARM)} harmonics): base={dc:.3f} mN.m, "
          f"slot ripple={ripple*100:.1f}% pk-pk, period close={period:.2e}; "
          f"AGE-vs-brute per-angle max={relmax:.2e}")

    assert abs(dc) > 1e-3, "operating torque should be physical"
    assert ripple > 0.3, f"slotted stator must add a real slot ripple, got {ripple*100:.1f}%"
    assert period < 1e-2, f"torque must repeat over one slot pitch, endpoint diff {period:.2e}"
    assert relmax < 1e-2, f"per-angle AGE (analytic gap) off brute (meshed gap) by {relmax:.2e}"


def main():
    test_slotted_salient_age_matches_brute()
    print("[OK] AGE sliding-band frontier: salient rotor + slotted stator, per-angle re-factorization "
          "(gap coupling reused), reproduces the brute meshed-gap torque + slot ripple mesh-free.")


if __name__ == "__main__":
    main()
