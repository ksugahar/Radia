# -*- coding: utf-8 -*-
r"""SYNCHRONOUS operating-point torque of an SPM via the AGE, vs the brute meshed-gap FE.

Extends the loaded-torque core (test_age_pmsm_physical) to a controlled SYNCHRONOUS operating
point: the stator current vector is PHASE-LOCKED to the rotor (electrical angle theta_e = p*theta_m;
2-pole p=1) at a commanded load angle delta = angle(stator-field d-axis, rotor d-axis).  Two
defining signatures of synchronous running, each gated and cross-checked AGE (mesh-free) vs brute:

  GATE A -- phase-lock gives CONSTANT torque: at the max-torque quadrature delta=90deg, advancing
    the rotor a full electrical turn keeps the torque steady (the stator field co-rotates with the
    rotor; the n=1 interaction is rotation-invariant).  Ripple ~ 0.

  GATE B -- torque-angle law: at fixed rotor, sweeping delta gives the textbook cross-field law
    T(delta) = T_max * sin(delta) -- a single odd sinusoid through 0, max at quadrature, with NO DC
    and NO 2*delta reluctance term (this is the non-salient Ld=Lq baseline; a salient IPM/SynRM
    would add a sin(2*delta) reluctance torque -- see the salient-rotor test).

The current sheet J=J0*cos(th - phi_s) makes its stator-field d-axis at phi_s+90deg, so
phi_s = theta_m + delta - 90deg gives the standard law above.  nu~ = 1/mur (gap nu~=1); one AGE
factorization (and one brute K-inverse) serves the whole sweep.
"""
import math
import os
import sys

import pytest
import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction,
                     grad, dx, x, y, sqrt, atan2, cos, sin, IfPos, Integrate, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

pytestmark = pytest.mark.xval

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import (
    airgap_coupling, airgap_factorize, airgap_solve, airgap_torque)

mm = 1e-3
R0, Rr, Rm, Rs, Rw, REXT = 3*mm, 8*mm, 11*mm, 11.5*mm, 12.5*mm, 22*mm
MUR_FE, MUR_PM, BR = 5.0e4, 1.05, 1.2
MU0 = 4e-7*math.pi
LSTK = 0.05
J0 = 6.0e6
HARM = [1, 3, 5]


def _geo(meshed_gap):
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=2.0*mm)
    g.AddCircle((0, 0), Rw, leftdomain=5, rightdomain=1, bc="wind_outer", maxh=0.4*mm)
    if meshed_gap:
        g.AddCircle((0, 0), Rs, leftdomain=4, rightdomain=5, bc="stator_ring", maxh=0.12*mm)
        g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=4, bc="rotor_ring", maxh=0.12*mm)
    else:
        g.AddCircle((0, 0), Rs, leftdomain=0, rightdomain=5, bc="stator_ring", maxh=0.12*mm)
        g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=0, bc="rotor_ring", maxh=0.12*mm)
    g.AddCircle((0, 0), Rr, leftdomain=2, rightdomain=3, bc="mag_inner", maxh=0.5*mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=2, bc="rotor_inner", maxh=0.6*mm)
    g.SetMaterial(1, "stator"); g.SetMaterial(2, "rotoriron")
    g.SetMaterial(3, "magnet"); g.SetMaterial(5, "winding")
    if meshed_gap:
        g.SetMaterial(4, "gap")
    return Mesh(g.GenerateMesh(maxh=1.0*mm))


def _nu(mesh):
    return mesh.MaterialCF({"stator": 1.0/MUR_FE, "rotoriron": 1.0/MUR_FE,
                            "magnet": 1.0/MUR_PM, "gap": 1.0, "winding": 1.0}, default=1.0)


def _source(fes, mesh, theta_m, delta):
    w = fes.TestFunction(); th = atan2(y, x)
    chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
    sgn = IfPos(cos(th - theta_m), 1.0, -1.0)
    brx, bry = BR*sgn*cos(th), BR*sgn*sin(th)
    chiw = mesh.MaterialCF({"winding": 1.0}, default=0.0)
    phi_s = theta_m + delta - math.pi/2
    L = LinearForm(fes)
    L += chim*(1.0/MUR_PM)*(brx*grad(w)[1] - bry*grad(w)[0])*dx
    L += chiw*MU0*(J0*cos(th - phi_s))*w*dx
    L.Assemble(); return L


def test_synchronous_operating_point_age_matches_brute():
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

        def t_age(tm, de):
            ga = airgap_solve(fac, fa, dirichlet_cf=None, source_lf=_source(fa, ma, tm, de))
            return airgap_torque(coup, ga, axial_length=LSTK).real

        def t_brute(tm, de):
            gb = GridFunction(fb); gb.vec.data = Kb*_source(fb, mb, tm, de).vec
            B = CoefficientFunction((grad(gb)[1], -grad(gb)[0]))
            Br_ = (B[0]*x+B[1]*y)/r_cf; Bth = (-B[0]*y+B[1]*x)/r_cf
            return (LSTK/(MU0*(Rs-Rm)))*Integrate(r_cf*Br_*Bth*gapdx, mb)

        # GATE A: phase-locked at quadrature (delta=90deg) -> torque constant over a full turn.
        thetas = np.linspace(0, 2*math.pi, 9)[:-1]
        TaA = np.array([t_age(tm, math.pi/2) for tm in thetas])
        TbA = np.array([t_brute(tm, math.pi/2) for tm in thetas])
        # GATE B: torque-angle law at fixed rotor.
        deltas = np.linspace(0, 2*math.pi, 17)[:-1]
        TaB = np.array([t_age(0.0, de) for de in deltas])
        TbB = np.array([t_brute(0.0, de) for de in deltas])

    ripple_b = (TbA.max()-TbA.min())/abs(TbA.mean())
    ripple_a = (TaA.max()-TaA.min())/abs(TaA.mean())
    relA = np.max(np.abs(TaA-TbA))/abs(TbA.mean())
    s = np.sin(deltas)
    Tmax = np.sum(TbB*s)/np.sum(s*s)
    Tpk = np.max(np.abs(TbB))
    dc = abs(np.mean(TbB))/Tpk
    rel2 = abs(np.sum(TbB*np.sin(2*deltas))/np.sum(np.sin(2*deltas)**2))/Tpk
    resid = np.sqrt(np.mean((TbB - Tmax*s)**2))/Tpk
    relB = np.max(np.abs(TaB-TbB))/Tpk
    print(f"GATE A (phase-lock): |T|={abs(TbA.mean())*1e3:.2f} mN.m steady; ripple brute={ripple_b:.2e} "
          f"AGE={ripple_a:.2e}; AGE-vs-brute={relA:.2e}")
    print(f"GATE B (T-angle): T_max={abs(Tmax)*1e3:.2f} mN.m; single-sin resid={resid:.2e}; "
          f"DC={dc:.2e}, 2delta(reluctance)={rel2:.2e} of peak; AGE-vs-brute={relB:.2e}")

    assert abs(TbA.mean()) > 1e-3, "synchronous torque should be physical (mN.m)"
    assert ripple_b < 1e-2 and ripple_a < 1e-2, "phase-locked torque must be steady (ripple ~ 0)"
    assert relA < 5e-3, f"AGE phase-locked torque off brute by {relA:.2e}"
    assert resid < 5e-3, f"torque-angle not a clean single sinusoid (resid {resid:.2e})"
    assert dc < 1e-2, f"torque-angle has spurious DC offset {dc:.2e}"
    assert rel2 < 1e-2, f"non-salient SPM must have ~0 reluctance (2delta) torque, got {rel2:.2e}"
    assert relB < 5e-3, f"AGE torque-angle off brute by {relB:.2e}"


def main():
    test_synchronous_operating_point_age_matches_brute()
    print("[OK] AGE synchronous operating point: phase-locked stator current -> constant torque, "
          "T(delta)=T_max*sin(delta) cross-field law (no reluctance term), AGE == brute meshed gap.")


if __name__ == "__main__":
    main()
