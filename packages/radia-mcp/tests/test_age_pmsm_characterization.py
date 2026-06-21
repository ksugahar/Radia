# -*- coding: utf-8 -*-
r"""PMSM characterisation suite via the AGE -- dq machine parameters straight from the FE field.

Extracts the dq parameters from the AGE FE solution (the bridge FIELD -> CIRCUIT) and locks the
textbook energy consistency, each gate cross-checked where possible AGE (mesh-free) vs brute:

 (1) BACK-EMF / flux linkage + THD: open-circuit (PM only), rotate the rotor; the phase flux
     linkage of a sinusoidally-distributed winding is  Lam(theta_m) = INT_winding cos(theta)*A_z dA
     (here ~ Lam1*sin(theta_m), peak at quadrature).  FFT -> fundamental Lam1 and THD.  A sinusoidal
     winding filters the PM square-wave / slot harmonics, so the back-EMF THD is small.

 (2) Kt == Ke ENERGY GATE (parameter-free): by reciprocity the winding turns-constant C cancels, so
     the loaded quadrature torque must equal  T = p*L*J0*Lam1  (p=1).  This links the open-circuit
     flux linkage (back-EMF constant Ke) to the loaded torque (torque constant Kt) with NO turns
     count -- a strong FE energy-consistency check.

 (3) Ld, Lq, dq decoupling (PM off): inject d- and q-axis current; the self flux-linkage ratio
     Ld/Lq = lam_dd/lam_qq (C cancels).  The uniform-ring surface PM has NO reluctance saliency
     (Ld/Lq = 1) and the cross term lam_qd ~ 0 (dq decoupled).  A salient IPM/SynRM breaks both
     (see the salient-rotor test) -- this is the non-salient baseline.

The current sheet J0*cos(theta-phi_s) has its stator-field d-axis at phi_s+90deg; d-axis current
=> phi_s=-90deg (weight -sin theta), q-axis current => phi_s=0 (weight cos theta).  nu~=1/mur.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, dx, x, y,
                     atan2, cos, sin, IfPos, Integrate, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

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
P = 1
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


def _pm(fes, mesh, theta_m):
    w = fes.TestFunction(); th = atan2(y, x)
    chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
    sgn = IfPos(cos(th - theta_m), 1.0, -1.0)
    L = LinearForm(fes)
    L += chim*(1.0/MUR_PM)*(BR*sgn*cos(th)*grad(w)[1] - BR*sgn*sin(th)*grad(w)[0])*dx
    L.Assemble(); return L


def _cur(fes, mesh, phi_s, theta_m=None):
    w = fes.TestFunction(); th = atan2(y, x)
    chiw = mesh.MaterialCF({"winding": 1.0}, default=0.0)
    L = LinearForm(fes)
    L += chiw*MU0*(J0*cos(th - phi_s))*w*dx
    if theta_m is not None:                                   # add PM (loaded case)
        chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
        sgn = IfPos(cos(th - theta_m), 1.0, -1.0)
        L += chim*(1.0/MUR_PM)*(BR*sgn*cos(th)*grad(w)[1] - BR*sgn*sin(th)*grad(w)[0])*dx
    L.Assemble(); return L


def test_pmsm_characterization_age():
    ma = _geo(False)
    fa = H1(ma, order=3, dirichlet="outer|rotor_inner", complex=True)
    ua, wa = fa.TnT(); aa = BilinearForm(fa); aa += _nu(ma)*grad(ua)*grad(wa)*dx
    mb = _geo(True)
    fb = H1(mb, order=3, dirichlet="outer|rotor_inner", complex=False)
    ub, wb = fb.TnT(); ab = BilinearForm(fb); ab += _nu(mb)*grad(ub)*grad(wb)*dx
    th = atan2(y, x)
    win_a = ma.MaterialCF({"winding": 1.0}, default=0.0)
    win_b = mb.MaterialCF({"winding": 1.0}, default=0.0)
    with TaskManager():
        aa.Assemble(); ab.Assemble()
        coup = airgap_coupling(fa, Rm, Rs, "rotor_ring", "stator_ring", HARM)
        fac = airgap_factorize(aa.mat, coup, fa.FreeDofs())
        Kb = ab.mat.Inverse(fb.FreeDofs(), inverse="umfpack")

        def lk_age(gfu, weight):
            return Integrate(win_a*weight*gfu*dx, ma).real

        # (1) back-EMF flux linkage sweep
        NTH = 24
        thetas = np.linspace(0, 2*math.pi, NTH, endpoint=False)
        Lam = np.array([lk_age(airgap_solve(fac, fa, source_lf=_pm(fa, ma, tm)), cos(th))
                        for tm in thetas])
        F = np.fft.rfft(Lam)/NTH
        Lam1 = 2*abs(F[1])
        thd = math.sqrt(np.sum(np.abs(2*F[2:])**2))/Lam1

        # lambda_pm AGE vs brute at the quadrature peak (theta_m=90deg)
        ga90 = airgap_solve(fac, fa, source_lf=_pm(fa, ma, math.pi/2))
        gb90 = GridFunction(fb); gb90.vec.data = Kb*_pm(fb, mb, math.pi/2).vec
        lam_a = lk_age(ga90, cos(th))
        lam_b = Integrate(win_b*cos(th)*gb90*dx, mb)
        rel_pm = abs(lam_a-lam_b)/abs(lam_b)

        # (2) Kt==Ke energy gate
        gL = airgap_solve(fac, fa, source_lf=_cur(fa, ma, math.pi/2 - math.pi/2, theta_m=0.0))
        T_load = abs(airgap_torque(coup, gL, axial_length=LSTK).real)
        T_pred = P*LSTK*J0*Lam1
        rel_KtKe = abs(T_load-T_pred)/T_pred

        # (3) Ld, Lq, cross-coupling (PM off)
        gd = airgap_solve(fac, fa, source_lf=_cur(fa, ma, -math.pi/2))
        lam_dd = lk_age(gd, -sin(th)); lam_qd = lk_age(gd, cos(th))
        gq = airgap_solve(fac, fa, source_lf=_cur(fa, ma, 0.0))
        lam_qq = lk_age(gq, cos(th)); lam_dq = lk_age(gq, -sin(th))

    saliency = lam_dd/lam_qq
    cross = max(abs(lam_qd/lam_dd), abs(lam_dq/lam_qq))
    print(f"(1) back-EMF: Lam1={Lam1:.4e} Wb.m, THD={thd:.2e}; lambda_pm AGE-vs-brute={rel_pm:.2e}")
    print(f"(2) Kt==Ke: T_load={T_load*1e3:.4f} vs p*L*J0*Lam1={T_pred*1e3:.4f} mN.m (rel {rel_KtKe:.2e})")
    print(f"(3) Ld/Lq={saliency:.4f} (non-salient SPM ~1), dq cross-coupling={cross:.2e}")

    assert Lam1 > 1e-9, "PM flux-linkage fundamental should be physical"
    assert thd < 1e-2, f"sinusoidal-winding back-EMF THD should be small, got {thd:.2e}"
    assert rel_pm < 5e-3, f"lambda_pm AGE vs brute off by {rel_pm:.2e}"
    assert rel_KtKe < 5e-3, f"Kt!=Ke energy inconsistency {rel_KtKe:.2e} (field<->torque)"
    assert 0.9 < saliency < 1.1, f"uniform-ring SPM must be non-salient, Ld/Lq={saliency:.3f}"
    assert cross < 1e-3, f"dq must be decoupled for the symmetric SPM, cross={cross:.2e}"


def main():
    test_pmsm_characterization_age()
    print("[OK] AGE PMSM characterisation: back-EMF (low-THD sinusoidal), Kt==Ke energy consistency "
          "(field<->circuit, parameter-free), Ld=Lq non-salient baseline + dq decoupling.")


if __name__ == "__main__":
    main()
