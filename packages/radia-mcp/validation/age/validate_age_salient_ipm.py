# -*- coding: utf-8 -*-
r"""SALIENT rotor (SynRM / IPM) via the AGE -- reluctance torque + dq saliency, vs brute meshed gap.

The non-salient SPM (validate_age_synchronous_torque / validate_age_pmsm_characterization) had Ld=Lq and a
pure sin(delta) torque with NO reluctance (2*delta) term.  Here the rotor is SALIENT -- a 2-pole
anisotropic core, low reluctivity (iron) on the d-axis, high (air barrier) on the q-axis, as a smooth
rotor-fixed modulation  nu~_rotor(theta) = 1/mur_fe + (1-1/mur_fe)*sin^2(theta).

The AGE stays ONE-FACTORIZATION for this interior-salient rotor: the rotor SURFACE (gap ring at Rm) is
still a circle, so airgap_coupling is unchanged; and because the STATOR is smooth we work in the ROTOR
FRAME (saliency fixed in the mesh = K assembled once) and rotate the stator CURRENT (excitation phase).
By linearity A^delta = sin(delta) A^q + cos(delta) A^d, and a pure-axis current on a salient rotor
makes no torque (T_qq=T_dd=0), so the reluctance torque is the cross term  T = sin(2 delta)*B(A^q,A^d)
-- a PURE 2nd harmonic, max at 45deg.  Adding a d-axis PM (IPM) adds 2*B(A^delta,A^PM) = a pure
FUNDAMENTAL (PM torque), independent of the reluctance 2nd harmonic.

Gates (phase-independent full-period magnitudes, no axis convention): (1) saliency Ld/Lq>1; (2) SynRM
torque is a pure 2nd harmonic (fundamental~0, clean sin2d) and AGE==brute; (3) the IPM PM switches the
fundamental on while leaving the reluctance 2nd harmonic unchanged (independent superposition).

HONEST limitation: one-factorization holds because the stator is smooth.  A salient rotor AND a slotted
stator simultaneously (both angle-structured) is the case needing a true per-angle re-assembly / sliding
band -- neither frame freezes the other; that is the frontier beyond this test.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction,
                     grad, dx, x, y, sqrt, atan2, cos, sin, Integrate, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import (
    airgap_coupling, airgap_factorize, airgap_solve, airgap_torque)

mm = 1e-3
R0, Rpm, Rm, Rs, Rw, REXT = 3*mm, 7*mm, 11*mm, 11.5*mm, 12.5*mm, 22*mm
MUR_FE, MUR_PM, BR = 1000.0, 1.05, 0.6
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
    g.AddCircle((0, 0), Rpm, leftdomain=2, rightdomain=3, bc="pm_inner", maxh=0.4*mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=2, bc="rotor_inner", maxh=0.6*mm)
    g.SetMaterial(1, "stator"); g.SetMaterial(3, "rotor")
    g.SetMaterial(2, "pmcore"); g.SetMaterial(5, "winding")
    if meshed_gap:
        g.SetMaterial(4, "gap")
    return Mesh(g.GenerateMesh(maxh=1.0*mm))


def _nu(mesh):
    th = atan2(y, x)
    sal = (1.0/MUR_FE) + (1.0 - 1.0/MUR_FE)*sin(th)**2          # d-axis iron, q-axis air barrier
    chi_rot = mesh.MaterialCF({"rotor": 1.0, "pmcore": 1.0}, default=0.0)
    nonrot = mesh.MaterialCF({"stator": 1.0/MUR_FE, "winding": 1.0, "gap": 1.0}, default=1.0)
    return chi_rot*sal + (1.0 - chi_rot)*nonrot


def _cur(fes, mesh, delta, with_pm=False):
    w = fes.TestFunction(); th = atan2(y, x)
    chiw = mesh.MaterialCF({"winding": 1.0}, default=0.0)
    phi_s = delta - math.pi/2
    L = LinearForm(fes)
    L += chiw*MU0*(J0*cos(th - phi_s))*w*dx
    if with_pm:
        chipm = mesh.MaterialCF({"pmcore": 1.0}, default=0.0)
        L += chipm*(1.0/MUR_PM)*(-BR*grad(w)[0])*dx            # M = Br*xhat (d-axis remanence)
    L.Assemble(); return L


def _fourier_mag(T, d):
    f1 = math.hypot(2*np.mean(T*np.cos(d)), 2*np.mean(T*np.sin(d)))
    f2 = math.hypot(2*np.mean(T*np.cos(2*d)), 2*np.mean(T*np.sin(2*d)))
    return f1, f2


def validate_salient_ipm_synrm_age_matches_brute():
    ma = _geo(False)
    fa = H1(ma, order=3, dirichlet="outer|rotor_inner", complex=True)
    ua, wa = fa.TnT(); aa = BilinearForm(fa); aa += _nu(ma)*grad(ua)*grad(wa)*dx
    mb = _geo(True)
    fb = H1(mb, order=3, dirichlet="outer|rotor_inner", complex=False)
    ub, wb = fb.TnT(); ab = BilinearForm(fb); ab += _nu(mb)*grad(ub)*grad(wb)*dx
    th = atan2(y, x)
    win = ma.MaterialCF({"winding": 1.0}, default=0.0)
    r_cf = sqrt(x*x+y*y)
    with TaskManager():
        aa.Assemble(); ab.Assemble()
        coup = airgap_coupling(fa, Rm, Rs, "rotor_ring", "stator_ring", HARM)
        fac = airgap_factorize(aa.mat, coup, fa.FreeDofs())
        Kb = ab.mat.Inverse(fb.FreeDofs(), inverse="umfpack")
        gapdx = dx(definedon=mb.Materials("gap"))

        def t_age(d, pm=False):
            return airgap_torque(coup, airgap_solve(fac, fa, source_lf=_cur(fa, ma, d, pm)),
                                 axial_length=LSTK).real

        def t_brute(d, pm=False):
            gb = GridFunction(fb); gb.vec.data = Kb*_cur(fb, mb, d, pm).vec
            B = CoefficientFunction((grad(gb)[1], -grad(gb)[0]))
            Br_ = (B[0]*x+B[1]*y)/r_cf; Bth = (-B[0]*y+B[1]*x)/r_cf
            return (LSTK/(MU0*(Rs-Rm)))*Integrate(r_cf*Br_*Bth*gapdx, mb)

        # (1) saliency: two orthogonal axis inductances (PM off)
        lam_a = Integrate(win*(-sin(th))*airgap_solve(fac, fa, source_lf=_cur(fa, ma, 0.0))*dx, ma).real
        lam_b = Integrate(win*cos(th)*airgap_solve(fac, fa, source_lf=_cur(fa, ma, math.pi/2))*dx, ma).real
        saliency = max(lam_a, lam_b)/min(lam_a, lam_b)

        # (2) SynRM reluctance torque (PM off), full period
        deltas = np.linspace(0, 2*math.pi, 24, endpoint=False)
        Ta = np.array([t_age(d)*1e3 for d in deltas])
        Tb = np.array([t_brute(d)*1e3 for d in deltas])
        # (3) IPM (+PM)
        Tap = np.array([t_age(d, True)*1e3 for d in deltas])
        Tbp = np.array([t_brute(d, True)*1e3 for d in deltas])

    f1, f2 = _fourier_mag(Tb, deltas)
    pk = np.max(np.abs(Tb))
    c2, s2 = 2*np.mean(Tb*np.cos(2*deltas)), 2*np.mean(Tb*np.sin(2*deltas))
    resid2 = np.sqrt(np.mean((Tb - (c2*np.cos(2*deltas)+s2*np.sin(2*deltas)))**2))/pk
    relAB = np.max(np.abs(Ta-Tb))/pk
    g1, g2 = _fourier_mag(Tbp, deltas)
    relABp = np.max(np.abs(Tap-Tbp))/np.max(np.abs(Tbp))
    print(f"(1) saliency Ld/Lq = {saliency:.3f}")
    print(f"(2) SynRM: fundamental={f1:.4f} 2nd-harm(reluctance)={f2:.4f} mN.m (fund/2nd={f1/f2:.2e}); "
          f"pure-2nd resid={resid2:.2e}; AGE-vs-brute={relAB:.2e}")
    print(f"(3) IPM: fundamental(PM)={g1:.4f} 2nd-harm(reluctance)={g2:.4f} mN.m; AGE-vs-brute={relABp:.2e}; "
          f"PM switches fundamental {f1:.4f}->{g1:.4f}, reluctance {f2:.4f}->{g2:.4f} (unchanged)")

    # saliency is real
    assert saliency > 1.2, f"rotor must be salient, Ld/Lq={saliency:.3f}"
    # SynRM reluctance torque = pure 2nd harmonic, no fundamental
    assert f2 > 1e-4, "reluctance (2nd-harmonic) torque should be physical"
    assert f1/f2 < 1e-2, f"PM-free SynRM must have no fundamental torque, fund/2nd={f1/f2:.2e}"
    assert resid2 < 5e-3, f"reluctance torque not a pure 2nd harmonic, resid={resid2:.2e}"
    assert relAB < 2e-2, f"AGE salient reluctance torque off brute by {relAB:.2e}"
    # IPM: PM switches the fundamental on, reluctance 2nd-harmonic unchanged (independent superposition)
    assert g1 > 20*f1 and g1 > 1.0, f"IPM PM must switch the fundamental torque on, got {g1:.3f}"
    assert abs(g2-f2)/f2 < 5e-2, f"PM should not change the reluctance 2nd harmonic ({f2:.4f}->{g2:.4f})"
    assert relABp < 2e-2, f"AGE IPM torque off brute by {relABp:.2e}"


def main():
    validate_salient_ipm_synrm_age_matches_brute()
    print("[OK] AGE salient rotor: SynRM pure-2nd-harmonic reluctance torque + IPM fundamental(PM), "
          "Ld!=Lq saliency, ONE factorization (rotor frame) reproducing the brute meshed-gap FE.")


if __name__ == "__main__":
    main()
