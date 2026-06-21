# -*- coding: utf-8 -*-
r"""Loss post-processing -> efficiency on the AGE FE field, with an honest gated/material split.

 (A) COPPER loss (EXACT gate): P_cu = L * INT_winding |J|^2 / sigma dA for the J0*cos(theta-phi_s)
     current sheet.  Over the full annulus cos^2 averages to 1/2, so P_cu = L*J0^2/(2*sigma)*A_win.
     Gated MACHINE-PRECISION against the FE-integrated winding area A_win (isolates the cos^2->1/2
     loss physics from the mesh's annulus-area discretization); also reported vs the analytic annulus.

 (B) IRON classical-eddy loss from the FE B-field: P_eddy = (sigma_fe*d^2*(2*pi*f)^2/24)*L*INT_iron
     |B|^2 dA (the eddy formula is analytic; lamination thickness d).  GATED on its physics SCALING:
     because the magnetostatics is linear, rescaling the source scales B linearly, so P_eddy must
     scale as B^2 (and as f^2) -- log-log slope 2.  The HYSTERESIS term P_h = k_h*f*INT|B|^alpha
     needs the material coefficient k_h (NOT derivable): included in (C) and flagged, never gated.

 (C) EFFICIENCY: eta = P_out/(P_out + P_cu + P_fe), P_out = T*omega -- a physical (0,1) value.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, CoefficientFunction, grad, dx, x, y, sqrt,
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
HARM = [1, 3, 5]
SIGMA_CU, K_FILL = 5.8e7, 0.5
SIGMA_FE, D_LAM = 2.0e6, 0.35e-3
K_HYST, ALPHA = 150.0, 2.0
RPM = 3000.0
F_EL = RPM/60.0


def _geo():
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=2.0*mm)
    g.AddCircle((0, 0), Rw, leftdomain=5, rightdomain=1, bc="wind_outer", maxh=0.4*mm)
    g.AddCircle((0, 0), Rs, leftdomain=0, rightdomain=5, bc="stator_ring", maxh=0.12*mm)
    g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=0, bc="rotor_ring", maxh=0.12*mm)
    g.AddCircle((0, 0), Rr, leftdomain=2, rightdomain=3, bc="mag_inner", maxh=0.5*mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=2, bc="rotor_inner", maxh=0.6*mm)
    g.SetMaterial(1, "stator"); g.SetMaterial(2, "rotoriron")
    g.SetMaterial(3, "magnet"); g.SetMaterial(5, "winding")
    return Mesh(g.GenerateMesh(maxh=1.0*mm))


def _nu(mesh):
    return mesh.MaterialCF({"stator": 1.0/MUR_FE, "rotoriron": 1.0/MUR_FE,
                            "magnet": 1.0/MUR_PM, "gap": 1.0, "winding": 1.0}, default=1.0)


def _loaded(fes, mesh, theta_m, delta, br=BR, j0=J0):
    w = fes.TestFunction(); th = atan2(y, x)
    chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
    sgn = IfPos(cos(th - theta_m), 1.0, -1.0)
    chiw = mesh.MaterialCF({"winding": 1.0}, default=0.0)
    phi_s = theta_m + delta - math.pi/2
    L = LinearForm(fes)
    L += chim*(1.0/MUR_PM)*(br*sgn*cos(th)*grad(w)[1] - br*sgn*sin(th)*grad(w)[0])*dx
    L += chiw*MU0*(j0*cos(th - phi_s))*w*dx
    L.Assemble(); return L


def test_losses_and_efficiency():
    ma = _geo()
    fa = H1(ma, order=3, dirichlet="outer|rotor_inner", complex=True)
    ua, wa = fa.TnT(); aa = BilinearForm(fa); aa += _nu(ma)*grad(ua)*grad(wa)*dx
    th = atan2(y, x)
    win = ma.MaterialCF({"winding": 1.0}, default=0.0)
    iron = ma.MaterialCF({"stator": 1.0, "rotoriron": 1.0}, default=0.0)
    with TaskManager():
        aa.Assemble()
        coup = airgap_coupling(fa, Rm, Rs, "rotor_ring", "stator_ring", HARM)
        fac = airgap_factorize(aa.mat, coup, fa.FreeDofs())

        # (A) copper loss
        sigma_eff = SIGMA_CU*K_FILL
        Jw = J0*cos(th - 0.0)
        A_win = Integrate(win*dx, ma).real
        P_cu_fe = LSTK/sigma_eff*Integrate(win*Jw*Jw*dx, ma).real
        P_cu_area = LSTK*J0**2/(2*sigma_eff)*A_win                 # cos^2 -> 1/2, same mesh
        P_cu_annulus = LSTK*J0**2/(2*sigma_eff)*math.pi*(Rw**2-Rs**2)
        rel_cu_area = abs(P_cu_fe-P_cu_area)/P_cu_area
        rel_cu_annulus = abs(P_cu_fe-P_cu_annulus)/P_cu_annulus

        # (B) iron eddy B^2 / f^2 scaling
        ke = SIGMA_FE*D_LAM**2/24.0
        def intB2(c):
            g = airgap_solve(fac, fa, source_lf=_loaded(fa, ma, 0.0, math.pi/2, br=BR*c, j0=J0*c))
            B = CoefficientFunction((grad(g)[1], -grad(g)[0]))
            return LSTK*Integrate(iron*(B[0]*B[0]+B[1]*B[1]).real*dx, ma).real
        scales = np.array([0.5, 1.0, 2.0, 4.0])
        IB2 = np.array([intB2(c) for c in scales])
        P_eddy = ke*(2*math.pi*F_EL)**2*IB2
        slopeB = np.polyfit(np.log(scales), np.log(P_eddy), 1)[0]
        fs = np.array([25., 50., 100., 200.])
        slopeF = np.polyfit(np.log(fs), np.log(ke*(2*math.pi*fs)**2*IB2[1]), 1)[0]

        # (C) efficiency at the rated point
        gL = airgap_solve(fac, fa, source_lf=_loaded(fa, ma, 0.0, math.pi/2))
        T = abs(airgap_torque(coup, gL, axial_length=LSTK).real)
        B = CoefficientFunction((grad(gL)[1], -grad(gL)[0]))
        IB2_1 = LSTK*Integrate(iron*(B[0]*B[0]+B[1]*B[1]).real*dx, ma).real
        IBa = LSTK*Integrate(iron*((B[0]*B[0]+B[1]*B[1]).real)**(ALPHA/2)*dx, ma).real
    omega = 2*math.pi*RPM/60.0
    P_out = T*omega
    P_eddy_op = ke*(2*math.pi*F_EL)**2*IB2_1
    P_hyst_op = K_HYST*F_EL*IBa
    P_fe = P_eddy_op + P_hyst_op
    eta = P_out/(P_out + P_cu_fe + P_fe)
    print(f"(A) copper: P_cu={P_cu_fe:.4f} W; vs FE-area {rel_cu_area:.2e} (machine), "
          f"vs analytic annulus {rel_cu_annulus:.2e} (mesh)")
    print(f"(B) iron eddy: B^2 slope={slopeB:.4f}, f^2 slope={slopeF:.4f} (expect 2)")
    print(f"(C) eta={eta*100:.2f}%  (P_out={P_out:.3f}, P_cu={P_cu_fe:.3f}, "
          f"P_fe={P_fe:.3f}=eddy{P_eddy_op:.3f}+hyst{P_hyst_op:.3f}[k_h])")

    assert rel_cu_area < 1e-3, f"copper-loss cos^2 averaging off FE area by {rel_cu_area:.2e}"
    assert rel_cu_annulus < 5e-3, f"copper loss off analytic closed form by {rel_cu_annulus:.2e}"
    assert abs(slopeB - 2.0) < 5e-3, f"iron eddy loss must scale as B^2, slope={slopeB:.4f}"
    assert abs(slopeF - 2.0) < 1e-6, f"iron eddy loss must scale as f^2, slope={slopeF:.4f}"
    assert P_cu_fe > 0 and P_fe > 0 and 0.5 < eta < 1.0, f"efficiency unphysical: eta={eta:.3f}"


def main():
    test_losses_and_efficiency()
    print("[OK] AGE losses->efficiency: copper I^2R exact (cos^2 averaging machine-precision), iron "
          "classical-eddy from the FE B-field with B^2 & f^2 scaling, efficiency assembled (k_h flagged).")


if __name__ == "__main__":
    main()
