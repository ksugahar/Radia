# -*- coding: utf-8 -*-
r"""MULTI-POLE machine + DISTRIBUTED-WINDING factor for the AGE.

(A) MULTI-POLE: a 4-pole (p=2) SPM.  The PM magnetisation alternates as sgn(cos(p*theta)) and the gap
    field lives on harmonics n = p*(2k+1) = 2, 6, 10, ..., so the AGE coupling just uses those
    harmonics.  Validated: open-circuit |B_r,p| and the loaded torque (a clean T=Tmax*sin(delta) law)
    reproduce the brute meshed-gap FE -- the AGE is pole-count agnostic (the rest of the suite is p=1).

(B) DISTRIBUTED-WINDING factor: a real 3-phase single-layer full-pitch winding (q slots/pole/phase)
    has phase-a slots in the electrical 60-deg belts [0,60) (+) and [180,240) (-).  The winding factor
    of MMF/EMF harmonic n is the normalised slot-phasor sum  kw_n = |sum_k s_k exp(-i n p theta_k)|/N_a,
    which must equal the textbook distribution factor  kd_n = |sin(n q gamma/2)/(q sin(n gamma/2))|
    (gamma = electrical slot pitch), including the harmonic suppression (kd_5, kd_7 << kd_1).
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction,
                     grad, dx, x, y, sqrt, atan2, cos, sin, IfPos, Integrate, Mesh, TaskManager)
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
P = 2
HARM = [2, 6, 10, 14]


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


def _src(fes, mesh, j0, delta):
    w = fes.TestFunction(); th = atan2(y, x)
    chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
    sgn = IfPos(cos(P*th), 1.0, -1.0)
    chiw = mesh.MaterialCF({"winding": 1.0}, default=0.0)
    phi = delta - math.pi/2
    L = LinearForm(fes)
    L += chim*(1.0/MUR_PM)*(BR*sgn*cos(th)*grad(w)[1] - BR*sgn*sin(th)*grad(w)[0])*dx
    L += chiw*MU0*(j0*cos(P*th - phi))*w*dx
    L.Assemble(); return L


def _brp(B, mesh, r, p, n=480):
    vals = []
    for k in range(n):
        ph = 2*math.pi*k/n
        bx, by = B(mesh(r*math.cos(ph), r*math.sin(ph)))
        bx = bx.real if hasattr(bx, "real") else bx
        by = by.real if hasattr(by, "real") else by
        vals.append(bx*math.cos(ph) + by*math.sin(ph))
    return 2*abs(np.fft.rfft(np.array(vals))[p])/n


def _winding_factor(p, q, n):
    Q = 2*p*3*q
    gamma = 2*math.pi*p/Q
    slots = []
    for i in range(Q):
        ed = round(360*p*i/Q) % 360
        if 0 <= ed < 60:
            slots.append((2*math.pi*i/Q, +1.0))
        elif 180 <= ed < 240:
            slots.append((2*math.pi*i/Q, -1.0))
    Na = len(slots)
    kw = abs(sum(s*np.exp(-1j*n*p*th) for th, s in slots))/Na
    kd = abs(math.sin(n*q*gamma/2)/(q*math.sin(n*gamma/2)))
    return kw, kd, Na, Q


def test_multipole_age_matches_brute():
    """4-pole (p=2) SPM: AGE (harmonics n=p(2k+1)) reproduces the brute meshed-gap field & torque."""
    ma = _geo(False)
    fa = H1(ma, order=3, dirichlet="outer|rotor_inner", complex=True)
    ua, wa = fa.TnT(); aa = BilinearForm(fa); aa += _nu(ma)*grad(ua)*grad(wa)*dx
    mb = _geo(True)
    fb = H1(mb, order=3, dirichlet="outer|rotor_inner", complex=False)
    ub, wb = fb.TnT(); ab = BilinearForm(fb); ab += _nu(mb)*grad(ub)*grad(wb)*dx
    r_cf = sqrt(x*x+y*y)
    with TaskManager():
        aa.Assemble(); ab.Assemble()
        coup = airgap_coupling(fa, Rm, Rs, "rotor_ring", "stator_ring", HARM)
        fac = airgap_factorize(aa.mat, coup, fa.FreeDofs())
        Kb = ab.mat.Inverse(fb.FreeDofs(), inverse="umfpack")
        gapdx = dx(definedon=mb.Materials("gap"))
        r_eval = Rm - 0.05*mm                          # rotor surface (meshed in both)
        gb = GridFunction(fb); gb.vec.data = Kb*_src(fb, mb, 0.0, 0.0).vec
        b_surf = _brp(CoefficientFunction((grad(gb)[1], -grad(gb)[0])), mb, r_eval, P)
        ga = airgap_solve(fac, fa, source_lf=_src(fa, ma, 0.0, 0.0))
        a_surf = _brp(CoefficientFunction((grad(ga)[1], -grad(ga)[0])), ma, r_eval, P)
        rel_field = abs(a_surf-b_surf)/b_surf

        deltas = np.linspace(0, 2*math.pi, 13)
        Ta, Tb = [], []
        for d in deltas:
            g = airgap_solve(fac, fa, source_lf=_src(fa, ma, 6e6, d))
            Ta.append(airgap_torque(coup, g, axial_length=LSTK).real)
            gbb = GridFunction(fb); gbb.vec.data = Kb*_src(fb, mb, 6e6, d).vec
            B = CoefficientFunction((grad(gbb)[1], -grad(gbb)[0]))
            Br_ = (B[0]*x+B[1]*y)/r_cf; Bth = (-B[0]*y+B[1]*x)/r_cf
            Tb.append((LSTK/(MU0*(Rs-Rm)))*Integrate(r_cf*Br_*Bth*gapdx, mb))
    Ta, Tb = np.array(Ta), np.array(Tb)
    s = np.sin(deltas); Tmax = np.sum(Tb*s)/np.sum(s*s)
    resid = np.sqrt(np.mean((Tb-Tmax*s)**2))/np.max(np.abs(Tb))
    rel_T = np.max(np.abs(Ta-Tb))/np.max(np.abs(Tb))
    print(f"(A) 4-pole: |B_r,p2| AGE-vs-brute={rel_field:.2e}; loaded torque peak={np.max(np.abs(Tb))*1e3:.2f} "
          f"mN.m, sin-law resid={resid:.2e}, AGE-vs-brute={rel_T:.2e}")
    assert rel_field < 1e-3, f"4-pole open-circuit field AGE off brute by {rel_field:.2e}"
    assert resid < 5e-3, f"4-pole torque not a clean sin(delta) law, resid={resid:.2e}"
    assert rel_T < 5e-3, f"4-pole loaded torque AGE off brute by {rel_T:.2e}"


def test_distributed_winding_factor():
    """Slot-phasor winding factor kw == analytic distribution factor kd for n=1,3,5,7."""
    worst = 0.0
    for (p, q) in [(2, 2), (2, 3), (1, 4)]:
        for n in (1, 3, 5, 7):
            kw, kd, Na, Q = _winding_factor(p, q, n)
            if n == 1:
                assert Na == 2*p*q, f"phase-a slot count {Na} != 2pq={2*p*q} (p={p},q={q})"
            rel = abs(kw-kd)/max(kd, 1e-12)
            worst = max(worst, rel)
            print(f"  p={p} q={q} n={n}: kw={kw:.5f} kd={kd:.5f} (rel {rel:.1e})")
    assert worst < 1e-9, f"winding factor kw off the analytic kd by {worst:.2e}"


def main():
    test_multipole_age_matches_brute()
    test_distributed_winding_factor()
    print("[OK] AGE multi-pole (4-pole p=2, n=p(2k+1) harmonics) == brute meshed-gap; distributed "
          "winding factor kw == analytic kd (harmonic suppression) to machine precision.")


if __name__ == "__main__":
    main()
