# -*- coding: utf-8 -*-
r"""3-PHASE TIME-WAVEFORM torque ripple via the AGE, vs the brute meshed-gap FE.

A real 3-phase distributed winding (60-deg belts -> MMF space harmonics 5, 7, 11, 13) carries balanced
sinusoidal currents i_a,b,c(t), with the square-wave PM rotor phase-locked at the max-torque operating
point (theta_m = omega_e t, p=1).  The belt MMF harmonics beat with the PM field to give the
characteristic 6k ELECTRICAL torque ripple on the DC operating torque -- the standard PMSM ripple that
a single-harmonic SINUSOIDAL winding does NOT show (test_age_synchronous_torque was ripple-free).
torque(t) over one electrical period -> FFT.

ONE factorization: the belts are fixed current SOURCES and the stator iron is smooth, so K is fixed;
only the PM source and the phase currents vary in time.  Validated: torque(t) AGE (mesh-free gap) vs
brute (meshed gap) per step; the ripple spectrum is the 6k harmonics (6, 12), with the non-6k orders
(1, 2, 3, ...) essentially zero.
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
HARM = list(range(1, 16))
C30 = math.cos(math.pi/6)
NT = 24


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


def _belt(th, shift):
    return IfPos(cos(th - shift) - C30, 1.0, 0.0) - IfPos(-cos(th - shift) - C30, 1.0, 0.0)


def _src(fes, mesh, wt):
    w = fes.TestFunction(); th = atan2(y, x)
    ia, ib, ic = math.cos(wt), math.cos(wt - 2*math.pi/3), math.cos(wt + 2*math.pi/3)
    Jt = _belt(th, 0.0)*ia + _belt(th, 2*math.pi/3)*ib + _belt(th, 4*math.pi/3)*ic
    chiw = mesh.MaterialCF({"winding": 1.0}, default=0.0)
    chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
    sgn = IfPos(cos(th - wt), 1.0, -1.0)               # PM phase-locked at theta_m = wt (max torque)
    L = LinearForm(fes)
    L += chiw*MU0*(J0*Jt)*w*dx
    L += chim*(1.0/MUR_PM)*(BR*sgn*cos(th)*grad(w)[1] - BR*sgn*sin(th)*grad(w)[0])*dx
    L.Assemble(); return L


def test_threephase_torque_ripple_age_matches_brute():
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
        fac = airgap_factorize(aa.mat, coup, fa.FreeDofs())          # ONE factorization
        Kb = ab.mat.Inverse(fb.FreeDofs(), inverse="umfpack")
        gapdx = dx(definedon=mb.Materials("gap"))
        wts = np.linspace(0, 2*math.pi, NT, endpoint=False)
        Ta, Tb = [], []
        for wt in wts:
            ga = airgap_solve(fac, fa, source_lf=_src(fa, ma, wt))
            Ta.append(airgap_torque(coup, ga, axial_length=LSTK).real)
            gb = GridFunction(fb); gb.vec.data = Kb*_src(fb, mb, wt).vec
            B = CoefficientFunction((grad(gb)[1], -grad(gb)[0]))
            Br_ = (B[0]*x+B[1]*y)/r_cf; Bth = (-B[0]*y+B[1]*x)/r_cf
            Tb.append((LSTK/(MU0*(Rs-Rm)))*Integrate(r_cf*Br_*Bth*gapdx, mb))
    Ta, Tb = np.array(Ta)*1e3, np.array(Tb)*1e3
    dc = np.mean(Tb)
    F = np.fft.rfft(Tb)/NT
    amp = 2*np.abs(F)
    ripple_pp = (Tb.max()-Tb.min())/abs(dc)
    relAB = np.max(np.abs(Ta-Tb))/np.max(np.abs(Tb))
    ac = amp[1:]; kdom = int(np.argmax(ac))+1
    non6k = max(amp[n] for n in range(1, len(amp)) if n % 6 != 0)
    print(f"3-phase ripple: DC={dc:.3f} mN.m, ripple={ripple_pp*100:.2f}% pk-pk, AGE-vs-brute={relAB:.2e}; "
          f"dominant ripple order={kdom} ({amp[6]/abs(dc)*100:.2f}% n=6, {amp[12]/abs(dc)*100:.2f}% n=12); "
          f"non-6k max={non6k/amp[6]:.2e} of n=6")

    assert abs(dc) > 1e-3, "operating-point mean torque should be physical"
    assert ripple_pp > 0.01, f"3-phase belt winding must show a real torque ripple, got {ripple_pp*100:.2f}%"
    assert kdom == 6, f"dominant torque ripple must be the 6th electrical harmonic, got {kdom}"
    assert non6k < 0.1*amp[6], f"ripple must be 6k-only (3-phase), non-6k={non6k:.2e} vs n=6={amp[6]:.2e}"
    assert relAB < 1e-2, f"AGE 3-phase torque waveform off brute by {relAB:.2e}"


def main():
    test_threephase_torque_ripple_age_matches_brute()
    print("[OK] AGE 3-phase time-waveform: real belt winding + phase-locked PM, one factorization, "
          "reproduces the brute meshed-gap torque(t) and its characteristic 6k-electrical ripple.")


if __name__ == "__main__":
    main()
