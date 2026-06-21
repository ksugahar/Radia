# -*- coding: utf-8 -*-
r"""NONLINEAR x SLIDING-BAND: salient rotor + slotted stator + SATURATING iron, per-angle x Picard.

The hardest regime in the suite: BOTH angle-structured (salient rotor nu~(theta-theta_m) rotating +
slotted stator nu~(theta) fixed -> per-angle re-assembly) AND nonlinear (saturating iron nu~(|B|) ->
Picard at each angle).  So at every rotor angle a Picard loop runs, each iteration re-assembling and
re-factorizing K with the combined structural-and-saturating reluctivity.  The AGE gap coupling is
built ONCE (rings stay circular).  Combined model: iron fraction t(theta;theta_m) interpolates
iron<->air (slot/barrier), and the iron part saturates

    nu~ = t * [1/mur + (1-1/mur)*s(|B|)] + (1-t)*1,   s(|B|) = |B|^M/(|B|^M + B_knee^M).

Validated: AGE (analytic gap, per-angle Picard) == brute (meshed gap, per-angle Picard) at each rotor
angle, both Picard solves converge, and the saturation visibly reduces the torque vs the linear
sliding-band machine.  Harmonic-hungry (slot x saliency) so several gap harmonics are used.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction, Norm,
                     grad, dx, x, y, sqrt, atan2, cos, sin, exp, Integrate, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import (
    airgap_coupling, airgap_factorize, airgap_solve, airgap_torque)
from radia_mcp.radia_ngsolve.airgap_motor_workflow import age_motor_nonlinear_solve

mm = 1e-3
R0, Rm, Rs, Rslot, REXT = 3*mm, 11*mm, 11.5*mm, 13.5*mm, 20.0*mm
MUR_FE = 1000.0
MU0 = 4e-7*math.pi
LSTK = 0.05
J0 = 2.0e7
NSLOTS = 6
HARM = list(range(1, 11))
DELTA = math.pi/4
B_KNEE, MSAT = 1.45, 3
RELAX, NITER, TOL = 0.35, 90, 1e-6


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


def _iron_frac(mesh, theta_m):
    th = atan2(y, x)
    t_slot = 1.0/(1.0 + exp(-3.0*cos(NSLOTS*th)))
    t_rot = cos(th - theta_m)**2
    chi_slot = mesh.MaterialCF({"slotzone": 1.0}, default=0.0)
    chi_rot = mesh.MaterialCF({"rotor": 1.0}, default=0.0)
    chi_stat = mesh.MaterialCF({"stator": 1.0}, default=0.0)
    return chi_slot*t_slot + chi_rot*t_rot + chi_stat*1.0


def _nu_lin(mesh, theta_m):
    t = _iron_frac(mesh, theta_m)
    return t*(1.0/MUR_FE) + (1.0 - t)*1.0


def _nu_nl(mesh, theta_m, gfu):
    t = _iron_frac(mesh, theta_m)
    s = Norm(CoefficientFunction((grad(gfu)[1], -grad(gfu)[0])))**MSAT
    s = s/(s + B_KNEE**MSAT)
    return t*((1.0/MUR_FE) + (1.0 - 1.0/MUR_FE)*s) + (1.0 - t)*1.0


def _cur(fes, mesh, theta_m):
    w = fes.TestFunction(); th = atan2(y, x)
    chiw = mesh.MaterialCF({"slotzone": 1.0}, default=0.0)
    L = LinearForm(fes)
    L += chiw*MU0*(J0*cos(th - (theta_m + DELTA)))*w*dx
    L.Assemble(); return L


def test_nonlinear_slotted_salient_age_matches_brute():
    ma = _geo(False)
    fa = H1(ma, order=3, dirichlet="outer|rotor_inner", complex=True); ua, wa = fa.TnT()
    mb = _geo(True)
    fb = H1(mb, order=3, dirichlet="outer|rotor_inner", complex=False); ub, wb = fb.TnT()
    r_cf = sqrt(x*x+y*y)
    with TaskManager():
        coup = airgap_coupling(fa, Rm, Rs, "rotor_ring", "stator_ring", HARM)   # ONCE
        gapdx = dx(definedon=mb.Materials("gap"))

        def bt(g):
            B = CoefficientFunction((grad(g)[1], -grad(g)[0]))
            Br_ = (B[0]*x+B[1]*y)/r_cf; Bth = (-B[0]*y+B[1]*x)/r_cf
            return (LSTK/(MU0*(Rs-Rm)))*Integrate(r_cf*Br_*Bth*gapdx, mb)*1e3

        Tl, Ta, Tb, conv = [], [], [], []
        for tm in np.linspace(0, 2*math.pi/NSLOTS, 2, endpoint=False):
            # AGE linear (reference)
            al = BilinearForm(fa); al += _nu_lin(ma, tm)*grad(ua)*grad(wa)*dx; al.Assemble()
            facl = airgap_factorize(al.mat, coup, fa.FreeDofs())
            gl = airgap_solve(facl, fa, source_lf=_cur(fa, ma, tm))
            Tl.append(airgap_torque(coup, gl, axial_length=LSTK).real*1e3)
            # AGE nonlinear (per-angle Picard)
            def bil(nu): a = BilinearForm(fa); a += nu*grad(ua)*grad(wa)*dx; a.Assemble(); return a
            gnl, info = age_motor_nonlinear_solve(fa, coup, bil, nu_cf_fn=lambda gf, _tm=tm: _nu_nl(ma, _tm, gf),
                                                  source_lf=_cur(fa, ma, tm), niter=NITER, tol=TOL, relax=RELAX)
            Ta.append(airgap_torque(coup, gnl, axial_length=LSTK).real*1e3)
            # brute nonlinear (per-angle Picard, meshed gap)
            Lb = _cur(fb, mb, tm); nu = _nu_lin(mb, tm); gp = GridFunction(fb); gc = GridFunction(fb)
            cb = False
            for k in range(NITER):
                a = BilinearForm(fb); a += nu*grad(ub)*grad(wb)*dx; a.Assemble()
                sol = a.mat.Inverse(fb.FreeDofs(), inverse="umfpack")*Lb.vec
                gc.vec.data = RELAX*sol + (1.0-RELAX)*gp.vec
                d = gc.vec.CreateVector(); d.data = gc.vec - gp.vec
                nd = np.linalg.norm(d.FV().NumPy()); nc = np.linalg.norm(gc.vec.FV().NumPy())
                gp.vec.data = gc.vec; nu = _nu_nl(mb, tm, gc)
                if nc > 0 and nd/nc < TOL: cb = True; break
            Tb.append(bt(gc)); conv.append(info["converged"] and cb)
    Tl, Ta, Tb = np.array(Tl), np.array(Ta), np.array(Tb)
    relAB = np.max(np.abs(Ta-Tb))/np.max(np.abs(Tb))
    ipk = int(np.argmax(np.abs(Tl)))                      # saturation matters at the loaded (peak) angle
    sat = Ta[ipk]/Tl[ipk]
    print(f"nonlinear x sliding-band: AGE={np.array2string(Ta, precision=3)} vs "
          f"brute={np.array2string(Tb, precision=3)} mN.m; AGE-vs-brute={relAB:.2e}; "
          f"saturation T_nl/T_lin={sat:.3f} @ peak angle; converged={all(conv)}")

    assert all(conv), "every per-angle Picard solve (AGE and brute) must converge"
    assert relAB < 4e-2, f"per-angle nonlinear AGE off brute by {relAB:.2e}"
    assert sat < 0.95, f"saturation must visibly reduce the loaded torque, T_nl/T_lin={sat:.3f}"
    assert np.max(np.abs(Tb)) > 1e-3, "torque should be physical"


def main():
    test_nonlinear_slotted_salient_age_matches_brute()
    print("[OK] AGE nonlinear x sliding-band: salient rotor + slotted stator + saturating iron, "
          "per-angle Picard with analytic gap, matches the brute meshed-gap nonlinear solve.")


if __name__ == "__main__":
    main()
