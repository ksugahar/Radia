# -*- coding: utf-8 -*-
r"""NONLINEAR (saturating) iron x the AGE, gated against a nonlinear brute meshed-gap Picard solve.

Adds a field-dependent iron reluctivity nu~(|B|) to the SPM and solves the machine nonlinearly with
the AGE Picard wrapper (age_motor_nonlinear_solve): the gap stays analytic (air, nu~=1, linear), only
the iron FE reluctivity updates each Picard sweep.  Saturating model (relative reluctivity, Hill knee)

    nu~_iron(|B|) = 1/mur_unsat + (1 - 1/mur_unsat) * |B|^M / (|B|^M + B_knee^M)

(unsaturated 1/mur at low |B|, -> 1 = air-like above the knee).  Validated: the AGE nonlinear torque
== an independent nonlinear brute Picard (meshed gap, same nu~(|B|)); both converge; and saturation
is REAL -- the nonlinear torque is well below the linear torque (the iron clamps the flux).
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction, Norm,
                     grad, dx, x, y, sqrt, atan2, cos, sin, IfPos, Integrate, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import (
    airgap_coupling, airgap_factorize, airgap_solve, airgap_torque)
from radia_mcp.radia_ngsolve.airgap_motor_workflow import age_motor_nonlinear_solve

mm = 1e-3
R0, Rr, Rm, Rs, Rw, REXT = 3*mm, 8*mm, 11*mm, 11.5*mm, 12.5*mm, 22*mm
MUR_FE, MUR_PM, BR = 2.0e3, 1.05, 0.9
MU0 = 4e-7*math.pi
LSTK = 0.05
HARM = [1, 3, 5]
B_KNEE, MSAT = 1.9, 4
RELAX, NITER, TOL = 0.3, 120, 1e-6
THETA, J0 = 0.6, 8e6


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


def _nu_lin(mesh):
    return mesh.MaterialCF({"stator": 1.0/MUR_FE, "rotoriron": 1.0/MUR_FE,
                            "magnet": 1.0/MUR_PM, "gap": 1.0, "winding": 1.0}, default=1.0)


def _nu_nl(mesh, gfu):
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    s = Norm(B)**MSAT / (Norm(B)**MSAT + B_KNEE**MSAT)
    nu_iron = (1.0/MUR_FE) + (1.0 - 1.0/MUR_FE)*s
    chi = mesh.MaterialCF({"stator": 1.0, "rotoriron": 1.0}, default=0.0)
    base = mesh.MaterialCF({"magnet": 1.0/MUR_PM, "gap": 1.0, "winding": 1.0}, default=1.0)
    return chi*nu_iron + (1.0 - chi)*base


def _src(fes, mesh, j0):
    w = fes.TestFunction(); th = atan2(y, x)
    chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
    sgn = IfPos(cos(th - THETA), 1.0, -1.0)
    chiw = mesh.MaterialCF({"winding": 1.0}, default=0.0)
    L = LinearForm(fes)
    L += chim*(1.0/MUR_PM)*(BR*sgn*cos(th)*grad(w)[1] - BR*sgn*sin(th)*grad(w)[0])*dx
    L += chiw*MU0*(j0*cos(th))*w*dx
    L.Assemble(); return L


def test_nonlinear_saturation_age_matches_brute():
    ma = _geo(False)
    fa = H1(ma, order=3, dirichlet="outer|rotor_inner", complex=True)
    ua, wa = fa.TnT()
    mb = _geo(True)
    fb = H1(mb, order=3, dirichlet="outer|rotor_inner", complex=False)
    ub, wb = fb.TnT()
    with TaskManager():
        coup = airgap_coupling(fa, Rm, Rs, "rotor_ring", "stator_ring", HARM)

        # AGE linear (reference) torque
        al = BilinearForm(fa); al += _nu_lin(ma)*grad(ua)*grad(wa)*dx; al.Assemble()
        facl = airgap_factorize(al.mat, coup, fa.FreeDofs())
        gl = airgap_solve(facl, fa, source_lf=_src(fa, ma, J0))
        T_lin = abs(airgap_torque(coup, gl, axial_length=LSTK).real)

        # AGE nonlinear (Picard)
        def bilinear_fn(nu_cf):
            a = BilinearForm(fa); a += nu_cf*grad(ua)*grad(wa)*dx; a.Assemble(); return a
        gnl, info = age_motor_nonlinear_solve(
            fa, coup, bilinear_fn, nu_cf_fn=lambda g: _nu_nl(ma, g),
            source_lf=_src(fa, ma, J0), niter=NITER, tol=TOL, relax=RELAX)
        T_age = abs(airgap_torque(coup, gnl, axial_length=LSTK).real)

        # brute nonlinear (meshed gap, own Picard)
        Lb = _src(fb, mb, J0); nu = _nu_lin(mb)
        gp = GridFunction(fb); gc = GridFunction(fb); conv_b = False; nit_b = NITER
        for k in range(NITER):
            a = BilinearForm(fb); a += nu*grad(ub)*grad(wb)*dx; a.Assemble()
            sol = a.mat.Inverse(fb.FreeDofs(), inverse="umfpack")*Lb.vec
            gc.vec.data = RELAX*sol + (1.0-RELAX)*gp.vec
            d = gc.vec.CreateVector(); d.data = gc.vec - gp.vec
            nd = np.linalg.norm(d.FV().NumPy()); nc = np.linalg.norm(gc.vec.FV().NumPy())
            gp.vec.data = gc.vec; nu = _nu_nl(mb, gc)
            if nc > 0 and nd/nc < TOL:
                conv_b = True; nit_b = k+1; break
        r_cf = sqrt(x*x+y*y); gap = dx(definedon=mb.Materials("gap"))
        B = CoefficientFunction((grad(gc)[1], -grad(gc)[0]))
        Br_ = (B[0]*x+B[1]*y)/r_cf; Bth = (-B[0]*y+B[1]*x)/r_cf
        T_brute = abs((LSTK/(MU0*(Rs-Rm)))*Integrate(r_cf*Br_*Bth*gap, mb))

    rel = abs(T_age-T_brute)/T_brute
    sat = T_age/T_lin
    print(f"nonlinear iron: AGE={T_age*1e3:.3f} (conv={info['converged']}, {info['niter_used']} it) "
          f"vs brute={T_brute*1e3:.3f} (conv={conv_b}, {nit_b} it) mN.m; rel={rel:.2e}; "
          f"T_nl/T_lin={sat:.3f} (linear {T_lin*1e3:.3f})")

    assert info["converged"] and conv_b, "both nonlinear Picard solves must converge"
    assert rel < 5e-3, f"AGE nonlinear torque off brute nonlinear by {rel:.2e}"
    assert sat < 0.95, f"iron saturation must reduce torque below linear, T_nl/T_lin={sat:.3f}"
    assert T_age > 1e-3, "nonlinear torque should be physical"


def main():
    test_nonlinear_saturation_age_matches_brute()
    print("[OK] AGE nonlinear iron: saturating nu~(|B|) Picard with analytic gap, reproduces the "
          "nonlinear brute meshed-gap torque; saturation cuts the torque below the linear machine.")


if __name__ == "__main__":
    main()
