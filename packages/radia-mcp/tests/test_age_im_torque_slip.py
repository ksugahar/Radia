# -*- coding: utf-8 -*-
r"""Induction-machine torque-slip via the eddy AGE, gated against the analytic Kloss equation.

Extends the AGE rotating-machine core from synchronous (test_age_pmsm_physical) to INDUCTION: a
polyphase stator MMF (forward-rotating 2-pole wave = complex spatial phasor J0*exp(i*theta)) drives a
CONDUCTING rotor (copper sleeve, sigma) across the UN-MESHED gap; the time-harmonic eddy solve at the
SLIP frequency f = s*f_sys gives the induced rotor currents and the drag torque. Sweeping slip traces
the torque-slip curve.

Physics gate (no commercial reference -- purely analytic): the single-rotor-harmonic eddy drag torque
is a single-time-constant response, so the torque-slip curve must follow the textbook KLOSS equation

    T(s) = 2 * T_max / (s/s_breakdown + s_breakdown/s)

(linear rise at low slip, breakdown peak, 1/s tail). The eddy AGE itself is validated to machine
precision against a fully-meshed complex reference in test_airgap_eddy_machine; here we lock that it
reproduces the IM torque-slip physics. Relative-reluctivity nu~=1/mur (gap nu~=1); eddy coefficient
jw*mu0*sigma.
"""
import math
import os
import sys

import numpy as np
from ngsolve import H1, BilinearForm, LinearForm, grad, dx, x, y, atan2, cos, sin, Mesh, TaskManager
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import (
    airgap_coupling, airgap_factorize, airgap_solve, airgap_torque)

mm = 1e-3
R0, Rr, Rm, Rs, Rw, REXT = 3*mm, 8*mm, 11*mm, 11.5*mm, 12.5*mm, 22*mm
MUR_FE = 1000.0
MU0 = 4e-7*math.pi
SIGMA_CU = 5.8e7
LSTK = 0.05
J0 = 5.0e6
F_SYS = 50.0
HARM = [1]


def _geo():
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=2*mm)
    g.AddCircle((0, 0), Rw, leftdomain=5, rightdomain=1, bc="wind_outer", maxh=0.4*mm)
    g.AddCircle((0, 0), Rs, leftdomain=0, rightdomain=5, bc="stator_ring", maxh=0.12*mm)
    g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=0, bc="rotor_ring", maxh=0.12*mm)
    g.AddCircle((0, 0), Rr, leftdomain=2, rightdomain=3, bc="sleeve_inner", maxh=0.3*mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=2, bc="rotor_inner", maxh=0.6*mm)
    g.SetMaterial(1, "stator_iron"); g.SetMaterial(2, "rotor_iron")
    g.SetMaterial(3, "sleeve"); g.SetMaterial(5, "winding")
    return Mesh(g.GenerateMesh(maxh=1*mm))


def test_im_torque_slip_follows_kloss():
    mesh = _geo()
    fes = H1(mesh, order=3, dirichlet="outer|rotor_inner", complex=True)
    u, v = fes.TnT()
    nu = mesh.MaterialCF({"stator_iron": 1/MUR_FE, "rotor_iron": 1/MUR_FE,
                          "sleeve": 1., "winding": 1.}, default=1.)
    sigma = mesh.MaterialCF({"sleeve": SIGMA_CU}, default=0.)
    th = atan2(y, x)
    chiw = mesh.MaterialCF({"winding": 1.}, default=0.)
    Jz = J0*(cos(th) + 1j*sin(th))                  # forward-rotating 2-pole stator MMF
    Lsrc = LinearForm(fes); Lsrc += chiw*MU0*Jz*v*dx; Lsrc.Assemble()
    coup = airgap_coupling(fes, Rm, Rs, "rotor_ring", "stator_ring", HARM)

    slips = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0])
    T = []
    with TaskManager():
        for s in slips:
            w = 2*math.pi*s*F_SYS
            a = BilinearForm(fes)
            a += nu*grad(u)*grad(v)*dx + 1j*w*MU0*sigma*u*v*dx
            a.Assemble()
            fac = airgap_factorize(a.mat, coup, fes.FreeDofs())
            gfu = airgap_solve(fac, fes, dirichlet_cf=None, source_lf=Lsrc)
            T.append(abs(airgap_torque(coup, gfu, axial_length=LSTK).real))
    T = np.array(T)

    # Kloss fit: 1/T = ca*s + cb/s  ->  s_bk = sqrt(cb/ca), Tmax = 1/(2 sqrt(ca cb))
    A = np.vstack([slips, 1.0/slips]).T
    (ca, cb), *_ = np.linalg.lstsq(A, 1.0/T, rcond=None)
    s_bk = math.sqrt(cb/ca); Tmax = 1.0/(2*math.sqrt(ca*cb))
    Tk = 2*Tmax/(slips/s_bk + s_bk/slips)
    resid = np.sqrt(np.mean((T - Tk)**2))/np.max(T)
    ipk = int(np.argmax(T))
    print(f"IM torque-slip: Tmax={Tmax*1e3:.3f} mN.m @ s_bk={s_bk:.3f}; Kloss resid={resid:.2e}; "
          f"T(low s)={T[0]*1e3:.3f}, peak@s={slips[ipk]:.2f}, T(s=1)={T[-1]*1e3:.3f} mN.m")

    assert T[0] < T[ipk], "torque must rise from low slip toward breakdown"
    assert 0 < ipk < len(slips)-1 or s_bk < slips[-1], "breakdown peak must be interior / in range"
    assert 0.03 < s_bk < 0.95, f"breakdown slip {s_bk:.3f} out of physical range"
    assert resid < 5e-2, f"torque-slip off the Kloss curve by {resid:.2e}"


def main():
    test_im_torque_slip_follows_kloss()
    print("[OK] AGE induction machine: forward-rotating stator MMF + conducting rotor across the "
          "un-meshed gap, slip-frequency eddy solve -> torque-slip curve follows the analytic Kloss law.")


if __name__ == "__main__":
    main()
