# -*- coding: utf-8 -*-
r"""HYSTERESIS-COUPLED AGE rotating-machine solve: the Play constitutive INSIDE the AGE rotor sweep.

The committed test_hysteresis_motor_loss.py post-processes a LINEAR AGE field through the play operator
(decoupled).  Here the iron's B-H response IS the B-input vector Play model DURING the solve, with the
per-element play states carried FORWARD across rotor angles (magnetic history) -- so the field, torque,
and iron loss are all hysteretic.  Driver: radia_ngsolve.airgap_motor_workflow.age_motor_hysteresis_sweep,
which couples the play constitutive into the analytic air-gap element via the fixed-point M-source
(excess-field) method, reusing ONE AGE factorization for the whole rotation sweep (the play co-energy is
convex, so with the optimal reference nu0=(nu_max+nu_min)/2 the M-source contracts in a handful of iters).

Verified: (1) every rotor angle converges in a handful of M-source iterations; (2) a physical,
finite torque; (3) a positive (dissipative) iron loss from the REAL per-element loop area; (4) that
loop loss EXCEEDS the scalar |B|pk-alternating estimate -- the rotating/elliptical stator B (which a
scalar Steinmetz cannot see) makes the coupled hysteretic loss strictly larger.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, LinearForm, GridFunction, CoefficientFunction, grad, dx, x, y,
                     atan2, cos, sin, IfPos, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import airgap_coupling
from radia_mcp.radia_ngsolve.airgap_motor_workflow import age_motor_hysteresis_sweep
from radia_mcp.radia_ngsolve.hysteresis import PlayHysteresis

mm = 1e-3
MU0 = 4e-7 * math.pi
R0, Rr, Rm, Rs, Rw, REXT = 3*mm, 8*mm, 11*mm, 11.5*mm, 12.5*mm, 22*mm
MUR_PM, BR = 1.05, 1.2
LSTK, F_MECH, HARM = 0.05, 50.0, [1, 2, 3, 4, 5]


def _geo():
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=3*mm)
    g.AddCircle((0, 0), Rw, leftdomain=5, rightdomain=1, bc="wo", maxh=0.8*mm)
    g.AddCircle((0, 0), Rs, leftdomain=0, rightdomain=5, bc="stator_ring", maxh=0.35*mm)
    g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=0, bc="rotor_ring", maxh=0.35*mm)
    g.AddCircle((0, 0), Rr, leftdomain=2, rightdomain=3, bc="mi", maxh=0.8*mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=2, bc="ri", maxh=1.0*mm)
    g.SetMaterial(1, "stator"); g.SetMaterial(2, "rotoriron"); g.SetMaterial(3, "magnet"); g.SetMaterial(5, "winding")
    return Mesh(g.GenerateMesh(maxh=1.5*mm))


def test_age_hysteresis_coupled_sweep():
    # normalized play iron: reversible mur ~ 800 (a0 = 1/800) + modest dissipative cells (a~ = mu0*nu)
    eta = np.array([0.0, 0.4, 0.8, 1.2])
    a = np.array([1.0/800, 0.3/800, 0.3/800, 0.3/800])
    model = PlayHysteresis(eta=eta, a=a)

    mesh = _geo()
    fes = H1(mesh, order=1, dirichlet="outer|ri", complex=True)

    def source_fn(theta_m):
        w = fes.TestFunction(); th = atan2(y, x)
        chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0); sgn = IfPos(cos(th-theta_m), 1.0, -1.0)
        L = LinearForm(fes)
        L += chim*(1.0/MUR_PM)*(BR*sgn*cos(th)*grad(w)[1] - BR*sgn*sin(th)*grad(w)[0])*dx
        L.Assemble(); return L

    NTH = 12
    thetas = np.linspace(0, 2*math.pi, NTH, endpoint=False)
    with TaskManager():
        coupling = airgap_coupling(fes, Rm, Rs, "rotor_ring", "stator_ring", HARM)
        res = age_motor_hysteresis_sweep(
            fes, coupling, model, source_fn, thetas,
            iron_materials="stator|rotoriron", other_nu={"magnet": 1.0/MUR_PM}, default_nu=1.0,
            axial_length=LSTK, freq=F_MECH, cycles=2, msource_iters=60, tol=1e-6,
        )

    T = res["torque"]; ms = res["ms_iters"]; Bpk = res["Bpk_iron"]; ld = res["loss_density"]
    # scalar |B|pk-alternating analytic loop area per element (normalized 4 a eta (Bpk-eta) -> /mu0)
    scalar = np.zeros(len(Bpk))
    for k in range(len(eta)):
        if eta[k] > 0:
            scalar += 4.0*a[k]*eta[k]*np.maximum(Bpk-eta[k], 0.0)
    scalar /= MU0
    P_play = res["iron_loss_W"]
    # robust rotational-capture check: among loss-bearing elements, the play loop exceeds the scalar
    sig = ld > 0.2*ld.max()
    ratio_play_scalar = ld[sig] / np.maximum(scalar[sig], 1e-30)
    print(f"AGE+hysteresis sweep: {res['n_iron_elements']} iron elems, M-source iters "
          f"{ms.min()}..{ms.max()} (nu0={res['nu0']:.4f}); torque mean={T.mean()*1e3:.3f} mN.m, "
          f"pk-pk={(T.max()-T.min())*1e3:.3f}; |B|iron pk {Bpk.min():.2f}..{Bpk.max():.2f} T; "
          f"iron loss P={P_play*1e3:.2f} mW; play/scalar ratio mean(loss-bearing)={ratio_play_scalar.mean():.3f}")

    assert res["n_iron_elements"] > 100, "iron must be meshed"
    assert ms.max() <= 25, ("the M-source (with the optimal convex reference nu0) must converge in a "
                            "handful of iterations per rotor angle -- the AGE one-factorization advantage")
    assert np.all(np.isfinite(T)) and abs(T.mean()) > 1e-4, "torque must be physical (finite, nonzero)"
    assert 0.5 < Bpk.max() < 5.0, "iron must magnetise to a physical flux density"
    assert P_play > 0.0, "the coupled hysteretic solve must dissipate a positive iron loss"
    assert ratio_play_scalar.mean() > 1.02, ("the rotating/elliptical iron B must make the coupled "
                                             "play-loop loss exceed the scalar |B|pk-alternating estimate")


def main():
    test_age_hysteresis_coupled_sweep()
    print("[OK] hysteresis-coupled AGE rotor sweep: B-input Play constitutive in the iron with per-element "
          "history carried across rotor angles, M-source reusing one AGE factorization -> hysteretic "
          "torque + real-loop iron loss (rotational/elliptical capture).")


if __name__ == "__main__":
    main()
