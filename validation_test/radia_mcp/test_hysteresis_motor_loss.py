# -*- coding: utf-8 -*-
r"""Motor iron loss from the REAL play-model loop area (no fitted k_h), on an AGE rotor sweep.

The per-point B(theta) from a rotor-angle sweep (here a PM-only SPM cogging-style sweep, the AGE
machinery of the committed motor tests) is fed through the B-input vector Play operator; the
hysteresis loss density at each point is the actual loop area ∮ H.dB.  Because the stator-iron B
ROTATES as the rotor turns (an elliptical B-path), this captures the rotational / minor-loop loss that
a scalar |B|pk-Steinmetz cannot -- so the play-waveform iron loss is strictly LARGER than the
scalar-alternating estimate.  Integrated over the iron, it is a physical iron-loss number from the
physical hysteresis loop, replacing the fitted k_h of coreloss.steinmetz_loss_density.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction,
                     grad, dx, x, y, sqrt, atan2, cos, sin, IfPos, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import (
    airgap_coupling, airgap_factorize, airgap_solve)
from radia_mcp.radia_ngsolve.hysteresis import steinmetz_cells

mm = 1e-3
R0, Rr, Rm, Rs, Rw, REXT = 3*mm, 8*mm, 11*mm, 11.5*mm, 12.5*mm, 22*mm
MUR_FE, MUR_PM, BR = 5.0e4, 1.05, 1.2
MU0 = 4e-7*math.pi
LSTK = 0.05
HARM = [1, 3, 5]
F_MECH = 50.0


def _geo():
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=2*mm)
    g.AddCircle((0, 0), Rw, leftdomain=5, rightdomain=1, bc="wo", maxh=0.4*mm)
    g.AddCircle((0, 0), Rs, leftdomain=0, rightdomain=5, bc="stator_ring", maxh=0.12*mm)
    g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=0, bc="rotor_ring", maxh=0.12*mm)
    g.AddCircle((0, 0), Rr, leftdomain=2, rightdomain=3, bc="mi", maxh=0.5*mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=2, bc="ri", maxh=0.6*mm)
    g.SetMaterial(1, "stator"); g.SetMaterial(2, "rotoriron"); g.SetMaterial(3, "magnet"); g.SetMaterial(5, "winding")
    return Mesh(g.GenerateMesh(maxh=1*mm))


def _pm(fes, mesh, theta_m):
    w = fes.TestFunction(); th = atan2(y, x)
    chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0); sgn = IfPos(cos(th-theta_m), 1.0, -1.0)
    L = LinearForm(fes)
    L += chim*(1.0/MUR_PM)*(BR*sgn*cos(th)*grad(w)[1] - BR*sgn*sin(th)*grad(w)[0])*dx
    L.Assemble(); return L


def test_motor_iron_loss_from_play_loop():
    ma = _geo()
    fa = H1(ma, order=3, dirichlet="outer|ri", complex=True)
    ua, wa = fa.TnT(); aa = BilinearForm(fa)
    aa += ma.MaterialCF({"stator": 1/MUR_FE, "rotoriron": 1/MUR_FE, "magnet": 1/MUR_PM, "winding": 1.}, default=1.)*grad(ua)*grad(wa)*dx
    with TaskManager():
        aa.Assemble()
        coup = airgap_coupling(fa, Rm, Rs, "rotor_ring", "stator_ring", HARM)
        fac = airgap_factorize(aa.mat, coup, fa.FreeDofs())
        radii = np.linspace(13.0*mm, 20.0*mm, 5); angs = np.linspace(0, 2*math.pi, 18, endpoint=False)
        pts = [(r*math.cos(a), r*math.sin(a)) for r in radii for a in angs]
        NTH = 24
        thetas = np.linspace(0, 2*math.pi, NTH, endpoint=False)
        Bx = np.zeros((len(pts), NTH)); By = np.zeros((len(pts), NTH))
        for j, tm in enumerate(thetas):
            ga = airgap_solve(fac, fa, source_lf=_pm(fa, ma, tm))
            B = CoefficientFunction((grad(ga)[1], -grad(ga)[0]))
            for i, (pxi, pyi) in enumerate(pts):
                bx, by = B(ma(pxi, pyi)); Bx[i, j] = bx.real; By[i, j] = by.real

    model = steinmetz_cells(eta_min=0.02, eta_max=1.5, K=40, a_each=150.0)
    play_dens = np.array([model.loss_from_waveform(np.append(Bx[i], Bx[i, 0]), np.append(By[i], By[i, 0]))
                          for i in range(len(pts))])                      # real loop area per point
    Bpk = np.array([np.max(np.sqrt(Bx[i]**2 + By[i]**2)) for i in range(len(pts))])
    scalar_dens = np.array([model.loss_per_cycle(b, n=1001) for b in Bpk])  # scalar |B|pk-alternating
    A_iron = math.pi*(REXT**2 - Rw**2)
    P_play = play_dens.mean()*A_iron*LSTK*F_MECH                          # W (mean-density estimate)
    print(f"motor iron loss: {len(pts)} pts, |B|pk {Bpk.min():.2f}..{Bpk.max():.2f} T; "
          f"play-loop P={P_play:.3f} W; play/scalar density ratio mean={np.mean(play_dens/scalar_dens):.3f}")

    assert np.all(play_dens > 0), "per-point hysteresis loss must be positive"
    assert 1e-3 < P_play < 1e3, "total iron loss should be a physical wattage"
    assert np.mean(play_dens / scalar_dens) > 1.02, ("the rotating/elliptical stator B must make the "
                                                     "play-loop loss exceed the scalar |B|pk estimate")


def main():
    test_motor_iron_loss_from_play_loop()
    print("[OK] motor iron loss from the play-model loop area: per-point B(theta) AGE sweep -> real "
          "loop-area hysteresis loss, capturing the rotational/elliptical loss a scalar Steinmetz misses.")


if __name__ == "__main__":
    main()
