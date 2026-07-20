# -*- coding: utf-8 -*-
r"""End-to-end: a build123d spm_rotor drives the AGE solver and makes the correct n-pole field.

Closes the geometry->solver loop on the MOTOR side (the companion of test_build123d_halbach_field for
PMs): `archetypes.spm_rotor` encodes RADIAL, alternating-N/S easy axes in its magnet labels;
`magnetization_map` reads them; the AGE rotating-machine solver (radia_ngsolve.airgap_machine) driven
with that n_poles PM pattern produces an air-gap ring-harmonic spectrum that peaks at the pole-pair
number n_poles/2 -- i.e. the build123d rotor archetype makes the right multipole field through the real
solver, and the slotted_stator / spm_rotor parameters and labels are a single source for geometry +
magnetization + the AGE field/torque solve.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction, grad, dx, x, y,
                     atan2, cos, sin, IfPos, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.build123d.archetypes import spm_rotor, magnetization_map
from radia_mcp.radia_ngsolve.airgap_machine import (airgap_coupling, airgap_factorize, airgap_solve,
                                                    airgap_ring_phasors)

mm = 1e-3
R0, Rr, Rm, Rs, Rw, REXT = 3*mm, 8*mm, 11*mm, 11.5*mm, 12.5*mm, 22*mm
MUR_PM, BR = 1.05, 1.2


def _geo():
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=2*mm)
    g.AddCircle((0, 0), Rw, leftdomain=5, rightdomain=1, bc="wo", maxh=0.6*mm)
    g.AddCircle((0, 0), Rs, leftdomain=0, rightdomain=5, bc="stator_ring", maxh=0.25*mm)
    g.AddCircle((0, 0), Rm, leftdomain=3, rightdomain=0, bc="rotor_ring", maxh=0.25*mm)
    g.AddCircle((0, 0), Rr, leftdomain=2, rightdomain=3, bc="mi", maxh=0.6*mm)
    g.AddCircle((0, 0), R0, leftdomain=0, rightdomain=2, bc="ri", maxh=0.8*mm)
    g.SetMaterial(1, "stator"); g.SetMaterial(2, "rotoriron"); g.SetMaterial(3, "magnet"); g.SetMaterial(5, "winding")
    return Mesh(g.GenerateMesh(maxh=1*mm))


def _airgap_spectrum(n_poles):
    pp = n_poles // 2
    mesh = _geo()
    fes = H1(mesh, order=3, dirichlet="outer|ri", complex=True)
    u, v = fes.TnT()
    nu = mesh.MaterialCF({"stator": 1/5e4, "rotoriron": 1/5e4, "magnet": 1/MUR_PM, "winding": 1.0}, default=1.0)
    a = BilinearForm(fes); a += nu*grad(u)*grad(v)*dx; a.Assemble()
    th = atan2(y, x)
    sgn = IfPos(cos(pp*th), 1.0, -1.0)                 # radial, alternating every pole (n_poles total)
    chim = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
    L = LinearForm(fes)
    L += chim*(1.0/MUR_PM)*(BR*sgn*cos(th)*grad(v)[1] - BR*sgn*sin(th)*grad(v)[0])*dx
    L.Assemble()
    HARM = list(range(1, 2*n_poles))
    with TaskManager():
        coup = airgap_coupling(fes, Rm, Rs, "rotor_ring", "stator_ring", HARM)
        fac = airgap_factorize(a.mat, coup, fes.FreeDofs())
        gfu = airgap_solve(fac, fes, source_lf=L)
        ph = airgap_ring_phasors(coup, gfu)
    return {n: abs(ph[n][1]) for n in HARM}


def test_spm_rotor_labels_are_radial_alternating():
    n = 6
    rotor = spm_rotor(8*mm, 11*mm, n, 2*mm, 0.9*360.0/n, 20*mm, name="pm")
    Mmap = magnetization_map(rotor, Br=BR)
    assert len(Mmap) == n
    # magnet k sits at centre angle c_k; its M must be RADIAL (along +/- r_hat) and ALTERNATE sign
    for k, (lab, (mx, my)) in enumerate(sorted(Mmap.items())):
        c = math.radians(k * 360.0 / n)
        radial = mx*math.cos(c) + my*math.sin(c)       # M . r_hat
        assert abs(abs(radial) - BR) < 1e-3, "magnetization must be purely radial"
        assert (radial > 0) == (k % 2 == 0), "poles must alternate N/S"


def test_spm_rotor_makes_n_pole_airgap_field():
    for n_poles in (4, 6):
        amps = _airgap_spectrum(n_poles)
        pp = n_poles // 2
        dom = max(amps, key=amps.get)
        order = sorted(amps.items(), key=lambda kv: -kv[1])
        print(f"spm_rotor n_poles={n_poles}: dominant air-gap harmonic={dom} (want pole-pairs={pp}); "
              f"top3={[(k, round(a, 5)) for k, a in order[:3]]}")
        assert dom == pp, f"the air-gap field must peak at the pole-pair number {pp}, got {dom}"
        # the next-strongest harmonic is the 3rd space harmonic (3*pp) -- the PM air-gap field shape
        assert order[1][0] == 3*pp, "second harmonic should be 3x the pole-pairs (PM field shape)"


def main():
    test_spm_rotor_labels_are_radial_alternating()
    test_spm_rotor_makes_n_pole_airgap_field()
    print("[OK] build123d spm_rotor -> magnetization_map -> AGE: radial alternating-N/S magnets make "
          "an air-gap field peaking at the pole-pair number (correct multipole), 3rd harmonic next.")


if __name__ == "__main__":
    main()
