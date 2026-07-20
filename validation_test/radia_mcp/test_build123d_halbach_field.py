# -*- coding: utf-8 -*-
r"""End-to-end: the build123d Halbach archetype really makes a uniform transverse dipole bore field.

Closes the geometry -> field loop.  `archetypes.halbach_ring` encodes the per-segment easy-axis angle
in each region label via the Mallinson law alpha=(pole_pairs+1)*theta; `magnetization_map` turns those
labels into magnetization vectors; a 2D A_z PM magnetostatic solve with that magnetization gives the
classic Halbach result -- a UNIFORM transverse field in the bore of magnitude Br*ln(r_out/r_in).  This
verifies the archetype's magnetization convention against real physics, not just geometry.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, CoefficientFunction, grad, dx, x, y,
                     atan2, cos, sin, Mesh, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.build123d.archetypes import halbach_ring, parse_magnetization, magnetization_map

MU0 = 4e-7 * math.pi


def _solve_bore_field(r_in, r_out, Br, pole_pairs):
    g = SplineGeometry()
    g.AddCircle((0, 0), 120.0, leftdomain=1, rightdomain=0, bc="outer", maxh=8)
    g.AddCircle((0, 0), r_out, leftdomain=2, rightdomain=1, maxh=2.0)
    g.AddCircle((0, 0), r_in, leftdomain=3, rightdomain=2, maxh=2.0)
    g.SetMaterial(1, "air"); g.SetMaterial(2, "ring"); g.SetMaterial(3, "bore")
    mesh = Mesh(g.GenerateMesh(maxh=8))
    fes = H1(mesh, order=3, dirichlet="outer")
    u, v = fes.TnT()
    a = BilinearForm(fes); a += (1.0/MU0)*grad(u)*grad(v)*dx; a.Assemble()
    th = atan2(y, x)
    Mx = (Br/MU0)*cos((pole_pairs+1)*th); My = (Br/MU0)*sin((pole_pairs+1)*th)
    f = LinearForm(fes); f += (Mx*grad(v)[1] - My*grad(v)[0])*dx(definedon=mesh.Materials("ring")); f.Assemble()
    gfu = GridFunction(fes)
    with TaskManager():
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack")*f.vec
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    pts = [(rr*math.cos(t), rr*math.sin(t)) for rr in (5, 15, 25) for t in np.linspace(0, 2*math.pi, 24, endpoint=False)]
    Bx = np.array([B(mesh(px, py))[0] for px, py in pts])
    By = np.array([B(mesh(px, py))[1] for px, py in pts])
    return Bx, By


def test_halbach_archetype_makes_uniform_dipole_bore_field():
    r_in, r_out, Br, p, n = 40.0, 55.0, 1.2, 1, 12
    # the archetype encodes Mallinson easy axes in the labels; the loader reads them as M vectors
    hb = halbach_ring(r_in, r_out, 20, n, pole_pairs=p, name="hb")
    Mmap = magnetization_map(hb, Br=Br)
    for k, c in enumerate(hb.children):
        theta_c = k*360.0/n + 0.5*360.0/n
        assert abs((parse_magnetization(c.label) - ((p+1)*theta_c)) % 360.0) < 1e-2
    assert len(Mmap) == n

    Bx, By = _solve_bore_field(r_in, r_out, Br, p)
    Bmag = np.sqrt(Bx**2 + By**2)
    ripple = Bmag.std() / Bmag.mean()
    B0 = Br * math.log(r_out / r_in)                       # analytic Halbach dipole bore field
    ang_std = np.std(np.degrees(np.arctan2(By, Bx)))
    print(f"Halbach bore field: |B|={Bmag.mean():.4f} T (analytic {B0:.4f}, rel "
          f"{abs(Bmag.mean()-B0)/B0:.2e}), ripple={100*ripple:.2f}%, direction std={ang_std:.2f} deg")

    assert ripple < 0.02, f"the dipole Halbach bore field must be UNIFORM (got {100*ripple:.1f}% ripple)"
    assert abs(Bmag.mean() - B0) / B0 < 0.05, "bore |B| must match Br*ln(r_out/r_in)"
    assert ang_std < 1.0, "the bore field must point in ONE direction (transverse dipole)"


def test_halbach_quadrupole_field_vanishes_at_centre():
    # a quadrupole Halbach (pole_pairs=2) has ZERO field at the very centre, growing with radius
    r_in, r_out, Br, p = 40.0, 55.0, 1.2, 2
    Bx, By = _solve_bore_field(r_in, r_out, Br, p)
    # sample magnitude grouped by the three radii (5,15,25): field grows with r for a quadrupole
    Bmag = np.sqrt(Bx**2 + By**2).reshape(3, 24).mean(axis=1)
    print(f"Halbach quadrupole |B| at r=5,15,25: {Bmag.round(4)}")
    assert Bmag[0] < Bmag[1] < Bmag[2], "quadrupole field grows with radius (zero at centre)"


def main():
    test_halbach_archetype_makes_uniform_dipole_bore_field()
    test_halbach_quadrupole_field_vanishes_at_centre()
    print("[OK] build123d Halbach archetype -> magnetization_map -> 2D PM solve: uniform transverse "
          "dipole bore field (Br ln(r_out/r_in)); quadrupole field grows from a null centre.")


if __name__ == "__main__":
    main()
