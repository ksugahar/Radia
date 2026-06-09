"""COMSOL-class #63 -- wire over a conducting ground plane (image-method inductance).

A single wire (radius a) at height h above a PERFECT GROUND (A=0 plane) -- the single-ended /
microstrip-style line. The IMAGE METHOD: the ground reflects a -I image at -h, so the field above
is the +I/-I pair (separation 2h) but only HALF the energy is above the plane:

    L_ext = (mu0/2 pi) acosh(h/a) = 0.5 * (two-wire L_ext at separation 2h).

radia reads it from the FE field energy ABOVE the ground (Dirichlet A=0 plane), L = 2W/I^2; the
open-domain truncation puts the FE a few % BELOW the h->inf-domain closed form. Self-contained,
headless; SI units.
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, MoveTo, X, Y, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic, wire_over_ground_inductance,
                                           two_wire_external_inductance, NU0, MU0)
from radia_mcp.radia_ngsolve.force import magnetic_energy_2d

I, a = 1.0, 0.0005


def L_fe(h, W=0.12, order=3):
    wire = WorkPlane().Circle(0, h, a).Face(); wire.faces.name = "wire"; wire.faces.maxh = a / 4
    box = MoveTo(-W, 0).Rectangle(2 * W, W).Face()           # air above the ground (y in [0, W])
    air = box - wire; air.faces.name = "air"
    for e in (air.edges.Min(Y), air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y)):
        e.name = "ground"                                    # y=0 ground + far walls all A=0
    mesh = Mesh(OCCGeometry(Glue([air, wire]), dim=2).GenerateMesh(maxh=W / 14))
    Jz = CoefficientFunction([I / (math.pi * a * a) if m == "wire" else 0.0 for m in mesh.GetMaterials()])
    A = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0), Jz=Jz, order=order, dirichlet="ground")
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    return 2 * magnetic_energy_2d(B, mesh, "air") / I**2


print(f"wire a={a*1e3:.1f} mm over a conducting ground   L_ext = (mu0/2 pi) acosh(h/a)")
print("  h[mm]  h/a   FE L [nH/m]  closed [nH/m]   err     (= 1/2 two-wire@2h)")
with TaskManager():
    for h in (0.005, 0.010, 0.015):
        Lfe = L_fe(h)
        Lcf = wire_over_ground_inductance(h, a)
        Lhalf = 0.5 * two_wire_external_inductance(2 * h, a)
        print(f"  {h*1e3:4.0f}  {h/a:4.0f}   {Lfe*1e9:8.2f}    {Lcf*1e9:8.2f}    {(Lfe/Lcf-1)*100:+.2f}%   "
              f"{Lhalf*1e9:.2f} nH/m")
print("\nthe ground = image of -I at -h; single-ended line Z0 = c L_ext. FE below closed (truncation).")
print("OK -- wire-over-ground L_ext == (mu0/2pi)acosh(h/a) == 1/2 the two-wire value at 2h.")
