"""COMSOL-class #55 -- deep-bar (slot-conductor) AC resistance & reactance.

The AC sequel to the slot-leakage inductance (#54): a SOLID conductor filling an iron slot, at
frequency, has its current pushed toward the slot opening by the slot-leakage field (a 1-D skin
effect over the slot DEPTH). Field's closed form:

    k_R = Rac/Rdc = xi (sinh2xi + sin2xi)/(cosh2xi - cos2xi)
    k_X = Lac/Ldc = (3/2xi)(sinh2xi - sin2xi)/(cosh2xi - cos2xi),   xi = h/delta.

This is the DEEP-BAR / double-cage effect: at high slip (high rotor frequency) Rac rises (~sqrt f)
-> high STARTING resistance/torque, while the slot-leakage reactance falls. Modelled with the
current-driven eddy solver on the iron-slot 1-D BVP (top edge Dirichlet yoke / Neumann tooth
walls -- the same BCs as #54, now time-harmonic). Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Y
from radia_mcp.radia_ngsolve.solve import (solve_planar_eddy, deep_bar_skin_ratio,
                                           deep_bar_resistance_factor, deep_bar_reactance_factor,
                                           slot_leakage_permeance, NU0, MU0)

SIGMA = 5.8e7          # copper [S/m]
w, h = 0.004, 0.020    # slot/bar width, depth
I = 1.0


def fe_RL(f, order=4):
    omega = 2 * math.pi * f
    bar = MoveTo(0, 0).Rectangle(w, h).Face(); bar.faces.name = "bar"
    bar.edges.Max(Y).name = "top"                       # yoke flux line (Dirichlet A=0)
    mesh = Mesh(OCCGeometry(bar, dim=2).GenerateMesh(maxh=w / 10))
    gfu = solve_planar_eddy(mesh, CoefficientFunction(NU0), CoefficientFunction(SIGMA), omega,
                            driven_region="bar", total_current=I, order=order, dirichlet="top")
    Vc = gfu.components[1].vec.FV().NumPy()[0]          # complex driving E-field per length
    Z = Vc / I
    return Z.real, Z.imag / omega                       # Rac [ohm/m], Lac [H/m]


Rdc = 1.0 / (SIGMA * w * h)
Ldc = MU0 * slot_leakage_permeance(h, w)               # = mu0 h/3w, the #54 DC slot leakage
print(f"copper deep bar  w={w*1e3:.0f} h={h*1e3:.0f} mm   Rdc={Rdc*1e6:.2f} uohm/m   Ldc={Ldc*1e6:.3f} uH/m")
print("  f[Hz]   xi     Rac/Rdc(FE / k_R)        Lac/Ldc(FE / k_X)")
with TaskManager():
    for f in [5, 20, 50, 100, 200]:
        xi = deep_bar_skin_ratio(h, f, SIGMA)
        Rac, Lac = fe_RL(f)
        kR, kX = deep_bar_resistance_factor(xi), deep_bar_reactance_factor(xi)
        print(f" {f:5.0f}  {xi:5.2f}    {Rac/Rdc:6.3f} / {kR:6.3f} ({(Rac/Rdc/kR-1)*100:+.2f}%)   "
              f"{Lac/Ldc:6.3f} / {kX:6.3f} ({(Lac/Ldc/kX-1)*100:+.2f}%)")

print("\ndeep-bar effect: high f -> Rac up (~sqrt f, starting torque), slot reactance down.")
print("OK -- current-driven eddy FE == Field's k_R/k_X to FE accuracy; xi->0 recovers #54 (mu0 h/3w).")
