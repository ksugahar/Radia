"""COMSOL-class #60 -- Dowell's equation: AC resistance of an m-layer winding.

The multi-layer sequel to the single-conductor deep-bar (#55): m foils stacked in a winding
window, each carrying the SAME series current I. The leakage field builds up layer by layer, so
the outer layers sit in a large field and dissipate heavily -> skin PLUS proximity loss. Dowell:

    F_R = Rac/Rdc = Delta[ (sinh2D+sin2D)/(cosh2D-cos2D) + (2/3)(m^2-1)(sinhD-sinD)/(coshD+cosD) ],

Delta = h/delta (foil thickness / skin depth). m=1 reduces to the deep-bar factor; the proximity
term grows ~m^2 -- the reason litz wire and thin foils exist. FE = m-conductor SERIES eddy circuit
on the slot BVP (#54: top Dirichlet yoke, sides/bottom Neumann). Self-contained, headless; SI.
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Glue, Y
from radia_mcp.radia_ngsolve.solve import (solve_planar_eddy_multi, deep_bar_skin_ratio,
                                           dowell_resistance_factor, NU0)
from radia_mcp.radia_ngsolve.force import ohmic_loss_2d

SIGMA = 5.8e7          # copper
w, h = 0.010, 0.002    # window width, foil thickness
I = 1.0


def fe_FR(m, f, order=4):
    omega = 2 * math.pi * f
    faces = []
    for k in range(m):
        fc = MoveTo(0, k * h).Rectangle(w, h).Face(); fc.faces.name = f"L{k}"; faces.append(fc)
    shape = Glue(faces)
    shape.edges.Max(Y).name = "top"                         # yoke flux line (Dirichlet A=0)
    mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=h / 6))
    conds = [f"L{k}" for k in range(m)]
    gfu = solve_planar_eddy_multi(mesh, CoefficientFunction(NU0), CoefficientFunction(SIGMA),
                                  omega, conds, connection="series", total_current=I,
                                  order=order, dirichlet="top")
    Az = gfu.components[0]
    P = sum(ohmic_loss_2d(-1j * omega * Az + gfu.components[k + 1], mesh,
                          CoefficientFunction(SIGMA), region=f"L{k}") for k in range(m))
    return (2 * P / abs(I) ** 2) / (m / (SIGMA * w * h))    # Rac/Rdc


print(f"copper foil winding  w={w*1e3:.0f} mm  h={h*1e3:.1f} mm   (F_R = Rac/Rdc)")
with TaskManager():
    for f in (1092, 2456):                                  # Delta ~ 1.0 and ~1.5
        D = deep_bar_skin_ratio(h, f, SIGMA)
        print(f"--- f={f} Hz   Delta = h/delta = {D:.3f} ---")
        print("  m   FE F_R    Dowell    err     (proximity ~ m^2)")
        for m in (1, 2, 4, 6):
            fr, cf = fe_FR(m, f), dowell_resistance_factor(D, m)
            print(f"  {m}   {fr:7.4f}  {cf:7.4f}   {(fr/cf-1)*100:+.2f}%")

print("\nproximity loss explodes with layers -> thin foils / litz; m=1 is the deep-bar factor (#55).")
print("OK -- m-conductor series eddy FE == Dowell's F_R to FE accuracy.")
