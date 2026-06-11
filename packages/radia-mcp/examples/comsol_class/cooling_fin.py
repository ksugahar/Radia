"""COMSOL-class #47 -- straight COOLING FIN (extended surface), live-3-way.

A heat-sink fin: a thin fin (length L, thickness t) at base rise dT_base, convecting (film h)
from both long faces, adiabatic tip. Conduction along the fin + distributed surface convection
give the 1-D fin equation, with

    dT(x) = dT_base cosh(m(L-x))/cosh(mL),   m = sqrt(2h/(k t)),
    fin efficiency  eta = tanh(mL)/(mL)  = actual heat / (heat if all at dT_base).

radia uses the Robin heat solver with a BASE Dirichlet (`solve_heat_steady_robin`, base
temperature + convective faces). Validated 3-way: NGSolve == this closed form == the reference
commercial FE code (LiveLink, ~0.14 %, internal). The conduction-convection / heat-sink design
workhorse; generalises the convective slab (#44, no conduction gradient along it). Headless.
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, Integrate, ds, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, X, Y
from radia_mcp.radia_ngsolve.multiphysics import (solve_heat_steady_robin, fin_parameter,
                                                  fin_tip_temperature_rise, fin_efficiency)

L, t, k, dTb = 0.05, 0.002, 200.0, 100.0      # fin length, thickness, k, base rise

fin = WorkPlane().MoveTo(0, -t/2).Rectangle(L, t).Face(); fin.faces.name = "fin"
fin.edges.Min(X).name = "base"
fin.edges.Max(Y).name = "conv"; fin.edges.Min(Y).name = "conv"
fin.edges.Max(X).name = "tip"                 # natural = adiabatic tip
mesh = Mesh(OCCGeometry(fin, dim=2).GenerateMesh(maxh=t/2))

print(f"{'h':>5} {'mL':>6} {'tip dT FE':>10} {'(analytic)':>11} {'eta FE':>8} {'(analytic)':>11}")
for h in (200.0, 500.0, 1500.0):
    m = fin_parameter(h, k, t)
    gT = solve_heat_steady_robin(mesh, CoefficientFunction(0.0), k, "fin", "conv", h,
                                 order=3, dirichlet="base", dirichlet_value=dTb)
    tip_fe = gT(mesh(L * 0.999, 0.0))
    Qconv = Integrate(h * gT * ds(definedon=mesh.Boundaries("conv")), mesh)
    eta_fe = Qconv / (h * (2 * L) * dTb)
    print(f"{h:5.0f} {m*L:6.3f} {tip_fe:10.4f} {fin_tip_temperature_rise(dTb,m,L):11.4f} "
          f"{eta_fe:8.4f} {fin_efficiency(m,L):11.4f}")
print("OK -- fin dT(x)=dT_base cosh(m(L-x))/cosh(mL); eta=tanh(mL)/(mL) (->1 short/fat, ->0 long/thin).")
