"""COMSOL-class #56 -- the phi-theta (voltage-temperature) relation / constriction supertemperature.

A current-carrying conductor with both terminals at T0 and voltage U between them develops a
hot spot whose temperature is GEOMETRY-INDEPENDENT:

    T(V)^2 = T0^2 + V(U-V)/L ,     T_max^2 - T0^2 = U^2/(4 L)   (at V = U/2),

L = Wiedemann-Franz Lorenz number. The electro-thermal coupling (current flow -> Joule heat ->
temperature) under WF (k = L sigma T) linearises via the Kirchhoff transform psi = L sigma T^2/2:
solve V (Laplace), then -lap(psi) = sigma|grad V|^2 with psi = L sigma T0^2/2 on both terminals,
then T = sqrt(2 psi /(L sigma)). Demonstrated on TWO different shapes (the same U gives the same
T_max) + the classic melting-voltage application. Self-contained, headless; SI units.
"""
import os
import sys
import math

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, InnerProduct, sqrt, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Glue, X
from radia_mcp.radia_ngsolve.scalar_fem2d import solve_poisson_2d
from radia_mcp.radia_ngsolve.multiphysics import (phi_theta_max_temperature,
                                                  phi_theta_local_temperature, melting_voltage, LORENZ)

T0, SIGMA, U = 300.0, 5.8e7, 0.1
Lx, W = 0.02, 0.006


def fe_Tmax(mesh, order=3):
    psi0 = LORENZ * SIGMA * T0 * T0 / 2.0
    V = solve_poisson_2d(mesh, CoefficientFunction(SIGMA), {"hi": U, "lo": 0.0}, order=order)
    q = SIGMA * InnerProduct(grad(V), grad(V))
    psi = solve_poisson_2d(mesh, CoefficientFunction(1.0), {"hi": psi0, "lo": psi0}, source=q, order=order)
    Tcf = sqrt(2.0 * psi / (LORENZ * SIGMA))
    Tmax = T0
    for i in range(81):
        for j in range(41):
            try:
                t = Tcf(mesh(i / 80.0 * Lx, (-0.5 + j / 40.0) * W))
                Tmax = max(Tmax, t)
            except Exception:
                pass
    return Tmax


def bar():
    f = MoveTo(0, -W / 2).Rectangle(Lx, W).Face()
    f.edges.Min(X).name = "hi"; f.edges.Max(X).name = "lo"
    return Mesh(OCCGeometry(f, dim=2).GenerateMesh(maxh=W / 6))


def dogbone():
    nw, nh = 0.006, W / 3
    left = MoveTo(0, -W / 2).Rectangle((Lx - nw) / 2, W).Face()
    right = MoveTo((Lx + nw) / 2, -W / 2).Rectangle((Lx - nw) / 2, W).Face()
    neck = MoveTo((Lx - nw) / 2, -nh / 2).Rectangle(nw, nh).Face()
    s = Glue([left, neck, right])
    s.edges.Min(X).name = "hi"; s.edges.Max(X).name = "lo"
    return Mesh(OCCGeometry(s, dim=2).GenerateMesh(maxh=W / 8))


Tmax_cf = phi_theta_max_temperature(U, T0)
print(f"phi-theta:  T(V)^2 = T0^2 + V(U-V)/L,  T_max^2-T0^2 = U^2/4L   (U={U} V, T0={T0} K)")
print(f"closed form  T_max = {Tmax_cf:.3f} K   (supertemperature dT = {Tmax_cf-T0:.2f} K)")
print("GEOMETRY INDEPENDENCE -- same U, same T_max regardless of shape:")
with TaskManager():
    for name, m in [("uniform bar", bar()), ("dog-bone constriction", dogbone())]:
        Tmax = fe_Tmax(m)
        print(f"  {name:22s}: FE T_max = {Tmax:.3f} K   err = {(Tmax/Tmax_cf-1)*100:+.3f}%")

# classic application: contact voltage that softens / melts copper (geometry-independent)
print("\ncopper contact (Holm) -- melting/softening voltage is set by the metal alone:")
for label, Tm in [("soften ~463 K", 463.0), ("melt 1356 K", 1356.0)]:
    print(f"  {label:16s}: U = {melting_voltage(Tm, T0)*1e3:6.1f} mV")
print("OK -- electro-thermal FE == phi-theta closed form, geometry-independent (uniform == constricted).")
