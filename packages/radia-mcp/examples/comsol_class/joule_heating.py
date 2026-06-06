"""COMSOL-class #19 -- ELECTRO-THERMAL (Joule heating) multiphysics.

The canonical two-physics coupling (COMSOL "Joule Heating" interface): current through a
resistive bar dissipates ohmic heat, which sets the temperature. Here we actually CHAIN
the two NGSolve solvers (not impose a heat source by hand):

    solve_current_flow  ->  joule_heat_source (q = sigma|grad V|^2)  ->  solve_heat_steady

A uniform bar (length L along x, voltage V across it, both ends held cold, sides
insulated) has uniform q = sigma(V/L)^2 and the exact parabolic temperature RISE

    dT(x) = (q/2k) x (L - x),   max dT = sigma V^2 / (8k)  at the centre  (length-independent).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, Integrate, dx, grad, InnerProduct
from netgen.occ import OCCGeometry, MoveTo, X, Y
from radia_mcp.radia_ngsolve.scalar_fem2d import solve_current_flow
from radia_mcp.radia_ngsolve.multiphysics import (solve_heat_steady, joule_heat_source,
                                                  joule_bar_temperature_rise)

L, h = 0.10, 0.02                          # bar length (x), height (y) [m]
sigma, k, V0 = 1.0e6, 50.0, 0.10           # conductivity, thermal cond., applied voltage

bar = MoveTo(0, 0).Rectangle(L, h).Face()
bar.faces.name = "bar"
bar.edges.Min(X).name = "left"; bar.edges.Max(X).name = "right"
bar.edges.Min(Y).name = "bottom"; bar.edges.Max(Y).name = "top"
mesh = Mesh(OCCGeometry(bar, dim=2).GenerateMesh(maxh=0.003))
print(f"mesh: {mesh.ne} elems, materials={mesh.GetMaterials()}, boundaries={mesh.GetBoundaries()}")

# --- physics 1: DC conduction  -div(sigma grad V) = 0, V across the ends ------
V = solve_current_flow(mesh, sigma, {"left": V0, "right": 0.0}, order=2)
q = joule_heat_source(V, sigma)                       # sigma |grad V|^2  [W/m^3]
q_uniform = sigma * (V0 / L) ** 2
q_fem = Integrate(q * dx, mesh) / (L * h)             # area-average (should be uniform)
print(f"\nJoule source q: FEM avg {q_fem:.4e}  analytic sigma(V/L)^2 {q_uniform:.4e} W/m^3")

# --- physics 2: steady heat  -div(k grad T) = q, ends cold (T=0 rise) ----------
dT = solve_heat_steady(mesh, q, k, conductor="bar", dirichlet="left|right", order=2)

print(f"\n x/L     dT_fem[K]   dT_exact[K]   err%")
worst = 0.0
for xl in (0.1, 0.25, 0.5, 0.75, 0.9):
    x = xl * L
    df = dT(mesh(x, h / 2))
    de = joule_bar_temperature_rise(sigma, V0, L, k, x)
    err = abs(df - de) / de * 100 if de else 0.0
    worst = max(worst, err)
    print(f"{xl:5.2f}   {df:9.4f}   {de:9.4f}   {err:5.2f}")
dT_max_fem = dT(mesh(L / 2, h / 2))
dT_max_exact = sigma * V0 ** 2 / (8.0 * k)
print(f"\npeak dT (centre): FEM {dT_max_fem:.4f} K   exact sigma V^2/(8k) {dT_max_exact:.4f} K  "
      f"({abs(dT_max_fem-dT_max_exact)/dT_max_exact*100:.2f}%)")
print("OK" if worst < 2.0 and abs(q_fem - q_uniform) / q_uniform < 0.02 else f"OFF ({worst:.2f}%)")
