"""COMSOL-class #20 -- ELECTRO-THERMO-MECHANICAL chain (3-physics multiphysics).

The full Joule-heating-with-thermal-stress workflow (a busbar bolted at both ends that
heats under current and is then put in compression by its own constrained expansion).
Three NGSolve solvers chained:

    solve_current_flow  ->  joule_heat_source  ->  solve_heat_steady  ->  solve_linear_elasticity

Exact references at each stage: peak rise dT_max = sigma V^2/(8k) (length-independent);
constrained-bar axial stress sigma_xx = -E alpha <dT>, with <dT> = (2/3) dT_max for the
parabolic Joule temperature profile. The ends are the electrodes AND heat sinks AND
mechanical clamps -- one consistent set of boundaries.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, Integrate, dx
from netgen.occ import OCCGeometry, MoveTo, X, Y
from radia_mcp.radia_ngsolve.scalar_fem2d import solve_current_flow
from radia_mcp.radia_ngsolve.multiphysics import solve_heat_steady, joule_heat_source
from radia_mcp.radia_ngsolve.elasticity import (solve_linear_elasticity, stress_2d,
                                                constrained_bar_thermal_stress)

L, h = 0.10, 0.02                          # bar length (x), height (y) [m]
sigma, k, V0 = 1.0e6, 50.0, 0.10           # electrical sigma, thermal k, voltage
E, nu, alpha = 200.0e9, 0.30, 1.2e-5       # Young, Poisson, thermal expansion (steel)

bar = MoveTo(0, 0).Rectangle(L, h).Face(); bar.faces.name = "bar"
bar.edges.Min(X).name = "left"; bar.edges.Max(X).name = "right"
bar.edges.Min(Y).name = "bottom"; bar.edges.Max(Y).name = "top"
mesh = Mesh(OCCGeometry(bar, dim=2).GenerateMesh(maxh=0.003))
print(f"mesh: {mesh.ne} elems")

# 1) ELECTRIC: -div(sigma grad V) = 0
V = solve_current_flow(mesh, sigma, {"left": V0, "right": 0.0}, order=2)
q = joule_heat_source(V, sigma)                                   # sigma |grad V|^2
# 2) THERMAL: -div(k grad T) = q, ends cold
dT = solve_heat_steady(mesh, q, k, conductor="bar", dirichlet="left|right", order=2)
dT_max, dT_mean = dT(mesh(L/2, h/2)), Integrate(dT * dx, mesh) / (L * h)
# 3) MECHANICAL: thermal pre-strain alpha*dT, axis (x) clamped both ends
u = solve_linear_elasticity(mesh, E, nu, dirichletx="left|right", dirichlety="bottom",
                            thermal=(alpha, dT), plane="stress")
sig = stress_2d(u, E, nu, "stress", thermal=(alpha, dT))
sxx = Integrate(sig[0] * dx, mesh) / (L * h)                       # axial stress (avg)

dTmax_ex = sigma * V0 ** 2 / (8.0 * k)
sxx_ex = constrained_bar_thermal_stress(E, alpha, dT_mean)        # -E alpha <dT>
print(f"\n[1->2] peak rise dT_max  : FEM {dT_max:.3f} K   exact sigma V^2/(8k) {dTmax_ex:.3f} K "
      f"({abs(dT_max-dTmax_ex)/dTmax_ex*100:.2f}%)")
print(f"       mean rise <dT>     : {dT_mean:.3f} K   (= 2/3 dT_max = {2/3*dT_max:.3f})")
print(f"[2->3] axial stress sxx   : FEM {sxx/1e6:.3f} MPa  exact -E alpha <dT> {sxx_ex/1e6:.3f} MPa "
      f"({abs(sxx-sxx_ex)/abs(sxx_ex)*100:.2f}%)")
ok = (abs(dT_max-dTmax_ex)/dTmax_ex < 0.01) and (abs(sxx-sxx_ex)/abs(sxx_ex) < 0.01)
print("\nOK (electro-thermo-mechanical chain validated)" if ok else "OFF")
