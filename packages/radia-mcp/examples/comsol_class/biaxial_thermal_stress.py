"""COMSOL-class #35 -- 2D BIAXIAL thermal stress of a fully-constrained plate (plane stress).

A plate fully constrained in-plane (u=0 on all edges) under a UNIFORM rise dT cannot expand ->
biaxial thermal stress; for PLANE STRESS

    sigma_x = sigma_y = -E alpha dT / (1 - nu)     (compressive),

the 2D constraint factor 1/(1-nu) vs the 1D uniaxial blocked bar -E alpha dT (#20). u=0 is the
exact solution (uniform eigenstrain + all edges clamped), so the FE is exact. radia
solve_linear_elasticity, all edges clamped + uniform thermal. Self-contained, headless,
tool-independent. (PLANE STRAIN -> -E alpha dT/(1-2 nu); radia's plane='strain' thermal load is
currently (1+nu) too small -- flagged for a separate fix -- so this covers plane STRESS.)
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, TaskManager
from netgen.occ import OCCGeometry, MoveTo
from radia_mcp.radia_ngsolve.elasticity import solve_linear_elasticity, stress_2d, biaxial_thermal_stress

L, E, nu, alpha, dT = 0.05, 200e9, 0.30, 12e-6, 50.0

plate = MoveTo(0, 0).Rectangle(L, L).Face(); plate.faces.name = "plate"
for e in plate.edges:
    e.name = "fix"
mesh = Mesh(OCCGeometry(plate, dim=2).GenerateMesh(maxh=L/12))

with TaskManager():
    u = solve_linear_elasticity(mesh, E, nu, dirichlet="fix",
                                thermal=(alpha, dT), plane="stress", order=2)
sig = stress_2d(u, E, nu, "stress", thermal=(alpha, dT))
sxx = sig[0](mesh(L/2, L/2)); syy = sig[3](mesh(L/2, L/2))
cf = biaxial_thermal_stress(E, alpha, dT, nu)
err = abs(sxx - cf) / abs(cf) * 100
print(f"plane stress biaxial: sxx={sxx/1e6:.3f} MPa  syy={syy/1e6:.3f} MPa")
print(f"analytic -E alpha dT/(1-nu) = {cf/1e6:.3f} MPa   err={err:.3f}%   (vs 1D uniaxial {-E*alpha*dT/1e6:.1f} MPa)")
print("OK" if err < 1.0 else f"OFF ({err:.2f}%)")
