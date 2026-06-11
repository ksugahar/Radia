"""COMSOL-class #32 -- BIMETALLIC STRIP (two bonded layers, different alpha -> bending).

The thermostat / thermal-switch / thermal-relay classic (Timoshenko 1925). Two bonded layers
of EQUAL thickness and modulus but different expansion (alpha1 bottom, alpha2 top) under a
UNIFORM rise dT curve into kappa = 3(alpha2-alpha1)dT/(2h) (h = total thickness); the cantilever
tip delta = kappa L^2/2. Distinct from #29 (ONE material, through-thickness GRADIENT): here a
uniform dT + a material Delta-alpha drives the bending. radia solve_linear_elasticity with a
region-wise alpha (mesh.MaterialCF). Self-contained, headless, tool-independent (Timoshenko ref).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, TaskManager
from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y
from radia_mcp.radia_ngsolve.elasticity import solve_linear_elasticity, bimetal_tip_deflection

L, h = 0.10, 0.004                     # length, total thickness (each layer h/2)
E, nu, dT = 70e9, 0.33, 50.0
a1, a2 = 10e-6, 23e-6                   # bottom (steel-ish), top (aluminium)

bot = MoveTo(0, 0).Rectangle(L, h/2).Face(); bot.faces.name = "lay1"
top = MoveTo(0, h/2).Rectangle(L, h/2).Face(); top.faces.name = "lay2"
strip = Glue([bot, top]); strip.edges.Min(X).name = "clamp"
mesh = Mesh(OCCGeometry(strip, dim=2).GenerateMesh(maxh=h/4))
alpha = mesh.MaterialCF({"lay1": a1, "lay2": a2}, default=0.0)

with TaskManager():
    u = solve_linear_elasticity(mesh, E, nu, dirichlet="clamp",
                                thermal=(alpha, dT), plane="stress", order=2)
tip_fe = abs(u(mesh(L, h/2))[1])
tip_cf = bimetal_tip_deflection(a1, a2, dT, h, L)
err = abs(tip_fe - tip_cf) / tip_cf * 100
print(f"kappa = 3 dAlpha dT/(2h) = {3*(a2-a1)*dT/(2*h):.4f} 1/m")
print(f"tip: FE={tip_fe*1e3:.4f} mm   Timoshenko (kappa L^2/2) = {tip_cf*1e3:.4f} mm   err={err:.2f}%")
print("OK" if err < 4.0 else f"OFF ({err:.2f}%)")
