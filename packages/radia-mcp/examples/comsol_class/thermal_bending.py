"""COMSOL-class #29 -- THERMAL BENDING: through-thickness gradient -> cantilever curvature.

The thermal-bimorph / bimetallic-strip actuator principle. A homogeneous cantilever (clamped
one end) with a LINEAR through-thickness temperature gradient T(y)=dT*(y/h - 1/2) bends
STRESS-FREE into the constant curvature kappa = alpha*dT/h; Euler-Bernoulli tip deflection

    delta_tip = kappa L^2/2 = alpha * dT * L^2 / (2 h).

radia solve_linear_elasticity with a y-linear thermal pre-strain. Self-contained, headless,
tool-independent (beam-theory reference). The BENDING thermo-mechanical case -- distinct from
the AXIAL blocked-bar thermal stress (#20, sigma_xx=-E alpha <dT>).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, TaskManager, y as Y_cf
from netgen.occ import OCCGeometry, MoveTo, X, Y
from radia_mcp.radia_ngsolve.elasticity import (solve_linear_elasticity,
                                                thermal_bending_tip_deflection)

L, h = 0.10, 0.004                       # slender cantilever, L/h = 25
E, nu, alpha, dT = 70e9, 0.33, 23e-6, 50.0      # aluminium; 50 K across the thickness

beam = MoveTo(0, 0).Rectangle(L, h).Face(); beam.faces.name = "beam"
beam.edges.Min(X).name = "clamp"; beam.edges.Max(X).name = "tip"
mesh = Mesh(OCCGeometry(beam, dim=2).GenerateMesh(maxh=h/3))

dT_field = dT * (Y_cf / h - 0.5)         # linear through thickness, zero mean (pure bending)
with TaskManager():
    u = solve_linear_elasticity(mesh, E, nu, dirichlet="clamp",
                                thermal=(alpha, dT_field), plane="stress", order=2)
tip_fe = abs(u(mesh(L, h / 2))[1])
tip_eb = thermal_bending_tip_deflection(alpha, dT, h, L)
err = abs(tip_fe - tip_eb) / tip_eb * 100
print(f"kappa = alpha*dT/h = {alpha*dT/h:.4f} 1/m")
print(f"tip: FE={tip_fe*1e3:.4f} mm   Euler-Bernoulli (alpha dT L^2/2h) = {tip_eb*1e3:.4f} mm   err={err:.2f}%")
print("OK" if err < 4.0 else f"OFF ({err:.2f}%)")
