"""COMSOL-class #22 -- MEMS ELECTRO-MECHANICAL chain (electrostatic force -> deflection).

The electric twin of #21 (magneto-mechanical): an applied field exerts a Maxwell-stress
force that deflects a structure.  Here the parallel-plate electrostatic pull (#10) drives
a clamped micro-beam (the movable electrode).  ONE-WAY coupling chained from existing
solvers -- no new physics, just wiring:

    solve_electrostatic  ->  electrostatic_surface_force_2d (P = 1/2 eps0 |E|^2)
                          ->  solve_linear_elasticity (traction)  ->  tip deflection

Two exact references:
  (a) electrostatic pull on the electrode  P = 1/2 eps0 (V0/d)^2   (1-D gap, no fringing
      because the gap's side walls are natural-Neumann symmetry planes).
  (b) cantilever tip deflection under that uniform pressure  delta = w L^4 / (8 E I),
      w = P (per-unit-depth line load), I = h^3/12  (Euler-Bernoulli, slender beam).

Two-way electrostatic pull-in (the deformation changes the gap -> changes the force) is
NONLINEAR and needs a Picard / moving-mesh loop on top; this is the linear small-signal
limit (the validated building block).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad
from netgen.occ import OCCGeometry, MoveTo, X, Y
from radia_mcp.radia_ngsolve.scalar_fem2d import solve_electrostatic
from radia_mcp.radia_ngsolve.force import electrostatic_eggshell_force_2d, EPS0
from radia_mcp.radia_ngsolve.elasticity import (solve_linear_elasticity,
                                                cantilever_tip_deflection)

# --- electrostatic gap (parallel plate, SI metres) -------------------------------
W, d, V0 = 20e-3, 1e-3, 500.0          # 20 mm wide, 1 mm gap, 500 V

face = MoveTo(0, 0).Rectangle(W, d).Face(); face.faces.name = "air"
face.edges.Max(Y).name = "top"; face.edges.Min(Y).name = "bottom"   # sides -> natural Neumann
mesh = Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=d / 6))

V = solve_electrostatic(mesh, EPS0, {"top": V0, "bottom": 0.0}, order=2)
E = CoefficientFunction((-grad(V)[0], -grad(V)[1]))
# weight g = y/d over the WHOLE source-free gap -> gradg = (0, 1/d) CONSTANT, and the
# ramp "edges" coincide with the (mesh-aligned) conductor surfaces.  This avoids the
# band-edge quadrature error you get from a discontinuous IfPos band straddling coarse
# elements (4.9% at maxh=d/4, 1.3% at d/12); a constant gradg integrates exactly.
# (Volume band -- a boundary trace of grad(V) reads ~0 since V=const on the conductor;
#  see the helper docstring + ngsolve_usage("electro_mechanical").)
gradg = CoefficientFunction((0.0, 1.0 / d))
Fx, Fy = electrostatic_eggshell_force_2d(E, mesh, gradg, air_region="air")  # N/m
P_fem = abs(Fy) / W                                        # electrostatic pressure [Pa]
P_exact = 0.5 * EPS0 * (V0 / d) ** 2
err_P = abs(P_fem - P_exact) / P_exact * 100
print(f"[electrostatic] P_fem={P_fem:.5e} Pa  P_exact={P_exact:.5e} Pa  err={err_P:.3f}%")

# --- movable electrode = clamped cantilever, loaded by the electrostatic pressure --
Lb, h, Emod, nu = 10e-3, 0.4e-3, 1.0e8, 0.3      # slender (L/h=25) soft beam
beam = MoveTo(0, 0).Rectangle(Lb, h).Face(); beam.faces.name = "beam"
beam.edges.Min(X).name = "clamp"; beam.edges.Max(Y).name = "loaded"
meshB = Mesh(OCCGeometry(beam, dim=2).GenerateMesh(maxh=h / 3))

u = solve_linear_elasticity(meshB, Emod, nu, dirichlet="clamp",
                            traction={"loaded": (0.0, -P_fem)}, plane="stress", order=2)
tip_fem = abs(u(meshB(Lb, h / 2))[1])
I_area = h ** 3 / 12.0
tip_euler = cantilever_tip_deflection(P_fem, Lb, Emod, I_area)   # w = P (line load, unit depth)
err_d = abs(tip_fem - tip_euler) / tip_euler * 100
print(f"[mechanical]    tip_fem={tip_fem:.5e} m  tip_euler={tip_euler:.5e} m  err={err_d:.2f}%")

ok = err_P < 1.0 and err_d < 4.0
print("OK" if ok else f"OFF (P {err_P:.2f}% / delta {err_d:.2f}%)")
