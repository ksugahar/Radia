"""COMSOL-class #21 -- MAGNETO-MECHANICAL: Lorentz force -> beam deflection.

The magnetic-actuator / loudspeaker / galvanometer principle: a current-carrying beam in
a transverse field B0 feels a distributed Lorentz force f = J x B0 that bends it. The EM
coupling term is the body force f_y = -J_x B0 [N/m^3]; fed to the linear-elasticity solver
it deflects a cantilever. Euler-Bernoulli reference (slender beam, clamped one end):

    delta_tip = w L^4 / (8 E I),   w = J_x B0 h = I B0 [N/m],   I = h^3/12 (per unit depth).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh
from netgen.occ import OCCGeometry, MoveTo, X, Y
from radia_mcp.radia_ngsolve.elasticity import solve_linear_elasticity, cantilever_tip_deflection

L, h = 0.20, 0.008                         # beam length (x), height (y) [m]  (slender L/h=25)
E, nu = 70.0e9, 0.33                        # aluminium
B0, I = 0.5, 100.0                          # transverse field [T], current [A]
Jx = I / h                                  # axial current density [A/m^2] (per unit depth)
fy = -Jx * B0                               # Lorentz body force [N/m^3] (downward)

beam = MoveTo(0, 0).Rectangle(L, h).Face(); beam.faces.name = "beam"
beam.edges.Min(X).name = "left"; beam.edges.Max(X).name = "right"
beam.edges.Min(Y).name = "bottom"; beam.edges.Max(Y).name = "top"
mesh = Mesh(OCCGeometry(beam, dim=2).GenerateMesh(maxh=0.002))
print(f"mesh: {mesh.ne} elems")

# clamp the left end (cantilever); Lorentz body force over the whole beam
u = solve_linear_elasticity(mesh, E, nu, dirichlet="left",
                            body_force=(0.0, fy), plane="stress", order=2)
tip_fem = u(mesh(L, h / 2))[1]                          # transverse tip deflection
w = abs(fy) * h                                         # distributed load [N/m] = I*B0
I_area = h ** 3 / 12.0
tip_eb = -cantilever_tip_deflection(w, L, E, I_area)    # downward (negative)
err = abs(tip_fem - tip_eb) / abs(tip_eb) * 100
print(f"\ndistributed load w = I*B0 = {w:.2f} N/m")
print(f"tip deflection: FEM {tip_fem*1e6:.4f} um   Euler-Bernoulli {tip_eb*1e6:.4f} um  ({err:.2f}%)")
print("OK" if err < 3.0 else f"OFF ({err:.2f}%)")
