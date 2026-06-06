"""COMSOL-class #10 -- ELECTROSTATIC FORCE (AC/DC electrostatic actuator / MEMS).

Attractive force between the plates of a parallel-plate capacitor. To make the check
EXACT (no fringing), the gap is a closed box whose four SIDE walls are natural Neumann
(symmetry) boundaries -> the field is perfectly 1-D (E = V/d), exactly as for infinite
plates. Exact:

    |F| = 1/2 * eps0 * E^2 * A = eps0 * A * V^2 / (2 d^2)

Method: 3D electrostatics (electrostatic3d) + the electrostatic Maxwell-stress eggshell
(force.electrostatic_eggshell_force) with an axis-aligned weight ramp across the gap.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, IfPos, z
from netgen.occ import Box, OCCGeometry, Z
from radia_mcp.radia_ngsolve.electrostatic3d import solve_electrostatic_3d, EPS0
from radia_mcp.radia_ngsolve.force import electrostatic_eggshell_force

V0, d, L = 100.0, 0.001, 0.005          # 100 V, 1 mm gap, 5 mm square plates

box = Box((-L/2, -L/2, 0), (L/2, L/2, d)); box.mat("air")
box.faces.Min(Z).name = "bottom"; box.faces.Max(Z).name = "top"
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.0002))
print(f"mesh: {mesh.ne} elems")

gfV = solve_electrostatic_3d(mesh, EPS0, {"top": V0, "bottom": 0.0}, order=2)
E = -grad(gfV)

# weight ramp g: 0 at z=0.3d -> 1 at z=0.7d (band in the upper gap, gradg in +z)
z0, z1 = 0.3 * d, 0.7 * d
gz = IfPos(z - z0, IfPos(z1 - z, 1.0 / (z1 - z0), 0.0), 0.0)
gradg = CoefficientFunction((0, 0, gz))

Fx, Fy, Fz = electrostatic_eggshell_force(E, mesh, gradg, air_region="air")
F_ref = EPS0 * (L * L) * V0**2 / (2.0 * d**2)        # |F| = eps0 A V^2 / (2 d^2)
err = abs(abs(Fz) - F_ref) / F_ref * 100
print(f"F_z_fem = {Fz:.4e} N   |F|_exact = {F_ref:.4e} N   err = {err:.1f}%")
print("OK" if err < 4 else f"OFF ({err:.1f}%)")
