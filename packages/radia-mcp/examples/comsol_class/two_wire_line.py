"""COMSOL-class #12 -- TWO-WIRE TRANSMISSION LINE capacitance (2D, per unit length).

Two parallel round conductors (radius a, centre separation D) driven differentially
(+-V0/2). Exact line capacitance per unit length:

    C = pi eps0 / arccosh(D / (2a))

Method: 2D electrostatics (scalar_fem2d.solve_electrostatic) with the two wires as
Dirichlet boundaries and a far ground ring; C = 2W/V0^2 (scalar_fem2d.capacitance).
The differential field is a 2D dipole (~1/r^2), so a modest ground ring suffices.
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh
from netgen.geom2d import SplineGeometry
from radia_mcp.radia_ngsolve.scalar_fem2d import solve_electrostatic, capacitance, EPS0

a, D, Rbox, V0 = 0.002, 0.01, 0.15, 1.0          # 2 mm wires, 10 mm apart, ground at 15*D

geo = SplineGeometry()
geo.AddCircle((0, 0), Rbox, leftdomain=1, rightdomain=0, bc="ground")
geo.AddCircle((-D / 2, 0), a, leftdomain=0, rightdomain=1, bc="wire1", maxh=0.0004)
geo.AddCircle((+D / 2, 0), a, leftdomain=0, rightdomain=1, bc="wire2", maxh=0.0004)
geo.SetMaterial(1, "air")
mesh = Mesh(geo.GenerateMesh(maxh=0.01))
print(f"mesh: {mesh.ne} elems")

V = solve_electrostatic(mesh, EPS0, {"wire1": V0 / 2, "wire2": -V0 / 2, "ground": 0.0}, order=2)
C_fem = capacitance(V, mesh, EPS0, V0)                  # F/m
C_ref = math.pi * EPS0 / math.acosh(D / (2 * a))
err = abs(C_fem - C_ref) / C_ref * 100
print(f"C_fem = {C_fem*1e12:.3f} pF/m   C_exact = {C_ref*1e12:.3f} pF/m   err = {err:.1f}%")
print("OK" if err < 3 else f"OFF ({err:.1f}%)")
