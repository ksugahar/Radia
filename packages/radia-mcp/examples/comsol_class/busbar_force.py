"""COMSOL-class #16 -- parallel BUSBAR / conductor Lorentz force (fault-force design).

The textbook force between two long parallel conductors a distance d apart, carrying
I1, I2: per unit length

    F/L = mu0 I1 I2 / (2 pi d)   (attractive if parallel, repulsive if anti-parallel).

This is the force that wrecks busbars during a short circuit. Method: 2D planar A_z
magnetostatic for the two wires, then the Lorentz volume integral on one conductor

    F = int J x B dA   ->   Fx = -int Jz By,  Fy = int Jz Bx     (force.lorentz_force_2d)

using the TOTAL field B (the wire's own self-field gives zero net force). Validated for
both the parallel (attractive) and anti-parallel (repulsive) wiring.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
import math
from ngsolve import Mesh, CoefficientFunction, grad, Integrate, dx
from netgen.geom2d import SplineGeometry
from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, NU0
from radia_mcp.radia_ngsolve.force import lorentz_force_2d, MU0

I1, I2, d, rw, Rbox = 1000.0, 1000.0, 0.04, 0.002, 0.40     # A, m
F_exact = MU0 * I1 * I2 / (2.0 * math.pi * d)               # N/m

geo = SplineGeometry()
geo.AddCircle((0, 0), Rbox, leftdomain=1, rightdomain=0, bc="outer")
geo.AddCircle((-d / 2, 0), rw, leftdomain=2, rightdomain=1); geo.SetMaterial(2, "w1")
geo.AddCircle((+d / 2, 0), rw, leftdomain=3, rightdomain=1); geo.SetMaterial(3, "w2")
geo.SetMaterial(1, "air")
mesh = Mesh(geo.GenerateMesh(maxh=0.02))
# area-normalize so int Jz = I exactly on each discrete wire region
A1 = Integrate(mesh.MaterialCF({"w1": 1.0}, default=0.0) * dx, mesh)
A2 = Integrate(mesh.MaterialCF({"w2": 1.0}, default=0.0) * dx, mesh)
print(f"mesh: {mesh.ne} elems ; F_exact = mu0 I1 I2/(2 pi d) = {F_exact:.4f} N/m")

print("\n case          Fx_fem[N/m]  Fy_fem[N/m]  |F|_exact   err%   sign")
worst = 0.0
for name, s2 in (("parallel (attract)", +1.0), ("anti-parallel (repel)", -1.0)):
    Jz = mesh.MaterialCF({"w1": I1 / A1, "w2": s2 * I2 / A2}, default=0.0)
    gfu = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0), Jz=Jz, order=3)
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    Fx, Fy = lorentz_force_2d(Jz, B, mesh, "w2")           # force on wire 2 (at +d/2)
    err = abs(abs(Fx) - F_exact) / F_exact * 100
    worst = max(worst, err)
    # parallel -> wire2 pulled toward wire1 (-x); anti-parallel -> pushed +x
    sign_ok = "OK" if ((s2 > 0 and Fx < 0) or (s2 < 0 and Fx > 0)) else "BAD"
    print(f" {name:22s} {Fx:+.4f}   {Fy:+.4f}   {F_exact:.4f}   {err:5.2f}  {sign_ok}")
print("\nOK" if worst < 3.0 else f"\nOFF ({worst:.2f}%)")
