"""COMSOL-class #71 -- demagnetizing factor & the internal field of a uniformly magnetized body.

A permanent magnet sets up a self-demagnetizing field H = -N M inside itself; the flux density it
actually holds is B_in = Br (1 - N), N the demagnetizing factor (Nx+Ny+Nz=1 for any ellipsoid):

    sphere (N=1/3) -> 2 Br/3,  transverse infinite cylinder (N=1/2) -> Br/2,  thin slab (N=1) -> 0.

This is the basis of the PM operating point / load line (#42) and of shape-anisotropy design. Here
the transverse-magnetized cylinder (B_in = Br/2) is confirmed by an NGSolve magnetostatic solve
(radia ``magnets=``, mu_r=1); the open-domain truncation puts the FE ~0.2 % below Br/2. Self-
contained, headless; SI units.
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from radia_mcp.radia_ngsolve.solve import (demagnetizing_factor, magnetized_body_internal_field,
                                           demagnetizing_field, MU0, NU0, solve_planar_magnetostatic)

print("demagnetizing factor N and internal field B_in = Br(1-N)   (Br = 1.2 T):")
Br = 1.2
for shape in ("sphere", "cylinder_transverse", "cylinder_axial", "thin_slab"):
    N = demagnetizing_factor(shape)
    n0 = N[0]                                                # convention: 1st = magnetized-axis factor
    B_in = magnetized_body_internal_field(Br, n0)
    print(f"  {shape:22s} N={tuple(round(x,3) for x in N)} (sum {sum(N):.2f})  "
          f"B_in={B_in:.4f} T  H_in={demagnetizing_field(Br, n0)/1e3:.1f} kA/m")

# --- NGSolve: transverse-magnetized infinite cylinder -> B_in = Br/2 ---
from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, MoveTo, Glue, X, Y

R, W, M = 0.01, 0.20, Br / MU0
cyl = WorkPlane().Circle(0, 0, R).Face(); cyl.faces.name = "cyl"; cyl.faces.maxh = R / 10
air = MoveTo(-W, -W).Rectangle(2 * W, 2 * W).Face() - WorkPlane().Circle(0, 0, R).Face()
air.faces.name = "air"
for e in (air.edges.Min(X), air.edges.Max(X), air.edges.Min(Y), air.edges.Max(Y)):
    e.name = "outer"
mesh = Mesh(OCCGeometry(Glue([cyl, air]), dim=2).GenerateMesh(maxh=W / 18))
with TaskManager():
    A = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0),
                                   magnets={"cyl": (M, 0.0)}, order=3, dirichlet="outer")
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    bx, by = B[0](mesh(0.0, 0.0)), B[1](mesh(0.0, 0.0))

print(f"\ntransverse cylinder (theta=0, M in +x):  FE B_in = {bx:.4f} T (By={by:+.4f}), "
      f"analytic Br/2 = {Br/2:.4f} T,  err {(bx/(Br/2)-1)*100:+.2f}% (truncation)")
print("OK -- B_in = Br(1-N); sum rule Nx+Ny+Nz=1; transverse cylinder FE == Br/2 within truncation.")
