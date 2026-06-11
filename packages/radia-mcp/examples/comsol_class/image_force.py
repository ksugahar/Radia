"""COMSOL-class #33 -- METHOD OF IMAGES: current wire above an iron plane.

A classic 2D-EM technique. A wire (current I) at height h above a high-mu iron half-plane
(mu->inf) is attracted to the iron; the image is a SAME-sign current at -h (field lines enter
the iron perpendicularly), so the force per unit length is (two parallel same-direction
currents at separation 2h)

    F = mu0 I^2 / (4 pi h)   [N/m]   (toward the iron).

radia: wire + large high-mu iron block; the force on the wire is the Lorentz body force
int J x B over it (its symmetric self-field nets to zero; the iron's image field gives F).
Self-contained, headless, tool-independent (image-method analytic). With a finite iron block
the FE force is ~4 % below the infinite-half-plane ideal -> exact as the block grows.
"""
import os, sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, MoveTo, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic, reluctivity,
                                           image_force_wire_iron, MU0)
from radia_mcp.radia_ngsolve.force import lorentz_force_2d

a, h, I, mur = 0.002, 0.02, 100.0, 5000.0
IW, ID, Rbox = 1.0, 0.5, 1.0                # iron block ~25h, far box


def build():
    wire = WorkPlane().Circle(0, h, a).Face(); wire.faces.name = "wire"
    iron = MoveTo(-IW/2, -ID).Rectangle(IW, ID).Face(); iron.faces.name = "iron"
    box = WorkPlane().Circle(0, 0, Rbox).Face(); box.edges.name = "outer"
    air = box - wire - iron; air.faces.name = "air"
    wire.faces.maxh = a/3; iron.faces.maxh = ID/8; air.faces.maxh = Rbox/10
    for e in wire.edges:
        e.maxh = a/3
    return Mesh(OCCGeometry(Glue([wire, iron, air]), dim=2).GenerateMesh(maxh=Rbox/8))


mesh = build()
j0 = I / (math.pi * a * a)
Jz = mesh.MaterialCF({"wire": j0}, default=0.0)
with TaskManager():
    A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {"iron": mur}), Jz=Jz,
                                   order=3, dirichlet="outer")
B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
Fx, Fy = lorentz_force_2d(Jz, B, mesh, "wire")
F_cf = image_force_wire_iron(I, h)
err = abs(abs(Fy) - F_cf) / F_cf * 100
print(f"F image = mu0 I^2/(4 pi h) = {F_cf:.5e} N/m  (infinite half-plane)")
print(f"radia: Fy={Fy:.5e} N/m (toward iron), Fx={Fx:.2e} (~0)   |Fy| err = {err:.2f}%")
print("OK" if err < 6.0 else f"OFF ({err:.2f}%)")
