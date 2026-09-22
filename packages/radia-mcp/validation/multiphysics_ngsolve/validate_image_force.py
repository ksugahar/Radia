r"""Method of images: current wire above an iron plane -- regression test.

A wire (current I) at height h above a high-mu iron half-plane is attracted with
F = mu0 I^2/(4 pi h) (image = same-sign current at -h; image_force_wire_iron). radia gets the
force as the Lorentz body force int J x B over the wire. With a large-but-FINITE iron block the
FE force is ~4 % below the infinite-half-plane ideal -> tol 6 %. Tool-independent (image method).
"""
import math
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, MoveTo, Glue
from radia_mcp.radia_ngsolve.solve import (solve_planar_magnetostatic, reluctivity,
                                           image_force_wire_iron)
from radia_mcp.radia_ngsolve.force import lorentz_force_2d

a, h, I, mur = 0.002, 0.02, 100.0, 5000.0
IW, ID, Rbox = 1.0, 0.5, 1.0


def validate_image_force_wire_above_iron():
    wire = WorkPlane().Circle(0, h, a).Face(); wire.faces.name = "wire"
    iron = MoveTo(-IW/2, -ID).Rectangle(IW, ID).Face(); iron.faces.name = "iron"
    box = WorkPlane().Circle(0, 0, Rbox).Face(); box.edges.name = "outer"
    air = box - wire - iron; air.faces.name = "air"
    wire.faces.maxh = a/3; iron.faces.maxh = ID/8; air.faces.maxh = Rbox/10
    for e in wire.edges:
        e.maxh = a/3
    mesh = Mesh(OCCGeometry(Glue([wire, iron, air]), dim=2).GenerateMesh(maxh=Rbox/8))
    j0 = I / (math.pi * a * a)
    Jz = mesh.MaterialCF({"wire": j0}, default=0.0)
    A = solve_planar_magnetostatic(mesh, reluctivity(mesh, {"iron": mur}), Jz=Jz,
                                   order=3, dirichlet="outer")
    B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
    Fx, Fy = lorentz_force_2d(Jz, B, mesh, "wire")
    F_cf = image_force_wire_iron(I, h)
    rel = abs(abs(Fy) - F_cf) / F_cf
    print(f"[image force] Fy={Fy:.4e} N/m  ideal={F_cf:.4e} N/m  ({100*rel:.2f}%)  Fx={Fx:.2e}")
    assert Fy < 0, "force must be TOWARD the iron (attraction, -y)"
    assert abs(Fx) < 0.02 * F_cf, "Fx should be ~0 by symmetry"
    assert rel < 0.06, f"|Fy| {abs(Fy):.3e} vs image ideal {F_cf:.3e} ({100*rel:.1f}%)"


if __name__ == "__main__":
    validate_image_force_wire_above_iron()
    print("[OK] method-of-images force (wire above iron) validated.")
