r"""Ampere's force law: force per unit length between two parallel wires -- regression test.

Two parallel wires a distance d apart carrying currents I1, I2 attract/repel with
F = mu0 I1 I2/(2 pi d) (two_wire_force_per_length); like currents attract. The analytic primitive of
the image-wire force (image_force_wire_iron = this with the image at 2h). The FE leg gets the force
as the Lorentz body force int J x B over one wire (only the OTHER wire's field nets a force; the
symmetric self-field integrates to zero). Tool-independent (Ampere's force law).
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (two_wire_force_per_length, image_force_wire_iron,
                                           MU0)

a, d, I = 0.002, 0.02, 1.0


def validate_two_wire_force_closed_form():
    # Ampere's force law + the sign/relationship conventions
    F = two_wire_force_per_length(I, I, d)
    assert F == pytest.approx(MU0 * I * I / (2 * math.pi * d))      # mu0 I1 I2/(2 pi d)
    assert F == pytest.approx(1.0e-5)                               # I=1, d=20mm -> 1e-5 N/m
    assert F > 0                                                    # like currents (I1 I2>0) attract
    assert two_wire_force_per_length(I, -I, d) < 0                  # opposite currents repel (sign flip)
    assert two_wire_force_per_length(2 * I, I, d) == pytest.approx(2 * F)   # linear in each current
    assert two_wire_force_per_length(I, I, 2 * d) == pytest.approx(F / 2)   # 1/d falloff
    # image-wire force is the special case I1=I2=I at separation 2h
    h = d / 2
    assert image_force_wire_iron(I, h) == pytest.approx(two_wire_force_per_length(I, I, 2 * h))


def validate_two_wire_force_fe_lorentz():
    # FE leg: two parallel like-current wires -> the Lorentz force on wire 2 reproduces mu0 I^2/(2 pi d),
    # attractive (toward wire 1 at -x), Fy ~ 0 by symmetry. The 3-way's NGSolve member.
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane, Glue
    from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, NU0
    from radia_mcp.radia_ngsolve.force import lorentz_force_2d

    Rbox = 0.30
    w1 = WorkPlane().Circle(-d / 2, 0, a).Face(); w1.faces.name = "w1"; w1.faces.maxh = a / 5
    w2 = WorkPlane().Circle(d / 2, 0, a).Face(); w2.faces.name = "w2"; w2.faces.maxh = a / 5
    box = WorkPlane().Circle(0, 0, Rbox).Face(); box.edges.name = "outer"
    air = box - w1 - w2; air.faces.name = "air"; air.faces.maxh = Rbox / 16
    j0 = I / (math.pi * a * a)
    with TaskManager():
        mesh = Mesh(OCCGeometry(Glue([w1, w2, air]), dim=2).GenerateMesh(maxh=Rbox / 10))
        Jz = mesh.MaterialCF({"w1": j0, "w2": j0}, default=0.0)
        A = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0), Jz=Jz, order=3, dirichlet="outer")
        B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
        Fx, Fy = lorentz_force_2d(Jz, B, mesh, "w2")
    F_cf = two_wire_force_per_length(I, I, d)
    rel = abs(abs(Fx) - F_cf) / F_cf
    # FE Lorentz converges to the Ampere force from below (mesh + finite box): ~1.5% here, ->0 on
    # refinement (3.6/1.5/0.5% at wm=a/3,a/5,a/8). The precise leg is the stored 2D solve (0.001%).
    print(f"[two-wire force] Fx={Fx:.4e} N/m  analytic={F_cf:.4e} N/m  ({100*rel:.2f}%)  Fy={Fy:.2e}")
    assert Fx < 0, "wire 2 must be pulled TOWARD wire 1 (-x): like currents attract"
    assert abs(Fy) < 0.02 * F_cf, "Fy ~ 0 by symmetry"
    assert rel < 0.02, f"|Fx| {abs(Fx):.3e} vs analytic {F_cf:.3e} ({100*rel:.1f}%)"


if __name__ == "__main__":
    validate_two_wire_force_closed_form()
    validate_two_wire_force_fe_lorentz()
    print("[OK] two-wire Ampere force validated (closed form + FE Lorentz).")
