r"""Wire over a conducting ground plane -- external inductance (image method) -- test (#63).

L_ext = (mu0/2 pi) acosh(h/a) = 1/2 the two-wire value at separation 2h (the ground reflects a -I
image at -h; half the energy is above the plane). Closed-form helper (tool-independent) + an FE
field-energy check (L = 2W/I^2 above the ground)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (wire_over_ground_inductance,
                                           two_wire_external_inductance, MU0)


def validate_wire_over_ground_closed_form():
    a = 0.0005
    for h in (0.005, 0.010, 0.020):
        L = wire_over_ground_inductance(h, a)
        assert math.isclose(L, MU0 / (2 * math.pi) * math.acosh(h / a), rel_tol=1e-12)
        # exactly HALF the two-wire external inductance at separation 2h
        assert math.isclose(L, 0.5 * two_wire_external_inductance(2 * h, a), rel_tol=1e-12)
    # monotone increasing with height; h >> a limit -> (mu0/2pi) ln(2h/a)
    assert wire_over_ground_inductance(0.005, a) < wire_over_ground_inductance(0.020, a)
    assert math.isclose(wire_over_ground_inductance(0.05, a),
                        MU0 / (2 * math.pi) * math.log(2 * 0.05 / a), rel_tol=2e-3)


def validate_wire_over_ground_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane, MoveTo, X, Y, Glue
    from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, NU0
    from radia_mcp.radia_ngsolve.force import magnetic_energy_2d

    a, I, W = 0.0005, 1.0, 0.12

    def L_fe(h, order=3):
        wire = WorkPlane().Circle(0, h, a).Face(); wire.faces.name = "wire"; wire.faces.maxh = a / 4
        box = MoveTo(-W, 0).Rectangle(2 * W, W).Face()
        air = box - wire; air.faces.name = "air"
        for e in (air.edges.Min(Y), air.edges.Max(X), air.edges.Min(X), air.edges.Max(Y)):
            e.name = "ground"
        mesh = Mesh(OCCGeometry(Glue([air, wire]), dim=2).GenerateMesh(maxh=W / 14))
        Jz = CoefficientFunction([I / (math.pi * a * a) if m == "wire" else 0.0 for m in mesh.GetMaterials()])
        A = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0), Jz=Jz, order=order, dirichlet="ground")
        B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
        return 2 * magnetic_energy_2d(B, mesh, "air") / I**2

    with TaskManager():
        for h in (0.005, 0.010):
            Lfe = L_fe(h)
            Lcf = wire_over_ground_inductance(h, a)
            assert math.isclose(Lfe, Lcf, rel_tol=5e-2)        # within the open-domain truncation
            assert Lfe < Lcf                                   # finite box undercuts the h->inf closed form


if __name__ == "__main__":
    validate_wire_over_ground_closed_form()
    validate_wire_over_ground_fe()
    print("[OK] wire over ground: L_ext = (mu0/2pi)acosh(h/a) = 1/2 two-wire@2h; FE within truncation.")
