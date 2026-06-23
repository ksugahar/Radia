r"""Demagnetizing factor & uniformly-magnetized-body internal field (#71) -- test.

N tuples sum to 1 (ellipsoid rule); B_in = Br(1-N) gives 2Br/3 (sphere) / Br/2 (transverse cylinder)
/ 0 (thin slab); H_in = -N Br/mu0. The transverse-cylinder B_in=Br/2 is confirmed by an NGSolve
solve within the open-domain truncation. Tool-independent."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (demagnetizing_factor, magnetized_body_internal_field,
                                           demagnetizing_field, MU0)


def validate_factors_sum_to_one_and_values():
    for shape in ("sphere", "cylinder_transverse", "cylinder_axial", "thin_slab"):
        N = demagnetizing_factor(shape)
        assert math.isclose(sum(N), 1.0, rel_tol=1e-12)
    assert demagnetizing_factor("sphere") == (1 / 3, 1 / 3, 1 / 3)
    assert demagnetizing_factor("cylinder_transverse") == (0.5, 0.5, 0.0)
    assert demagnetizing_factor("thin_slab") == (1.0, 0.0, 0.0)
    with pytest.raises(ValueError):
        demagnetizing_factor("banana")


def validate_internal_field_and_demag_field():
    Br = 1.2
    # B_in = Br(1-N): sphere 2/3, transverse cyl 1/2, slab 0
    assert math.isclose(magnetized_body_internal_field(Br, 1 / 3), 2 * Br / 3, rel_tol=1e-12)
    assert math.isclose(magnetized_body_internal_field(Br, 0.5), Br / 2, rel_tol=1e-12)
    assert math.isclose(magnetized_body_internal_field(Br, 1.0), 0.0, abs_tol=1e-15)
    # axial cylinder (N=0) holds the full remanence
    assert math.isclose(magnetized_body_internal_field(Br, 0.0), Br, rel_tol=1e-12)
    # demag field H_in = -N Br/mu0 (negative, opposes M)
    assert math.isclose(demagnetizing_field(Br, 0.5), -0.5 * Br / MU0, rel_tol=1e-12)
    assert demagnetizing_field(Br, 1 / 3) < 0.0


def validate_transverse_cylinder_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane, MoveTo, Glue, X, Y
    from radia_mcp.radia_ngsolve.solve import NU0, solve_planar_magnetostatic

    Br, R, W = 1.2, 0.01, 0.20
    cyl = WorkPlane().Circle(0, 0, R).Face(); cyl.faces.name = "cyl"; cyl.faces.maxh = R / 10
    air = MoveTo(-W, -W).Rectangle(2 * W, 2 * W).Face() - WorkPlane().Circle(0, 0, R).Face()
    air.faces.name = "air"
    for e in (air.edges.Min(X), air.edges.Max(X), air.edges.Min(Y), air.edges.Max(Y)):
        e.name = "outer"
    mesh = Mesh(OCCGeometry(Glue([cyl, air]), dim=2).GenerateMesh(maxh=W / 16))
    with TaskManager():
        A = solve_planar_magnetostatic(mesh, CoefficientFunction(NU0),
                                       magnets={"cyl": (Br / MU0, 0.0)}, order=3, dirichlet="outer")
        B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
        bx, by = B[0](mesh(0.0, 0.0)), B[1](mesh(0.0, 0.0))
    assert math.isclose(bx, Br / 2, rel_tol=1.5e-2)        # transverse cylinder N=1/2 -> Br/2
    assert bx < Br / 2                                      # finite box undercuts (truncation)
    assert abs(by) < 1e-2 * (Br / 2)                       # internal field is purely along M (x)


if __name__ == "__main__":
    validate_factors_sum_to_one_and_values()
    validate_internal_field_and_demag_field()
    validate_transverse_cylinder_fe()
    print("[OK] demag factor: sum=1; B_in=Br(1-N) (sphere 2/3, cyl 1/2, slab 0); FE cyl==Br/2.")
