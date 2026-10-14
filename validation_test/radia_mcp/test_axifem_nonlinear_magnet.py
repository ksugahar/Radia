from __future__ import annotations

import math

import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")


def test_axisymmetric_nonlinear_solver_keeps_permanent_magnet_load() -> None:
    pytest.importorskip("radia.axifem")
    from netgen.occ import Glue, MoveTo, OCCGeometry, X, Y
    from ngsolve import Integrate, Mesh, grad, x

    from radia_mcp.radia_ngsolve.solve import solve_axi_magnetostatic_nonlinear

    mu0 = 4.0e-7 * math.pi
    magnet = MoveTo(0.0, -0.2).Rectangle(0.3, 0.4).Face()
    magnet.faces.name = "magnet"
    magnet.edges.Min(X).name = "axis"
    outer = MoveTo(0.0, -1.0).Rectangle(1.0, 2.0).Face()
    air = outer - magnet
    air.faces.name = "air"
    air.edges.Min(X).name = "axis"
    air.edges.Max(X).name = "outer"
    air.edges.Min(Y).name = "outer"
    air.edges.Max(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([air, magnet]), dim=2).GenerateMesh(maxh=0.12))

    base = 1.0 / (mu0 * mesh.MaterialCF({"magnet": 2.0}, default=1.0))

    def nu_of_b(_field):
        return base

    solution = solve_axi_magnetostatic_nonlinear(
        mesh,
        nu_of_b,
        magnets={"magnet": (3.0e5, 90.0)},
        order=2,
        relax=1.0,
        max_iter=4,
        min_iter=2,
        tol=1.0e-12,
    )
    bz = grad(solution)[0] + solution / x
    weighted_volume = Integrate(x, mesh, definedon=mesh.Materials("magnet"))
    average_bz = Integrate(
        x * bz, mesh, definedon=mesh.Materials("magnet")
    ) / weighted_volume

    assert solution.vec.Norm() > 0.0
    assert average_bz > 0.0
