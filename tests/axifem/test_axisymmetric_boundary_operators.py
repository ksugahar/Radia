import os
import sys

import numpy as np
from netgen.geom2d import SplineGeometry
from ngsolve import BilinearForm, CoefficientFunction, LinearForm, Mesh, ds, dx, grad, x

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
for _path in (
    os.path.join(_ROOT, "src"),
    os.path.join(_ROOT, "packages", "radia-mcp", "src"),
):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from radia.axifem import H1Henrotte
from radia_mcp.radia_ngsolve.solve import (
    axisymmetric_constraint_residual_contract,
    solve_axi_magnetostatic,
    solve_axi_magnetostatic_dual_boundary_average,
)


def _mesh():
    geometry = SplineGeometry()
    points = [
        geometry.AppendPoint(0.2, 0.0),
        geometry.AppendPoint(1.0, 0.0),
        geometry.AppendPoint(1.0, 1.0),
        geometry.AppendPoint(0.2, 1.0),
    ]
    geometry.Append(["line", points[0], points[1]], bc="bottom")
    geometry.Append(["line", points[1], points[2]], bc="right")
    geometry.Append(["line", points[2], points[3]], bc="top")
    geometry.Append(["line", points[3], points[0]], bc="left")
    return Mesh(geometry.GenerateMesh(maxh=0.2))


def _assembled(fes, mixed=None):
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += (1.0 / x) * (x * grad(u)[0] + u) * (x * grad(v)[0] + v) * dx
    a += x * grad(u)[1] * grad(v)[1] * dx
    f = LinearForm(fes)
    f += x * v * dx
    if mixed is not None:
        c0, c1 = mixed
        a += c0 * x * u.Trace() * v.Trace() * ds("right")
        f += c1 * x * v.Trace() * ds("right")
    a.Assemble()
    f.Assemble()
    return a, f


def test_signed_dof_identification_is_exact_and_reduced_residual_is_small():
    mesh = _mesh()
    fes = H1Henrotte(mesh, order=2, dirichlet="left")
    free_vertices = []
    from ngsolve import NodeId, VERTEX

    for vertex in mesh.vertices:
        dofs = [dof for dof in fes.GetDofNrs(NodeId(VERTEX, vertex.nr)) if dof >= 0]
        if dofs and fes.FreeDofs()[dofs[0]]:
            free_vertices.append(dofs[0])
    constraints = [(free_vertices[0], free_vertices[-1], 1.0)]
    solution = solve_axi_magnetostatic(
        mesh,
        CoefficientFunction(1.0),
        Jr=CoefficientFunction(1.0),
        order=2,
        dirichlet="left",
        dof_constraints=constraints,
    )
    a, f = _assembled(solution.space)
    contract = axisymmetric_constraint_residual_contract(
        solution.space,
        mesh,
        a.mat,
        f.vec,
        solution.vec,
        dof_constraints=constraints,
    )
    assert contract["constraint_count"] == 1
    assert contract["max_identification_abs_error"] < 1e-12
    assert contract["relative_reduced_residual"] < 1e-10


def test_antiperiodic_self_identification_pins_the_shared_trace_vertex():
    mesh = _mesh()
    fes = H1Henrotte(mesh, order=2, dirichlet="left")
    from ngsolve import NodeId, VERTEX

    dof = next(
        candidate
        for vertex in mesh.vertices
        for candidate in fes.GetDofNrs(NodeId(VERTEX, vertex.nr))
        if candidate >= 0 and fes.FreeDofs()[candidate]
    )
    constraint = [(dof, dof, -1.0)]
    solution = solve_axi_magnetostatic(
        mesh,
        CoefficientFunction(1.0),
        Jr=CoefficientFunction(1.0),
        order=2,
        dirichlet="left",
        dof_constraints=constraint,
    )
    assert abs(float(solution.vec[dof])) < 1e-14


def test_mixed_boundary_terms_are_included_in_the_solved_system():
    mesh = _mesh()
    mixed = (3.0, 0.25)
    solution = solve_axi_magnetostatic(
        mesh,
        CoefficientFunction(1.0),
        Jr=CoefficientFunction(1.0),
        order=2,
        dirichlet="left",
        mixed_boundaries={"right": mixed},
    )
    a, f = _assembled(solution.space, mixed)
    contract = axisymmetric_constraint_residual_contract(
        solution.space, mesh, a.mat, f.vec, solution.vec
    )
    assert contract["relative_reduced_residual"] < 1e-10


def test_dual_boundary_average_returns_verifiable_components():
    mesh = _mesh()
    averaged, natural, essential = solve_axi_magnetostatic_dual_boundary_average(
        mesh,
        CoefficientFunction(1.0),
        dual_boundary="right",
        Jr=CoefficientFunction(1.0),
        order=2,
        dirichlet="left",
    )
    average_values = averaged.vec.FV().NumPy()
    expected = 0.5 * (
        natural.vec.FV().NumPy() + essential.vec.FV().NumPy()
    )
    assert np.max(np.abs(average_values - expected)) < 1e-14
    assert np.linalg.norm(natural.vec.FV().NumPy() - essential.vec.FV().NumPy()) > 1e-6
