"""Fast same-process state-isolation contract for the HDiv-VIM production solve.

The C++ Gram, PARDISO, TaskManager, and Python backend registry must survive
alternating solves without a rerun or a fresh Python interpreter.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import Box, Pnt, OCCGeometry  # noqa: E402
from radia.vim import Solve  # noqa: E402


def _solve(mesh, direction, mu_r):
    with ng.TaskManager():
        result = Solve(
            mesh, mu_r=mu_r,
            H_ext=ng.CoefficientFunction(tuple(direction)),
            order=1)
    result["_coefficients"] = result["gfM"].vec.FV().NumPy().copy()
    return result


def _solve_nonlinear(mesh, direction):
    mu0 = 4.0e-7 * np.pi
    linear_bh = np.asarray([
        [0.0, 0.0],
        [1.0e3, mu0 * 80.0 * 1.0e3],
        [1.0e5, mu0 * 80.0 * 1.0e5],
    ])
    with ng.TaskManager():
        result = Solve(
            mesh, bh_table=linear_bh,
            H_ext=ng.CoefficientFunction(tuple(direction)),
            order=1, tol=1.0e-9, nl_maxit=20)
    result["_coefficients"] = result["gfM"].vec.FV().NumPy().copy()
    return result


def test_alternating_linear_solves_are_same_process_deterministic():
    mesh_a = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.7))
    mesh_b = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1.2, 0.8, 1.0))).GenerateMesh(maxh=0.7))

    first = _solve(mesh_a, (0.0, 0.0, 1.0), 100.0)
    _solve(mesh_b, (1.0, 0.0, 0.0), 500.0)
    second = _solve(mesh_a, (0.0, 0.0, 1.0), 100.0)
    _solve(mesh_b, (0.0, 1.0, 0.0), 50.0)
    third = _solve(mesh_a, (0.0, 0.0, 1.0), 100.0)

    for current in (second, third):
        np.testing.assert_allclose(
            current["_coefficients"], first["_coefficients"], rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(current["M_avg"], first["M_avg"], rtol=2e-12, atol=2e-12)
        assert current["iters"] == first["iters"]


def test_alternating_nonlinear_solves_are_bitwise_deterministic():
    mesh_a = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.7))
    mesh_b = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1.2, 0.8, 1.0))).GenerateMesh(maxh=0.7))

    first = _solve_nonlinear(mesh_a, (0.0, 0.0, 1000.0))
    _solve_nonlinear(mesh_b, (800.0, 0.0, 0.0))
    second = _solve_nonlinear(mesh_a, (0.0, 0.0, 1000.0))

    np.testing.assert_array_equal(second["_coefficients"], first["_coefficients"])
    np.testing.assert_array_equal(second["M_avg"], first["M_avg"])
    assert second["iters"] == first["iters"]
    assert second["linear_solver"] == first["linear_solver"] == "energy-newton-cpp"
