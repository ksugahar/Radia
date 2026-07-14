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
