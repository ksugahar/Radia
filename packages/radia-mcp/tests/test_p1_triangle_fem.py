"""Clean-room P1 triangle FEM element gates.

These are the tiny local matrices behind scalar Poisson / heat / current-flow
assembly.  They mirror the MATLAB P1 FEM prototype path: constant shape gradients,
exact stiffness, consistent mass, and constant-source load.
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.scalar_fem2d import (
    p1_triangle_constant_load,
    p1_triangle_geometry,
    p1_triangle_mass,
    p1_triangle_stiffness,
)


def _matvec(M, x):
    return [sum(M[i][j] * x[j] for j in range(3)) for i in range(3)]


def _quad(x, M):
    y = _matvec(M, x)
    return sum(xi * yi for xi, yi in zip(x, y))


def test_reference_triangle_gradients_and_stiffness():
    tri = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    geo = p1_triangle_geometry(tri)
    assert math.isclose(geo["area"], 0.5, rel_tol=1e-15)
    assert geo["gradients"] == [(-1.0, -1.0), (1.0, 0.0), (0.0, 1.0)]

    K = p1_triangle_stiffness(tri)
    expected = [[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]]
    for row, exp in zip(K, expected):
        for got, want in zip(row, exp):
            assert math.isclose(got, want, abs_tol=1e-15)
    assert all(abs(sum(row)) < 1e-15 for row in K)      # constants are in the nullspace


def test_affine_function_energy_matches_exact_gradient_integral():
    tri = [(-0.2, 0.1), (1.3, 0.4), (0.2, 1.1)]
    K = p1_triangle_stiffness(tri, coeff=2.0)
    area = p1_triangle_geometry(tri)["area"]
    nodal = [x + 2.0 * y for x, y in tri]
    exact = 2.0 * area * (1.0 ** 2 + 2.0 ** 2)
    assert math.isclose(_quad(nodal, K), exact, rel_tol=1e-14)


def test_mass_and_constant_load_integrals():
    tri = [(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)]
    area = p1_triangle_geometry(tri)["area"]
    M = p1_triangle_mass(tri, density=3.0)
    ones = [1.0, 1.0, 1.0]
    assert math.isclose(_quad(ones, M), 3.0 * area, rel_tol=1e-15)

    F = p1_triangle_constant_load(tri, source=5.0)
    assert all(math.isclose(v, 5.0 * area / 3.0, rel_tol=1e-15) for v in F)
    assert math.isclose(sum(F), 5.0 * area, rel_tol=1e-15)


def test_orientation_reversal_and_degenerate_guard():
    tri = [(0.1, 0.0), (1.0, 0.2), (0.3, 1.4)]
    K1 = p1_triangle_stiffness(tri)
    K2 = p1_triangle_stiffness(list(reversed(tri)))
    # Reversing the local node order permutes the matrix by the same reversal.
    for i in range(3):
        for j in range(3):
            assert math.isclose(K1[i][j], K2[2 - i][2 - j], rel_tol=1e-14, abs_tol=1e-14)

    with pytest.raises(ValueError):
        p1_triangle_geometry([(0, 0), (1, 1), (2, 2)])
