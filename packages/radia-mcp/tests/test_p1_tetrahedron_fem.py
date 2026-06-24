"""Clean-room P1 tetrahedron FEM element gates.

These are the tiny local matrices behind scalar 3D Poisson / heat /
current-flow assembly on `.vol` tetrahedra.  The formulas are deliberately
plain so they can be mirrored in MATLAB FEM teaching code.
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.scalar_fem3d import (
    assemble_p1_tet_robin_system,
    p1_tetrahedron_constant_load,
    p1_tetrahedron_geometry,
    p1_tetrahedron_mass,
    p1_tetrahedron_stiffness,
)
from radia_mcp.radia_ngsolve.netgen_vol import NetgenSurfaceTriangle, NetgenTetrahedron


def _matvec(M, x):
    return [sum(M[i][j] * x[j] for j in range(4)) for i in range(4)]


def _quad(x, M):
    y = _matvec(M, x)
    return sum(xi * yi for xi, yi in zip(x, y))


def _matvec_n(M, x):
    return [sum(row[j] * x[j] for j in range(len(x))) for row in M]


def test_reference_tetrahedron_gradients_and_stiffness():
    tet = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    geo = p1_tetrahedron_geometry(tet)
    assert math.isclose(geo["volume"], 1.0 / 6.0, rel_tol=1e-15)
    assert geo["gradients"] == [(-1.0, -1.0, -1.0), (1.0, 0.0, 0.0),
                                (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]

    K = p1_tetrahedron_stiffness(tet)
    expected = [
        [0.5, -1.0 / 6.0, -1.0 / 6.0, -1.0 / 6.0],
        [-1.0 / 6.0, 1.0 / 6.0, 0.0, 0.0],
        [-1.0 / 6.0, 0.0, 1.0 / 6.0, 0.0],
        [-1.0 / 6.0, 0.0, 0.0, 1.0 / 6.0],
    ]
    for row, exp in zip(K, expected):
        for got, want in zip(row, exp):
            assert math.isclose(got, want, abs_tol=1e-15)
    assert all(abs(sum(row)) < 1e-15 for row in K)


def test_affine_function_energy_matches_exact_gradient_integral():
    tet = [(-0.1, 0.2, 0.0), (1.1, 0.3, 0.2), (0.2, 1.4, 0.1), (0.1, 0.4, 1.2)]
    coeff = 4.0
    K = p1_tetrahedron_stiffness(tet, coeff=coeff)
    volume = p1_tetrahedron_geometry(tet)["volume"]
    nodal = [x + 2.0 * y - 0.5 * z for x, y, z in tet]
    exact = coeff * volume * (1.0 ** 2 + 2.0 ** 2 + (-0.5) ** 2)
    assert math.isclose(_quad(nodal, K), exact, rel_tol=1e-14)


def test_mass_and_constant_load_integrals():
    tet = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 4.0)]
    volume = p1_tetrahedron_geometry(tet)["volume"]
    assert math.isclose(volume, 4.0, rel_tol=1e-15)

    M = p1_tetrahedron_mass(tet, density=2.5)
    ones = [1.0, 1.0, 1.0, 1.0]
    assert math.isclose(_quad(ones, M), 2.5 * volume, rel_tol=1e-15)

    F = p1_tetrahedron_constant_load(tet, source=7.0)
    assert all(math.isclose(v, 7.0 * volume / 4.0, rel_tol=1e-15) for v in F)
    assert math.isclose(sum(F), 7.0 * volume, rel_tol=1e-15)


def test_orientation_reversal_and_degenerate_guard():
    tet = [(0.1, 0.0, 0.2), (1.0, 0.2, 0.0), (0.3, 1.4, 0.1), (0.2, 0.4, 1.5)]
    K1 = p1_tetrahedron_stiffness(tet)
    K2 = p1_tetrahedron_stiffness(list(reversed(tet)))
    for i in range(4):
        for j in range(4):
            assert math.isclose(K1[i][j], K2[3 - i][3 - j], rel_tol=1e-14, abs_tol=1e-14)

    with pytest.raises(ValueError):
        p1_tetrahedron_geometry([(0, 0, 0), (1, 0, 0), (0, 1, 0), (2, 2, 0)])
    with pytest.raises(ValueError):
        p1_tetrahedron_geometry([(0, 0, 0), (1, 0, 0), (0, 1, 0)])


def test_p1_tet_robin_assembly_constant_solution_gate():
    points = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    tetrahedra = [NetgenTetrahedron(1, (1, 2, 3, 4))]
    surface_triangles = [
        NetgenSurfaceTriangle(1, 1, 1, 0, (1, 2, 3)),
        NetgenSurfaceTriangle(1, 2, 1, 0, (1, 4, 2)),
        NetgenSurfaceTriangle(1, 2, 1, 0, (2, 4, 3)),
        NetgenSurfaceTriangle(1, 2, 1, 0, (3, 4, 1)),
    ]
    robin = {1: 1.0, 2: 2.0}
    system = assemble_p1_tet_robin_system(
        points, tetrahedra, surface_triangles,
        volume_coeff={1: 1.0},
        source={1: 0.0},
        robin_coeff=robin,
        boundary_flux=robin,
    )
    ones = [1.0] * 4
    residual = [a - b for a, b in zip(_matvec_n(system["matrix"], ones), system["rhs"])]
    expected_base_area = 0.5
    expected_side_area = 1.0 + 0.5 * math.sqrt(3.0)
    expected_robin_weight = expected_base_area + 2.0 * expected_side_area

    assert max(abs(value) for value in residual) < 1.0e-15
    assert system["volume_by_material"][1] == pytest.approx(1.0 / 6.0)
    assert system["boundary_area_by_number"][1] == pytest.approx(expected_base_area)
    assert system["boundary_area_by_number"][2] == pytest.approx(expected_side_area)
    assert system["robin_area_weight"] == pytest.approx(expected_robin_weight)
    assert system["flux_area_weight"] == pytest.approx(expected_robin_weight)
    assert _quad(ones, system["matrix"]) == pytest.approx(expected_robin_weight)
    assert sum(system["rhs"]) == pytest.approx(expected_robin_weight)
    assert len(system["matrix_triplets"]) == 14


def test_p1_tet_robin_assembly_source_and_neumann_limit():
    points = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 4.0)]
    tetrahedra = [(1, 2, 3, 4)]
    surface_triangles = [(1, 2, 3), (1, 4, 2), (2, 4, 3), (3, 4, 1)]
    source = 7.0
    system = assemble_p1_tet_robin_system(points, tetrahedra, surface_triangles, source=source)
    volume = p1_tetrahedron_geometry(points)["volume"]

    assert sum(system["rhs"]) == pytest.approx(source * volume)
    assert max(abs(sum(row)) for row in system["matrix"]) < 1.0e-15
    assert system["volume_by_material"][0] == pytest.approx(volume)
    assert system["robin_area_weight"] == 0.0
    assert system["flux_area_weight"] == 0.0

    with pytest.raises(ValueError):
        assemble_p1_tet_robin_system(points, [(1, 2, 3, 5)], [])
