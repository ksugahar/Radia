"""Clean-room P1 surface-triangle gates for readable FEM/BEM coupling."""

import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.scalar_fem3d import (
    laplace_single_layer_far_potential,
    p1_surface_triangle_constant_load,
    p1_surface_triangle_density_moments,
    p1_surface_triangle_geometry,
    p1_surface_triangle_mass,
    p1_surface_triangle_stiffness,
)


def _matvec(M, x):
    return [sum(M[i][j] * x[j] for j in range(len(x))) for i in range(len(x))]


def _quad(x, M):
    y = _matvec(M, x)
    return sum(xi * yi for xi, yi in zip(x, y))


def _dot(a, b):
    return sum(ai * bi for ai, bi in zip(a, b))


def test_reference_surface_triangle_geometry_and_stiffness():
    tri = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    geo = p1_surface_triangle_geometry(tri)

    assert geo["area"] == pytest.approx(0.5)
    assert geo["area_vector"] == pytest.approx((0.0, 0.0, 0.5))
    assert geo["unit_normal"] == pytest.approx((0.0, 0.0, 1.0))
    assert geo["gradients"] == pytest.approx([
        (-1.0, -1.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
    ])

    K = p1_surface_triangle_stiffness(tri)
    expected = [[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]]
    for row, exp in zip(K, expected):
        assert row == pytest.approx(exp)
    assert all(abs(sum(row)) < 1.0e-15 for row in K)


def test_surface_mass_and_constant_load_integrals():
    tri = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    area = p1_surface_triangle_geometry(tri)["area"]
    ones = [1.0, 1.0, 1.0]

    M = p1_surface_triangle_mass(tri, density=3.0)
    assert _quad(ones, M) == pytest.approx(3.0 * area)

    F = p1_surface_triangle_constant_load(tri, source=5.0)
    assert F == pytest.approx([5.0 * area / 3.0] * 3)
    assert sum(F) == pytest.approx(5.0 * area)


def test_surface_density_moments_for_single_layer_bem_gate():
    tri = [(1.0, 0.0, 0.0), (1.0, 2.0, 0.0), (1.0, 0.0, 3.0)]
    area = p1_surface_triangle_geometry(tri)["area"]
    centroid = (1.0, 2.0 / 3.0, 1.0)

    const = p1_surface_triangle_density_moments(tri, [2.0, 2.0, 2.0])
    assert const["total_source"] == pytest.approx(2.0 * area)
    assert const["first_moment"] == pytest.approx(tuple(2.0 * area * c for c in centroid))

    zero_net = p1_surface_triangle_density_moments(tri, [1.0, -1.0, 0.0])
    shifted = [(x + 5.0, y - 4.0, z + 2.0) for x, y, z in tri]
    zero_net_shifted = p1_surface_triangle_density_moments(shifted, [1.0, -1.0, 0.0])
    assert zero_net["total_source"] == pytest.approx(0.0)
    assert zero_net_shifted["first_moment"] == pytest.approx(zero_net["first_moment"])

    far = laplace_single_layer_far_potential(
        (10.0, 0.0, 0.0),
        const["total_source"],
        const["first_moment"],
    )
    assert far["monopole_potential"] == pytest.approx(const["total_source"] / (40.0 * math.pi))
    assert far["dipole_potential"] == pytest.approx(const["first_moment"][0] / (400.0 * math.pi))
    assert far["far_potential"] == pytest.approx(
        far["monopole_potential"] + far["dipole_potential"]
    )


def test_tilted_surface_affine_energy_uses_tangential_gradient():
    tri = [(0.2, -0.1, 0.3), (1.1, 0.4, 0.6), (0.1, 1.2, 1.3)]
    coeff = 2.0
    field = (1.0, 2.0, 3.0)
    nodal = [_dot(field, p) for p in tri]

    geo = p1_surface_triangle_geometry(tri)
    n = geo["unit_normal"]
    projected = tuple(field[i] - _dot(field, n) * n[i] for i in range(3))
    exact = coeff * geo["area"] * _dot(projected, projected)
    K = p1_surface_triangle_stiffness(tri, coeff=coeff)

    assert _quad(nodal, K) == pytest.approx(exact, rel=1.0e-14)


def test_orientation_reversal_permutates_surface_matrix_and_normal():
    tri = [(0.1, 0.0, 0.2), (1.0, 0.2, -0.1), (0.3, 1.4, 0.4)]
    rev = list(reversed(tri))
    geo = p1_surface_triangle_geometry(tri)
    geo_rev = p1_surface_triangle_geometry(rev)
    assert geo_rev["area"] == pytest.approx(geo["area"])
    assert geo_rev["unit_normal"] == pytest.approx(tuple(-v for v in geo["unit_normal"]))

    K1 = p1_surface_triangle_stiffness(tri)
    K2 = p1_surface_triangle_stiffness(rev)
    for i in range(3):
        for j in range(3):
            assert K1[i][j] == pytest.approx(K2[2 - i][2 - j])

    with pytest.raises(ValueError):
        p1_surface_triangle_geometry([(0, 0, 0), (1, 1, 1), (2, 2, 2)])
    with pytest.raises(ValueError):
        p1_surface_triangle_density_moments(tri, [1.0, 2.0])
    with pytest.raises(ValueError):
        laplace_single_layer_far_potential((0.0, 0.0, 0.0), 1.0, (1.0, 0.0, 0.0))
