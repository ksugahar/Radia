"""Fast regressions for PEEC panel potential coefficients."""

import numpy as np
import pytest

peec_matrices = pytest.importorskip("radia.peec_matrices")
PEECBuilder = peec_matrices.PEECBuilder

EPSILON_0 = 8.854187817e-12


def _potential_matrix(*panels):
    builder = PEECBuilder()
    for panel in panels:
        builder.add_panel(panel)
    _, _, potential, _ = builder.build(include_star=True)
    return np.asarray(potential)


def test_triangle_self_potential_is_physical():
    side = 0.01
    height = side * np.sqrt(3.0) / 2.0
    triangle = [[0.0, 0.0, 0.0], [side, 0.0, 0.0], [side / 2.0, height, 0.0]]

    self_potential = _potential_matrix(triangle)[0, 0]
    characteristic_size = np.sqrt(0.5 * side * height)
    dimensional_scale = 1.0 / (4.0 * np.pi * EPSILON_0 * characteristic_size)

    assert self_potential > 0.0
    assert 0.5 < self_potential / dimensional_scale < 2.0


def test_near_panel_mutual_potential_is_reciprocal_and_decays():
    triangle = [[0.0, 0.0, 0.0], [0.01, 0.0, 0.0], [0.005, 0.0087, 0.0]]
    mutual = []
    for distance in (0.002, 0.005, 0.01, 0.02, 0.05):
        shifted = [[x, y, z + distance] for x, y, z in triangle]
        potential = _potential_matrix(triangle, shifted)
        np.testing.assert_allclose(potential[0, 1], potential[1, 0], rtol=1.0e-12)
        mutual.append(potential[0, 1])

    assert np.all(np.asarray(mutual) > 0.0)
    assert np.all(np.diff(mutual) < 0.0)


def test_quad_self_potential_matches_triangle_split():
    quad = [
        [0.0, 0.0, 0.0],
        [0.01, 0.0, 0.0],
        [0.01, 0.01, 0.0],
        [0.0, 0.01, 0.0],
    ]
    triangle_1 = [quad[0], quad[1], quad[2]]
    triangle_2 = [quad[0], quad[2], quad[3]]

    quad_self = _potential_matrix(quad)[0, 0]
    split_self = 0.5 * (
        _potential_matrix(triangle_1)[0, 0] + _potential_matrix(triangle_2)[0, 0]
    )

    np.testing.assert_allclose(quad_self, split_self, rtol=1.0e-12)
