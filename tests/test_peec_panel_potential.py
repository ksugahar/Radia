"""Fast regressions for PEEC panel potential coefficients."""

import numpy as np
import pytest

peec_matrices = pytest.importorskip("radia.peec_matrices")
PEECBuilder = peec_matrices.PEECBuilder

def _potential_matrix(*panels):
    builder = PEECBuilder()
    for panel in panels:
        builder.add_panel(panel)
    _, _, potential, _ = builder.build(include_star=True)
    return np.asarray(potential)


def _scaled(panel, factor):
    return [[factor * coordinate for coordinate in point] for point in panel]


def test_triangle_self_potential_obeys_inverse_length_scaling():
    side = 0.01
    height = side * np.sqrt(3.0) / 2.0
    triangle = [[0.0, 0.0, 0.0], [side, 0.0, 0.0], [side / 2.0, height, 0.0]]

    self_potential = _potential_matrix(triangle)[0, 0]
    doubled_self_potential = _potential_matrix(_scaled(triangle, 2.0))[0, 0]

    assert np.isfinite(self_potential) and self_potential > 0.0
    assert self_potential / doubled_self_potential == pytest.approx(2.0, rel=5.0e-6)


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


def test_quad_self_potential_obeys_inverse_length_scaling():
    quad = [
        [0.0, 0.0, 0.0],
        [0.01, 0.0, 0.0],
        [0.01, 0.01, 0.0],
        [0.0, 0.01, 0.0],
    ]
    quad_self = _potential_matrix(quad)[0, 0]
    doubled_quad_self = _potential_matrix(_scaled(quad, 2.0))[0, 0]

    assert np.isfinite(quad_self) and quad_self > 0.0
    assert quad_self / doubled_quad_self == pytest.approx(2.0, rel=5.0e-6)
