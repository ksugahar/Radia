from __future__ import annotations

import numpy as np
import pytest


from radia.bem.electrostatic_p1 import solve_prescribed_potential_p1


def _two_triangles():
    vertices = np.array(
        [
            [-0.5, -0.5, 0.5],
            [0.5, -0.5, 0.5],
            [0.0, 0.5, 0.5],
            [-0.5, -0.5, -0.5],
            [0.5, -0.5, -0.5],
            [0.0, 0.5, -0.5],
        ]
    )
    triangles = np.array([[0, 1, 2], [3, 5, 4]])
    potential = np.array([1.0, 1.0, 1.0, -1.0, -1.0, -1.0])
    return vertices, triangles, potential


def test_prescribed_potential_has_equal_and_opposite_charge():
    vertices, triangles, potential = _two_triangles()
    result = solve_prescribed_potential_p1(vertices, triangles, potential)
    positive = result.charge_on_vertices(np.array([0, 1, 2]))
    negative = result.charge_on_vertices(np.array([3, 4, 5]))
    assert positive > 0.0
    assert negative == pytest.approx(-positive, rel=1e-12, abs=1e-24)
    assert result.total_charge_c == pytest.approx(0.0, abs=1e-24)


@pytest.mark.parametrize(
    ("vertices", "triangles", "potential", "message"),
    [
        (np.zeros((3, 2)), np.array([[0, 1, 2]]), np.zeros(3), "vertices"),
        (np.zeros((3, 3)), np.empty((0, 3), dtype=int), np.zeros(3), "triangles"),
        (np.zeros((3, 3)), np.array([[0, 1, 3]]), np.zeros(3), "out of range"),
    ],
)
def test_prescribed_potential_rejects_invalid_mesh(
    vertices, triangles, potential, message
):
    with pytest.raises(ValueError, match=message):
        solve_prescribed_potential_p1(vertices, triangles, potential)
