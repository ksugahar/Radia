"""Native volume-current sources must agree in direct and NGSolve field paths."""

import os
import sys

import numpy as np
import pytest


sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src"))

import radia as rad


VERTICES = [
    [0.0, 0.0, 0.0],
    [0.1, 0.0, 0.0],
    [0.0, 0.1, 0.0],
    [0.0, 0.0, 0.1],
]
PROBE = [0.035, 0.025, 0.18]


@pytest.fixture(autouse=True)
def _clear_radia_objects():
    rad.UtiDelAll()
    yield
    rad.UtiDelAll()


def _field_for_current(current_density):
    source = rad.ObjTetrahedronCurrent(VERTICES, current_density)
    return np.asarray(rad.Fld(source, "b", PROBE), dtype=float)


def test_tetrahedron_current_is_nonzero_and_linear():
    field = _field_for_current([0.0, 0.0, 1.0e6])
    assert np.all(np.isfinite(field))
    assert np.linalg.norm(field) > 1.0e-8

    rad.UtiDelAll()
    doubled = _field_for_current([0.0, 0.0, 2.0e6])
    np.testing.assert_allclose(doubled, 2.0 * field, rtol=1.0e-11, atol=1.0e-13)


def test_tetrahedron_current_b_matches_curl_of_a():
    source = rad.ObjTetrahedronCurrent(VERTICES, [0.0, 0.0, 1.0e6])
    point = np.asarray(PROBE, dtype=float)
    step = 1.0e-5

    def vector_potential(offset):
        return np.asarray(rad.Fld(source, "a", (point + offset).tolist()), dtype=float)

    curl_a = np.asarray(
        [
            (vector_potential([0.0, step, 0.0])[2]
             - vector_potential([0.0, -step, 0.0])[2]) / (2.0 * step),
            -(vector_potential([step, 0.0, 0.0])[2]
              - vector_potential([-step, 0.0, 0.0])[2]) / (2.0 * step),
            0.0,
        ],
        dtype=float,
    )
    direct_b = np.asarray(rad.Fld(source, "b", point.tolist()), dtype=float)
    np.testing.assert_allclose(direct_b, curl_a, rtol=1.0e-7, atol=1.0e-11)


def test_tetrahedron_current_b_and_h_use_the_public_si_contract():
    """Current-source internal H must receive mu_0 exactly once for B."""
    source = rad.ObjTetrahedronCurrent(VERTICES, [0.0, 0.0, 1.0e6])
    field_b = np.asarray(rad.Fld(source, "b", PROBE), dtype=float)
    field_h = np.asarray(rad.Fld(source, "h", PROBE), dtype=float)
    mu_0 = 4.0e-7 * np.pi
    np.testing.assert_allclose(field_b, mu_0 * field_h, rtol=1.0e-11, atol=1.0e-13)


def test_tetrahedron_current_radiafield_matches_direct_field():
    ng = pytest.importorskip("ngsolve")
    from netgen.occ import Box, OCCGeometry, Pnt

    source = rad.ObjTetrahedronCurrent(VERTICES, [0.0, 0.0, 1.0e6])
    with ng.TaskManager():
        mesh = ng.Mesh(
            OCCGeometry(Box(Pnt(-0.1, -0.1, -0.1), Pnt(0.3, 0.3, 0.3))).GenerateMesh(
                maxh=0.2
            )
        )
        coefficient = rad.RadiaField(source, "h")
        from_coefficient = np.asarray(coefficient(mesh(*PROBE)), dtype=float)
    direct = np.asarray(rad.Fld(source, "h", PROBE), dtype=float)
    np.testing.assert_allclose(from_coefficient, direct, rtol=1.0e-11, atol=1.0e-9)


@pytest.mark.parametrize("current_density", ([1.0, 2.0], [1.0, 2.0, 3.0, 4.0]))
def test_tetrahedron_current_requires_three_components(current_density):
    with pytest.raises(RuntimeError, match="current_density must have 3 elements"):
        rad.ObjTetrahedronCurrent(VERTICES, np.asarray(current_density, dtype=float))
