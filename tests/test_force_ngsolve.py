"""radia.force_ngsolve must agree with the sampled radia.force kernel.

The two routes exist because some callers hold an NGSolve CoefficientFunction
and some hold quadrature samples.  They are only worth having if they give the
same number, so this locks the symbolic boundary integral against an explicit
Gauss-Legendre evaluation of the same analytic field on the same patch.
"""

import math

import numpy as np
import pytest

ngsolve = pytest.importorskip("ngsolve")
occ = pytest.importorskip("netgen.occ")

from ngsolve import CoefficientFunction, Mesh, TaskManager, x, y, z  # noqa: E402

from radia.force import (  # noqa: E402
    MU0,
    integrate_maxwell_surface_force,
    integrate_time_average_maxwell_surface_force,
)
from radia.force_ngsolve import (  # noqa: E402
    maxwell_surface_force,
    time_average_maxwell_surface_force,
)


def _unit_box_with_named_right_face():
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    for face in box.faces:
        face.name = "other"
    box.faces.Max(occ.X).name = "right"
    return Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=0.34))


def _analytic_field(points):
    """B(x, y, z) = (x, -y, 0) -- divergence free and curl free, so the stress
    integrand is smooth and the two routes must match to quadrature accuracy."""
    points = np.asarray(points, dtype=float)
    field = np.zeros_like(points)
    field[..., 0] = points[..., 0]
    field[..., 1] = -points[..., 1]
    return field


def _right_face_quadrature(order=12):
    """Gauss-Legendre samples on the x = 1 face of the unit box, normal +x."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    coords = 0.5 * (nodes + 1.0)
    span = 0.5 * weights
    points, normals, areas = [], [], []
    for cy, wy in zip(coords, span):
        for cz, wz in zip(coords, span):
            points.append([1.0, cy, cz])
            normals.append([1.0, 0.0, 0.0])
            areas.append(wy * wz)
    return np.array(points), np.array(normals), np.array(areas)


def test_static_surface_force_matches_the_sampled_kernel():
    mesh = _unit_box_with_named_right_face()
    field = CoefficientFunction((x, -y, 0.0))
    with TaskManager():
        symbolic = maxwell_surface_force(field, mesh, mesh.Boundaries("right"))

    points, normals, areas = _right_face_quadrature()
    sampled = integrate_maxwell_surface_force(_analytic_field(points), normals, areas)

    assert np.allclose(symbolic, sampled, rtol=1e-8, atol=1e-8 * np.abs(sampled).max())


def test_orientation_flag_negates_the_traction_integral():
    mesh = _unit_box_with_named_right_face()
    field = CoefficientFunction((x, -y, 0.0))
    with TaskManager():
        outward = maxwell_surface_force(field, mesh, mesh.Boundaries("right"))
        inward = maxwell_surface_force(
            field, mesh, mesh.Boundaries("right"), normal_points_out_of_body=False
        )
    assert np.allclose(inward, [-value for value in outward], rtol=1e-12, atol=1e-12)


def test_peak_phasor_force_is_half_the_static_force_of_the_same_real_field():
    mesh = _unit_box_with_named_right_face()
    field = CoefficientFunction((x, -y, 0.0))
    phasor = CoefficientFunction((x + 0j, -y + 0j, 0.0 + 0j))
    with TaskManager():
        static = maxwell_surface_force(field, mesh, mesh.Boundaries("right"))
        averaged = time_average_maxwell_surface_force(
            phasor, mesh, mesh.Boundaries("right"), amplitude="peak"
        )
        rms = time_average_maxwell_surface_force(
            phasor, mesh, mesh.Boundaries("right"), amplitude="rms"
        )
    assert np.allclose(averaged, [0.5 * value for value in static], rtol=1e-10, atol=1e-12)
    assert np.allclose(rms, static, rtol=1e-10, atol=1e-12)


def test_phasor_surface_force_matches_the_sampled_phasor_kernel():
    mesh = _unit_box_with_named_right_face()
    phasor = CoefficientFunction((x + 0.5j * y, -y + 0j, 0.25j * z))
    with TaskManager():
        symbolic = time_average_maxwell_surface_force(
            phasor, mesh, mesh.Boundaries("right"), amplitude="peak"
        )

    points, normals, areas = _right_face_quadrature()
    sampled_field = np.zeros(points.shape, dtype=complex)
    sampled_field[:, 0] = points[:, 0] + 0.5j * points[:, 1]
    sampled_field[:, 1] = -points[:, 1]
    sampled_field[:, 2] = 0.25j * points[:, 2]
    sampled = integrate_time_average_maxwell_surface_force(
        sampled_field, normals, areas, amplitude="peak"
    )

    assert np.allclose(symbolic, sampled, rtol=1e-8, atol=1e-8 * np.abs(sampled).max())


def test_nonpositive_permeability_fails_loudly():
    mesh = _unit_box_with_named_right_face()
    field = CoefficientFunction((x, -y, 0.0))
    with pytest.raises(ValueError):
        maxwell_surface_force(
            field, mesh, mesh.Boundaries("right"), permeability_H_per_m=0.0
        )
    with pytest.raises(ValueError):
        time_average_maxwell_surface_force(
            field, mesh, mesh.Boundaries("right"), amplitude="neither"
        )


def test_module_constant_matches_the_sampled_kernel():
    assert MU0 == pytest.approx(4.0e-7 * math.pi, rel=0.0, abs=0.0)
