"""Curved RT2 full-domain/IMA production contract."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from netgen.meshing import (  # noqa: E402
    Element2D,
    Element3D,
    FaceDescriptor,
    Mesh as NetgenMesh,
    MeshPoint,
    Pnt,
)

from radia import vim  # noqa: E402


def _mirror_pair(full):
    mesh = NetgenMesh(dim=3)
    coordinates = [(0, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0)]
    if full:
        coordinates.append((-1, 0, 0))
    points = [mesh.Add(MeshPoint(Pnt(*point))) for point in coordinates]
    mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    mesh.Add(Element3D(1, [points[i] for i in (0, 1, 2, 3)]))
    for face in ((0, 1, 3), (1, 2, 3), (2, 0, 3)):
        mesh.Add(Element2D(1, [points[i] for i in face]))
    if full:
        mesh.Add(Element3D(1, [points[i] for i in (0, 2, 1, 4)]))
        for face in ((0, 2, 4), (2, 1, 4), (1, 0, 4)):
            mesh.Add(Element2D(1, [points[i] for i in face]))
    else:
        mesh.Add(Element2D(1, [points[i] for i in (0, 2, 1)]))
    mesh.SetMaterial(1, "iron")
    return ng.Mesh(mesh)


def test_curved_rt2_material_solve_and_field_match_full_model_to_roundoff():
    observations = np.array([
        [1.5, 0.2, 0.2], [-1.5, 0.2, 0.2],
        [0.2, 1.5, 0.2], [-0.2, -1.5, 0.2],
        [0.2, 0.2, 1.5], [-0.2, 0.2, -1.5],
        [2.0, 2.0, 2.0], [-2.0, -2.0, -2.0],
    ])
    applied_scale = 5000.0
    kwargs = dict(
        mu_r=100.0,
        H_ext=ng.CoefficientFunction((0.0, 0.0, applied_scale)),
        order=2,
        curve_order=2,
        gram_eps=1e-12,
        tol=1e-12,
        ho_far_factor=float("inf"),
    )
    with ng.TaskManager():
        full = vim.Solve(_mirror_pair(True), **kwargs)
        half = vim.Solve(_mirror_pair(False), image="+x", **kwargs)
        full_field = vim.FieldFromSolution(full, observations)
        half_field = vim.FieldFromSolution(half, observations)

    eps10 = 10.0*np.finfo(float).eps
    field_scale = np.maximum(np.linalg.norm(full_field, axis=1), applied_scale)
    field_delta = full_field-half_field
    field_component_error = np.max(np.abs(field_delta)/field_scale[:, None])
    field_vector_error = np.max(np.linalg.norm(field_delta, axis=1)/field_scale)
    magnetization_error = (
        np.linalg.norm(full["M_avg"]-half["M_avg"])
        / np.linalg.norm(full["M_avg"])
    )
    assert abs(full["demag"]-half["demag"]) < eps10
    # The full and reduced CG systems have different dimensions and therefore
    # different reduction/update orders.  Their solved means stay in a small
    # roundoff band.  Each field component remains below 10 eps; the Euclidean
    # norm consequently has the sharp three-component bound sqrt(3)*10 eps.
    assert magnetization_error < 32.0*np.finfo(float).eps
    assert field_component_error < eps10
    assert field_vector_error < np.sqrt(3.0)*eps10
    assert full["field_evaluator_stats"]["source_kind"] == "curved-element-exact-bdm2"
    assert half["field_evaluator_stats"]["source_kind"] == "curved-element-exact-bdm2"
