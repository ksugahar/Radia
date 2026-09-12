"""Vertex reflection audit regressions without Cubit or NGSolve."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest


path = Path(__file__).resolve().parents[1] / "validation_test/c_type_three_engine/build_cubit_meshes.py"
spec = importlib.util.spec_from_file_location("c_type_reflection_builder", path)
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def reference(coordinates):
    keys = {tuple(np.round(point, 14)) for point in coordinates}
    missing, error = 0, 0.0
    for point in coordinates:
        reflected = point * (1, 1, -1)
        if tuple(np.round(reflected, 14)) not in keys:
            missing += 1
            continue
        candidates = coordinates[np.all(np.isclose(
            coordinates, reflected, rtol=0, atol=1e-13), axis=1)]
        error = max(error, float(np.min(np.linalg.norm(candidates - reflected, axis=1))))
    return {"vertex_count": len(coordinates), "missing_reflected_vertices": missing,
            "maximum_reflected_vertex_error_m": error}


@pytest.mark.parametrize("coordinates", [
    [[0, 0, -1], [4e-15, 0, 1], [0, 0, 1]],
    [[4.9e-15, 0, -1], [-4.9e-15, 0, 1], [5.1e-15, 0, 1]],
    [[0, 0, -1], [6e-15, 0, 1]],
    [[0, 0, 0], [0, 0, 0]],
])
def test_collisions_rounding_boundaries_and_duplicates_match_reference(coordinates):
    points = np.asarray(coordinates, dtype=float)
    assert builder._reflected_vertex_inventory(points) == reference(points)
    assert builder._reflected_vertex_inventory(points[::-1]) == reference(points)


def test_perturbed_symmetric_cloud_matches_reference():
    rng = np.random.default_rng(927)
    upper = rng.uniform(0.01, 0.1, (80, 3))
    lower = upper * (1, 1, -1) + rng.uniform(-4e-15, 4e-15, upper.shape)
    points = np.vstack([upper, lower])
    assert builder._reflected_vertex_inventory(points) == reference(points)


def test_exact_mirror_does_not_construct_spatial_tree(monkeypatch):
    import scipy.spatial

    def forbidden(*args, **kwargs):
        raise AssertionError("exact mirrors need no tree")

    monkeypatch.setattr(scipy.spatial, "cKDTree", forbidden)
    result = builder._reflected_vertex_inventory([[0, 0, -1], [0, 0, 1]])
    assert result["maximum_reflected_vertex_error_m"] == 0


@pytest.mark.parametrize("points", [[[0, 0, float("nan")]], [[0, 0, float("inf")]], [[0, 0]]])
def test_invalid_coordinates_are_rejected(points):
    with pytest.raises(ValueError):
        builder._reflected_vertex_inventory(points)


def test_empty_mesh_vertices():
    assert builder._reflected_vertex_inventory(np.empty((0, 3))) == reference(np.empty((0, 3)))
