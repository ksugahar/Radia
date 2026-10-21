"""Persistent ESIM cache invalidation; synthetic material only."""
from dataclasses import replace
import json
import numpy as np
import pytest
from radia.vim._eddy_hybrid import LocalESIMSurfaceModel, LocalESIMSurfaceLUT


def model():
    return LocalESIMSurfaceModel(
        bh_curve=np.array([[0., 0.], [100., .1], [1000., .5]]), sigma=2.e6)


def table(m):
    return LocalESIMSurfaceLUT(
        frequency_nodes_hz=np.array([1.e4, 2.e4]),
        h_nodes_A_per_m=np.array([1., 10.]),
        impedance_ohm=np.full((2, 2), .01+.01j),
        cell_iterations=np.ones((2, 2), dtype=int),
        model_signature=m.lut_signature())


def test_old_persistent_schema_is_rejected(tmp_path):
    current = table(model())
    path = current.save(tmp_path / "current.npz")
    with np.load(path, allow_pickle=False) as archive:
        contents = {name: archive[name] for name in archive.files}
    metadata = json.loads(str(contents["metadata"].item()))
    metadata["schema_version"] = 1
    metadata["model_signature"]["schema_version"] = 1
    for name in ("cell_algorithm", "mesh_tol", "max_nodes", "max_refinements"):
        metadata["model_signature"].pop(name)
    contents["metadata"] = np.array(json.dumps(metadata))
    old = tmp_path / "old.npz"
    np.savez(old, **contents)
    with pytest.raises(ValueError, match="unsupported local ESIM LUT schema"):
        LocalESIMSurfaceLUT.load(old)


@pytest.mark.parametrize("change", [
    {"mesh_tol": 2.e-3}, {"max_nodes": 16385}, {"max_refinements": 8}])
def test_depth_controls_invalidate_lut_and_reach_solver(change):
    original = model()
    changed = replace(original, **change)
    lut = table(original)
    assert changed.lut_signature() != original.lut_signature()
    with pytest.raises(ValueError, match="model signature does not match"):
        changed.with_lut(lut)
    solver = changed._cell_solver(1.e4)
    for name in ("mesh_tol", "max_nodes", "max_refinements"):
        assert getattr(solver, name) == getattr(changed, name)


def test_algorithm_identity_and_new_schema_roundtrip(tmp_path):
    m = model()
    signature = m.lut_signature()
    assert signature["schema_version"] == 2
    assert signature["cell_algorithm"] == "conservative-cv-adaptive-depth-exp-positive-iwt-v1"
    lut = table(m)
    loaded = LocalESIMSurfaceLUT.load(lut.save(tmp_path / "new.npz"))
    loaded.validate_model(m)
    assert loaded.model_key == lut.model_key
    old_algorithm = dict(signature, cell_algorithm="geometric-1.1")
    with pytest.raises(ValueError, match="model signature does not match"):
        replace(lut, model_signature=old_algorithm).validate_model(m)


@pytest.mark.parametrize("change", [
    {"mesh_tol": 0.}, {"mesh_tol": float("nan")},
    {"max_nodes": 7}, {"max_refinements": 0}, {"max_refinements": 1.5}])
def test_invalid_depth_controls_fail_before_cell_evaluation(change):
    with pytest.raises(ValueError):
        replace(model(), **change)
