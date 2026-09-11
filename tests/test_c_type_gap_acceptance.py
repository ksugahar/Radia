"""Fast probe-gate checks without Cubit or a solver mesh."""
import importlib.util
from pathlib import Path

import pytest


path = Path(__file__).resolve().parents[1] / "validation_test/c_type_three_engine/build_cubit_meshes.py"
spec = importlib.util.spec_from_file_location("c_type_gap_builder", path)
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def inventory():
    return {"elements": 12, "gap_height_m": 0.01, "line_profile": {
        "lines": [{"curved": {"covered_m": 0.01}}],
        "minimum_segments": 6, "curved_minimum_segments": 6,
        "clean_traversal": True, "sampling_is_stable": True,
        "maximum_segment_m": 0.0018, "curved_maximum_segment_m": 0.0018}}


def test_complete_probe_passes():
    passed, limit = builder._gap_acceptance(inventory(), 6, 1.2)
    assert passed
    assert limit == pytest.approx(0.002)


@pytest.mark.parametrize("covered", [0.009, 0.011, float("nan")])
def test_curved_coverage_failure_rejected(covered):
    data = inventory()
    data["line_profile"]["lines"][0]["curved"]["covered_m"] = covered
    assert not builder._gap_acceptance(data, 6, 1.2)[0]


@pytest.mark.parametrize("field,value", [
    ("lines", []), ("clean_traversal", False), ("sampling_is_stable", False),
    ("minimum_segments", 5), ("curved_minimum_segments", 5),
    ("maximum_segment_m", 0.003), ("curved_maximum_segment_m", 0.003)])
def test_incomplete_probe_rejected(field, value):
    data = inventory()
    data["line_profile"][field] = value
    assert not builder._gap_acceptance(data, 6, 1.2)[0]


@pytest.mark.parametrize("required,factor", [(0, 1.2), (-1, 1.2),
    (1.5, 1.2), (6, 0), (6, float("nan")), (6, float("inf"))])
def test_invalid_requirements_rejected(required, factor):
    with pytest.raises(ValueError):
        builder._gap_acceptance(inventory(), required, factor)


GAP = 0.01
BAND = 2.0e-7


class _BiasedLocatorMesh:
    """A stack of tetrahedra in z whose point locator mimics NGSolve's tolerance.

    Elements are scanned from the top down and accepted up to ``BAND`` outside
    their faces, so just below every face the locator returns the element
    ABOVE it with a negative reference margin -- the behaviour measured on the
    C-type family meshes (iron returned 0.2 um inside the air gap).
    """

    faces = (-1.0, -GAP / 2, -GAP / 6, GAP / 6, GAP / 2, 1.0)
    element_material = (0, 1, 1, 1, 0)

    def GetMaterials(self):
        return ("iron", "air")

    def __call__(self, x, y, z):
        from types import SimpleNamespace
        for number in reversed(range(len(self.faces) - 1)):
            low, high = self.faces[number], self.faces[number + 1]
            if low - BAND <= z <= high + BAND:
                margin = min(z - low, high - z) / (high - low)
                return SimpleNamespace(nr=number,
                                       pnt=(margin, (1 - margin) / 3, (1 - margin) / 3))
        return SimpleNamespace(nr=-1, pnt=(0.0, 0.0, 0.0))

    def __getitem__(self, element_id):
        import ngsolve as ng
        from types import SimpleNamespace
        return SimpleNamespace(index=self.element_material[element_id.nr],
                               type=ng.ET.TET)


def test_fake_locator_reproduces_the_band():
    mesh = _BiasedLocatorMesh()
    number, margin, material = builder._containing_element(
        mesh, mesh.GetMaterials(), 0.0, 0.0, GAP / 2 - BAND / 2)
    assert material == "iron" and margin < 0.0


def test_curved_walk_puts_switches_on_the_faces_not_the_band():
    result = builder._curved_line_segments(_BiasedLocatorMesh(), 0.0, 0.0, GAP / 2)
    assert result["segments"] == 3
    assert abs(result["covered_m"] - GAP) < 1.0e-12
    assert result["maximum_segment_m"] == pytest.approx(GAP / 3, abs=2.0e-9)
    assert result["minimum_segment_m"] == pytest.approx(GAP / 3, abs=2.0e-9)
    assert result["locator_band_samples"] >= 1


def _sliver_mesh():
    """Insert a 3 um air element strictly between two coarse walk samples."""
    import numpy as np

    grid = np.linspace(-GAP / 2 + 1e-9, GAP / 2 - 1e-9, 2001)
    start = float(grid[700]) + 1.0e-6
    mesh = _BiasedLocatorMesh()
    mesh.faces = (-1.0, -GAP / 2, start, start + 3.0e-6, GAP / 6, GAP / 2, 1.0)
    mesh.element_material = (0, 1, 1, 1, 1, 0)
    return mesh


def test_unseeded_walk_steps_over_a_thin_element():
    result = builder._curved_line_segments(_sliver_mesh(), 0.0, 0.0, GAP / 2)
    assert result["segments"] == 3


def test_seeds_make_the_walk_visit_every_element():
    mesh = _sliver_mesh()
    faces = mesh.faces[1:-1]
    seeds = [0.5 * (low + high) for low, high in zip(faces, faces[1:])]
    result = builder._curved_line_segments(mesh, 0.0, 0.0, GAP / 2, seeds=seeds)
    assert result["segments"] == 4
    assert result["minimum_segment_m"] == pytest.approx(3.0e-6, abs=2.0e-9)
    assert abs(result["covered_m"] - GAP) < 1.0e-12


def test_point_outside_the_mesh_has_no_element():
    mesh = _BiasedLocatorMesh()
    assert builder._containing_element(mesh, mesh.GetMaterials(), 0.0, 0.0, 2.0) \
        == (None, None, None)
