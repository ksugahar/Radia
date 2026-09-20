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
