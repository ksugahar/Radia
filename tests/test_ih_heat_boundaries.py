import importlib.util
import json
import sys
from pathlib import Path

import pytest

SOURCE = Path(__file__).parents[1] / "src" / "radia" / "ih_heat_boundaries.py"
SPEC = importlib.util.spec_from_file_location("ih_heat_boundaries_under_test", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)
load_convection_map = MODULE.load_convection_map
normalize_convection_map = MODULE.normalize_convection_map


def test_convection_map_normalizes_exact_labels_and_insulation():
    records = normalize_convection_map({
        "top": {"h_W_m2K": 25, "ambient_C": 20},
        "shield[1]": {"h_W_m2K": 0, "ambient_C": -5},
    })
    assert [(r.label, r.h_W_m2K, r.ambient_C) for r in records] == [
        ("top", 25.0, 20.0), ("shield[1]", 0.0, -5.0)]
    assert normalize_convection_map({}) == ()


@pytest.mark.parametrize("mapping, message", [
    ({"top": {"h_W_m2K": -1, "ambient_C": 20}}, "nonnegative"),
    ({"top": {"h_W_m2K": float("nan"), "ambient_C": 20}}, "finite"),
    ({"top": {"h_W_m2K": 1, "ambient_C": float("inf")}}, "finite"),
    ({"top": {"h_W_m2K": 1, "ambient_C": -273.16}}, "absolute zero"),
    ({"top": {"h_W_m2K": 1}}, "requires exactly"),
    ({"top": {"h_W_m2K": 1, "ambient_C": 20, "extra": 0}},
     "requires exactly"),
])
def test_convection_map_rejects_invalid_coefficients(mapping, message):
    with pytest.raises(ValueError, match=message):
        normalize_convection_map(mapping)


def test_json_loader_rejects_duplicate_labels(tmp_path):
    path = tmp_path / "convection.json"
    path.write_text(
        '{"top":{"h_W_m2K":10,"ambient_C":20},'
        '"top":{"h_W_m2K":30,"ambient_C":40}}', encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate convection boundary label"):
        load_convection_map(path)


@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig"])
def test_json_loader_returns_canonical_mapping(tmp_path, encoding):
    path = tmp_path / "convection.json"
    path.write_text(json.dumps({
        "top": {"h_W_m2K": 10, "ambient_C": 20}}), encoding=encoding)
    assert load_convection_map(path) == {
        "top": {"h_W_m2K": 10.0, "ambient_C": 20.0}}
