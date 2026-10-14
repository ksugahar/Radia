from pathlib import Path

from validation_test.c_type_three_engine.cubit_reflection_mesh import (
    _material_boundary_surfaces,
)


ROOT = Path(__file__).resolve().parents[1]


class _CubitTopology:
    def __init__(self):
        self._volume_surfaces = {
            1: [10, 11, 12],
            2: [12, 13, 14],
        }
        self._surface_volumes = {
            10: [1],
            11: [1, 3],
            12: [1, 2],
            13: [2, 4],
            14: [2],
        }

    def get_relatives(self, source, entity, target):
        if (source, target) == ("volume", "surface"):
            return self._volume_surfaces[entity]
        if (source, target) == ("surface", "volume"):
            return self._surface_volumes[entity]
        raise AssertionError((source, entity, target))


def test_material_boundary_excludes_only_same_material_internal_seam():
    assert _material_boundary_surfaces(_CubitTopology(), [1, 2]) == {
        10, 11, 13, 14,
    }


def test_exporter_fails_loudly_instead_of_deleting_explicit_internal_label():
    source = (
        ROOT / "src" / "cubit_plugin" / "ExportNetgenCommand.cpp"
    ).read_text(encoding="utf-8")
    curver = (
        ROOT / "src" / "cubit_plugin" / "NetgenCurver.cpp"
    ).read_text(encoding="utf-8")

    assert "labelled_same_material_surfaces" in source
    assert 'label == "Surface " + std::to_string(surface_id)' in source
    assert "Same-material internal surface %d is explicitly labelled" in source
    assert "or model the intended interface/cut with distinct domains" in source
    assert source.index("labelled_same_material_surfaces") < source.index(
        "DeleteElement(sei)")
    assert "same_domain_internal_surfaces" not in curver
