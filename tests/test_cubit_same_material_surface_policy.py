from validation_test.c_type_three_engine.cubit_reflection_mesh import (
    _material_boundary_surfaces,
)


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
