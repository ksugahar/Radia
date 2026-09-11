"""add_kelvin_cubit on a real Cubit model with gap slabs inside the air.

Requires a Cubit seat (LAB or 100 only); run headless:

    python -m pytest validation_test/c_type_three_engine/test_add_kelvin_sphere_faces.py

The air domain is built by the same helper as the C-type validation meshes
(cubit_reflection_mesh.build_reflection_invariant_physical_mesh with gap
slabs), around a small iron brick instead of the full C-yoke, so the test
checks the production path in seconds.  Before the fix the slab faces went
into kelvin_int; at 24 slabs the Netgen export then failed.
"""

import importlib.util
import math
import os
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# The same Cubit the mesh builder drives.
_cubit_path = str(_load("c_type_builder_under_test",
                        HERE / "build_cubit_meshes.py").CUBIT_DEFAULT.parent)
if not (Path(_cubit_path) / "cubit.py").is_file():
    pytest.skip("Coreform Cubit is not installed on this host", allow_module_level=True)
sys.path.insert(0, _cubit_path)
# Under pytest, validation_test/ is on sys.path and validation_test/cubit is a
# regular package that shadows Cubit's own module, so load it by path and
# register it; the helpers under test do `import cubit` inside functions.
sys.modules.pop("cubit", None)
_spec = importlib.util.spec_from_file_location("cubit", Path(_cubit_path) / "cubit.py")
cubit = importlib.util.module_from_spec(_spec)
sys.modules["cubit"] = cubit
_spec.loader.exec_module(cubit)

cubit.init(["cubit", "-nojournal", "-batch", "-nographics"])

reflection = _load("c_type_reflection_under_test",
                   REPO / "validation_test" / "c_type_three_engine" / "cubit_reflection_mesh.py")
add_kelvin = _load("add_kelvin_under_test",
                   REPO / "packages" / "cubit-mesh-export" / "src" / "cubit_mesh_export"
                   / "cubit_helpers" / "add_kelvin.py")

R = 0.22


# Cubit's start-up reads the deployed Radia toolbar script and leaves the file
# open; the ResourceWarning surfaces at the first teardown and is not about
# the code under test.
pytestmark = [pytest.mark.filterwarnings("ignore::ResourceWarning"),
              pytest.mark.filterwarnings("ignore::pytest.PytestUnraisableExceptionWarning")]


# Probed on Cubit 2025.12: get_block_name and get_exodus_entity_name return
# the names; there is no get_sideset_name, and get_entity_name rejects blocks.
def _block_id(name):
    return next(b for b in cubit.get_block_id_list()
                if (cubit.get_block_name(b) or "").lower() == name)


def _sideset_id(name):
    return next(s for s in cubit.get_sideset_id_list()
                if (cubit.get_exodus_entity_name("sideset", s) or "").lower() == name)


def _layered_air(gap_layers):
    cubit.cmd("reset")
    # a 40 x 30 mm iron pole with a yoke block on top; the pole's bottom face
    # (1.2e-3 m^2) is the pole face the helper refines
    cubit.cmd("create brick x 0.040 y 0.030 z 0.015")
    cubit.cmd(f"move volume {cubit.get_last_id('volume')} z 0.0125 include_merged")
    cubit.cmd("create brick x 0.040 y 0.030 z 0.015")
    cubit.cmd(f"move volume {cubit.get_last_id('volume')} z 0.0275 include_merged")
    return reflection.build_reflection_invariant_physical_mesh(
        iron_size=0.02, air_size=0.08, gap_size=0.005, kelvin_radius=R,
        gap_layers=gap_layers)


@pytest.mark.parametrize("gap_layers", [0, 4])
def test_kelvin_int_holds_only_the_air_sphere(gap_layers):
    physical = _layered_air(gap_layers)
    assert len(physical["gap_layers_up"]) == gap_layers // 2
    air_volumes = cubit.parse_cubit_list("volume", f"in block {_block_id('air')}")
    assert len(air_volumes) == 2 + gap_layers

    add_kelvin.add_kelvin_cubit(R=R, air_block="air", symmetry=["z"], mesh_size=0.15)

    faces = cubit.get_sideset_surfaces(_sideset_id("kelvin_int"))
    assert len(faces) == 2
    for face in faces:
        assert cubit.get_surface_type(face) == "sphere surface"
        near = cubit.surface(face).closest_point_trimmed([0.0, 0.0, 0.0])
        assert abs(math.dist(near, (0.0, 0.0, 0.0)) - R) <= 1.0e-9
    kelvin = cubit.parse_cubit_list("volume", f"in block {_block_id('kelvin')}")
    assert len(kelvin) == 2


def test_air_without_a_sphere_face_is_rejected_before_any_geometry():
    cubit.cmd("reset")
    cubit.cmd("create brick x 0.1 y 0.1 z 0.1")
    air = cubit.get_last_id("volume")
    cubit.cmd(f"volume {air} scheme tetmesh")
    cubit.cmd(f"mesh volume {air}")
    cubit.cmd(f"block 1 add volume {air}")
    cubit.cmd('block 1 name "air"')
    before = set(cubit.parse_cubit_list("volume", "all"))
    with pytest.raises(RuntimeError, match="No face"):
        add_kelvin.add_kelvin_cubit(R=R, air_block="air", symmetry=["z"])
    assert set(cubit.parse_cubit_list("volume", "all")) == before


def test_air_sphere_off_the_origin_is_rejected_before_any_geometry():
    cubit.cmd("reset")
    cubit.cmd(f"create sphere radius {R}")
    air = cubit.get_last_id("volume")
    cubit.cmd(f"move volume {air} x 0.01 include_merged")
    cubit.cmd(f"volume {air} scheme tetmesh")
    cubit.cmd(f"volume {air} size 0.08")
    cubit.cmd(f"mesh volume {air}")
    cubit.cmd(f"block 1 add volume {air}")
    cubit.cmd('block 1 name "air"')
    before = set(cubit.parse_cubit_list("volume", "all"))
    with pytest.raises(NotImplementedError, match="origin"):
        add_kelvin.add_kelvin_cubit(R=R, air_block="air")
    assert set(cubit.parse_cubit_list("volume", "all")) == before
