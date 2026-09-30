from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    ROOT / "validation_test" / "team21" / "geometry_source.py"
)


def _load_team21a():
    module_name = "radia_team21a_validation"
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(module_name, MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def test_team21a_public_spec_and_slit_layout_are_locked():
    team21a = _load_team21a()

    assert team21a.SPEC.frequency_hz == 50.0
    assert team21a.SPEC.coil_turns == 300
    assert team21a.SPEC.coil_current_rms_A == 10.0
    assert team21a.SPEC.conductivity_S_per_m == 1.3889e6
    assert team21a.slit_centers_m(0) == ()
    assert team21a.slit_centers_m(1) == (0.0,)
    assert team21a.slit_centers_m(2) == (-0.060, 0.060)
    assert team21a.slit_centers_m(3) == (-0.090, 0.0, 0.090)
    assert team21a.SPEC.skin_depth_m == pytest.approx(0.06039481294661581)
    assert team21a.SPEC.plate_thickness_m / team21a.SPEC.skin_depth_m < 0.17


def test_team21a_radia_mother_source_has_expected_opposed_coil_field():
    team21a = _load_team21a()
    import radia

    source, ports = team21a.build_mother_source(
        radial_order=2,
        axial_order=3,
        arc_segments=3,
    )
    try:
        assert len(ports) == 2
        center = np.asarray(radia.Fld(source, "b", [0.0, 0.0, 0.0]))
        assert center[0] == pytest.approx(-1.0e-2, rel=0.08)
        assert abs(center[1]) < 1.0e-4
        assert abs(center[2]) < 1.0e-7
    finally:
        radia.UtiDelAll()


@pytest.mark.parametrize("model", [0, 1, 2, 3])
def test_team21a_slitted_conductor_mesh_preserves_volume_and_topology(model):
    ng = pytest.importorskip("ngsolve")
    team21a = _load_team21a()
    import radia.vim as vim

    mesh = team21a.build_conductor_mesh(model, maxh=0.12)
    assert mesh.GetMaterials() == ("cond",)
    assert "skin" in mesh.GetBoundaries()
    expected_volume = (
        team21a.SPEC.plate_width_m
        * team21a.SPEC.plate_length_m
        * team21a.SPEC.plate_thickness_m
        - model
        * team21a.SPEC.slit_length_m
        * team21a.SPEC.slit_width_m
        * team21a.SPEC.plate_thickness_m
    )
    actual_volume = ng.Integrate(1.0, mesh, definedon=mesh.Materials("cond"))
    assert actual_volume == pytest.approx(expected_volume, rel=1.0e-10)

    topology = vim.ClassifyNgsolveEddyTopology(mesh, conductive_materials="cond")
    diagnostics = topology.diagnostics()
    assert diagnostics["conductive_component_count"] == 1
    assert diagnostics["sibc_face_count"] > 0
    assert diagnostics["conductor_conductor_face_count"] > 0

