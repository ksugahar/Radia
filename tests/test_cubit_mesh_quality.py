"""Fast production gates for curved Cubit/Netgen volume meshes."""

from pathlib import Path
import importlib.util

import pytest


pytest.importorskip("ngsolve")

PACKAGE_SRC = Path(__file__).resolve().parents[1] / "packages" / "cubit-mesh-export" / "src"
CHECK_FILE = PACKAGE_SRC / "cubit_mesh_export" / "check.py"
SPEC = importlib.util.spec_from_file_location("cubit_mesh_export_check_current", CHECK_FILE)
CHECK_MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(CHECK_MODULE)
check_mesh_quality = CHECK_MODULE.check_mesh_quality
from netgen.csg import CSGeometry, OrthoBrick, Pnt, Sphere  # noqa: E402
import ngsolve as ng  # noqa: E402


def _conductor_box():
    geometry = CSGeometry()
    geometry.Add(
        OrthoBrick(Pnt(0.0, 0.0, 0.0), Pnt(1.0, 1.0, 1.0))
        .mat("cond")
        .bc("sibc")
    )
    return ng.Mesh(geometry.GenerateMesh(maxh=0.8))


def test_curved_mapping_gate_samples_actual_high_order_transformation():
    geometry = CSGeometry()
    geometry.Add(Sphere(Pnt(0.0, 0.0, 0.0), 1.0).mat("body").bc("outer"))
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=0.7))
    mesh.Curve(3)

    result = check_mesh_quality(
        mesh,
        min_curve_order=3,
        require_tetrahedra=True,
        required_materials=("body",),
        required_boundaries=("outer",),
    )

    assert result["passed"]
    assert result["curve_order"] == 3
    assert result["integration_order"] == 8
    assert result["mapping_sample_count"] > result["volume_element_count"]
    assert result["minimum_jacobian"] > 0.0
    assert result["minimum_scaled_jacobian"] > 0.0


def test_conductor_face_gate_separates_sibc_surface_from_loop_bridges():
    result = check_mesh_quality(
        _conductor_box(),
        conductive_materials=("cond",),
        sibc_boundaries=("sibc",),
        require_all_sibc_labeled=True,
    )

    adjacency = result["adjacency"]
    assert result["passed"]
    assert adjacency["sibc_candidate_face_count"] > 0
    assert adjacency["loop_bridge_face_count"] > 0
    assert adjacency["face_role_counts"]["conductor-exterior"] > 0
    assert adjacency["face_role_counts"]["conductor-conductor"] > 0
    assert not adjacency["marked_non_sibc_faces"]
    assert not adjacency["unlabeled_sibc_faces"]


def test_mesh_gate_reports_quality_label_and_sibc_failures():
    result = check_mesh_quality(
        _conductor_box(),
        min_scaled_jacobian=0.99,
        min_curve_order=2,
        required_materials=("cond", "missing_material"),
        required_boundaries=("sibc", "missing_boundary"),
        conductive_materials=("cond",),
        sibc_boundaries=("not_the_surface",),
        require_all_sibc_labeled=True,
    )

    assert not result["passed"]
    assert result["low_scaled_jacobian_sample_count"] > 0
    assert result["missing_materials"] == ["missing_material"]
    assert result["missing_boundaries"] == ["missing_boundary"]
    assert result["adjacency"]["unlabeled_sibc_faces"]
    assert any("curve order" in warning for warning in result["warnings"])
