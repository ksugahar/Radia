"""Focused tests for the solver-ready Cubit .vol smoke-test gate."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

pytest.importorskip("ngsolve")

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_SRC = ROOT / "packages" / "cubit-mesh-export" / "src"
if str(PACKAGE_SRC) not in sys.path:
    sys.path.insert(0, str(PACKAGE_SRC))

import ngsolve as ng
from cubit_mesh_export import smoke_test
from netgen.meshing import (
    Element2D,
    Element3D,
    FaceDescriptor,
    MeshPoint,
    Pnt,
)
from netgen.meshing import (
    Mesh as NetgenMesh,
)


@pytest.fixture(scope="module", autouse=True)
def _taskmanager():
    with ng.TaskManager():
        yield


def _save_labeled_tetrahedron(tmp_path, *, order=1):
    # A fixed tetrahedron avoids Netgen's randomized CSG mesher making this
    # contract test intermittently degenerate. Complex real-Cubit geometry is
    # covered separately by test_cubit_vol_accuracy_dataset.py.
    ngmesh = NetgenMesh(dim=3)
    ngmesh.SetMaterial(1, "body")
    boundary = ngmesh.Add(
        FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1)
    )
    ngmesh.SetBCName(0, "outer")
    points = [
        ngmesh.Add(MeshPoint(Pnt(*coordinate)))
        for coordinate in (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        )
    ]
    ngmesh.Add(Element3D(1, points))
    for face in ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)):
        ngmesh.Add(Element2D(boundary, [points[index] for index in face]))
    ngmesh.Update()
    mesh = ng.Mesh(ngmesh)
    if order > 1:
        mesh.Curve(order)
    vol_path = tmp_path / "tetrahedron.vol"
    mesh.ngmesh.Save(str(vol_path))
    sidecar = Path(str(vol_path) + ".json")
    sidecar.write_text(
        json.dumps({
            "materials": {"body": 1.0 / 6.0},
            "mesh_only_materials": [],
            "boundaries": {"outer": 1.5 + (3.0 ** 0.5) / 2.0},
            "edges": {},
            "n_elements": mesh.ne,
            "n_points": mesh.nv,
            "order": mesh.GetCurveOrder(),
            "export_time_s": 0.01,
        }),
        encoding="utf-8",
    )
    return mesh, vol_path, sidecar


def test_solver_ready_validation_loads_vol_and_checks_complete_contract(tmp_path):
    mesh, vol_path, sidecar = _save_labeled_tetrahedron(tmp_path)

    result = smoke_test._validate_exported_vol(
        vol_path,
        order=1,
        expect=["outer"],
        expect_materials=["body"],
        threshold=1.0,
    )

    assert result["passed"], result["issues"]
    assert result["sidecar_path"] == sidecar
    assert result["report_path"].is_file()
    assert set(result["bcnames"]) == {"outer"}
    assert result["materials"] == ["body"]
    report = result["report"]
    assert report["mesh"]["n_elements"] == mesh.ne
    assert report["mesh"]["n_points"] == mesh.nv
    assert report["quality"]["passed"]
    assert report["quality"]["tetrahedron_count"] == mesh.ne
    assert report["boundary_domain_ownership"]["passed"]
    assert report["boundary_domain_ownership"][
        "unreferenced_volume_domains"
    ] == []
    assert report["materials"][0]["ng_volume"] == pytest.approx(1.0 / 6.0)
    assert report["boundaries"][0]["ng_area"] == pytest.approx(
        1.5 + (3.0 ** 0.5) / 2.0
    )


def test_solver_ready_validation_fails_stale_sidecar_metadata(tmp_path):
    _, vol_path, sidecar = _save_labeled_tetrahedron(tmp_path)
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    payload["order"] = 2
    sidecar.write_text(json.dumps(payload), encoding="utf-8")

    result = smoke_test._validate_exported_vol(
        vol_path,
        order=1,
        expect=["outer"],
        expect_materials=["body"],
        threshold=1.0,
    )

    assert not result["passed"]
    assert result["report_path"].is_file()
    assert any("CAD metadata order" in issue for issue in result["issues"])


def test_solver_ready_validation_samples_high_order_map_per_element(tmp_path):
    mesh, vol_path, _ = _save_labeled_tetrahedron(tmp_path, order=2)

    result = smoke_test._validate_exported_vol(
        vol_path,
        order=2,
        expect=["outer"],
        expect_materials=["body"],
        threshold=1.0,
    )

    assert result["passed"], result["issues"]
    quality = result["report"]["quality"]
    assert result["report"]["mesh"]["curve_order"] == 2
    assert quality["volume_element_count"] == mesh.ne
    assert quality["mapping_sample_count"] > mesh.ne
    assert quality["invalid_jacobian_sample_count"] == 0


@pytest.mark.parametrize(
    ("body", "message"),
    [
        ("mesh3d\nmaterials\n1\n1 body\nendmesh\n", "bcnames.*missing"),
        ("mesh3d\nbcnames\n2\n1 source\nendmesh\n", "truncated"),
        (
            "mesh3d\nbcnames\n2\n1 source\n1 sink\nendmesh\n",
            "repeats id",
        ),
    ],
)
def test_named_section_parser_rejects_silent_partial_or_duplicate_data(
    tmp_path, body, message,
):
    vol_path = tmp_path / "bad.vol"
    vol_path.write_text(body, encoding="ascii")

    with pytest.raises(ValueError, match=message):
        smoke_test._read_vol_bcnames(vol_path)


def test_named_section_parser_allows_one_name_on_multiple_face_descriptors(
    tmp_path,
):
    vol_path = tmp_path / "shared_name.vol"
    vol_path.write_text(
        "mesh3d\nbcnames\n2\n1 outer\n2 outer\nendmesh\n",
        encoding="ascii",
    )

    assert smoke_test._read_vol_bcnames(vol_path) == ["outer", "outer"]


def test_solver_ready_validation_requires_companion_cad_sidecar(tmp_path):
    _, vol_path, sidecar = _save_labeled_tetrahedron(tmp_path)
    sidecar.unlink()

    with pytest.raises(RuntimeError, match="required CAD sidecar"):
        smoke_test._validate_exported_vol(
            vol_path,
            order=1,
            expect=["outer"],
            expect_materials=["body"],
            threshold=1.0,
        )


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda payload: payload.pop("n_points"), "missing required keys.*n_points"),
        (
            lambda payload: payload.__setitem__("export_time_s", "NaN"),
            "export_time_s is invalid",
        ),
    ],
)
def test_solver_ready_validation_rejects_invalid_sidecar_contract(
    tmp_path, mutate, message,
):
    _, vol_path, sidecar = _save_labeled_tetrahedron(tmp_path)
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    mutate(payload)
    sidecar.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(RuntimeError, match=message):
        smoke_test._validate_exported_vol(
            vol_path,
            order=1,
            expect=["outer"],
            expect_materials=["body"],
            threshold=1.0,
        )
