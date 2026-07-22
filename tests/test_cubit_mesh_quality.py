"""Fast production gates for curved Cubit/Netgen volume meshes."""

import json
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
VERIFY_FILE = (
    Path(__file__).resolve().parents[1]
    / "src" / "radia" / "panels" / "calc_verify_vol.py"
)
VERIFY_SPEC = importlib.util.spec_from_file_location("radia_calc_verify_vol_current", VERIFY_FILE)
VERIFY_MODULE = importlib.util.module_from_spec(VERIFY_SPEC)
assert VERIFY_SPEC.loader is not None
VERIFY_SPEC.loader.exec_module(VERIFY_MODULE)
LABEL_CONTRACT_SCHEMA = CHECK_MODULE.LABEL_CONTRACT_SCHEMA
REPORT_SCHEMA = CHECK_MODULE.REPORT_SCHEMA
check_consistency = CHECK_MODULE.check_consistency
check_label_contract = CHECK_MODULE.check_label_contract
check_mesh_quality = CHECK_MODULE.check_mesh_quality
check_main = CHECK_MODULE.main
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


def _save_box(tmp_path, *, material="body", boundary="outer"):
    geometry = CSGeometry()
    geometry.Add(
        OrthoBrick(Pnt(0.0, 0.0, 0.0), Pnt(1.0, 1.0, 1.0))
        .mat(material)
        .bc(boundary)
    )
    mesh = ng.Mesh(geometry.GenerateMesh(maxh=0.8))
    vol_path = tmp_path / "box.vol"
    mesh.ngmesh.Save(str(vol_path))
    return mesh, vol_path


def _reverse_volume_orientations(source, target):
    lines = Path(source).read_text(encoding="utf-8").splitlines()
    output = []
    in_elements = False
    remaining = None
    for line in lines:
        if line.strip() == "volumeelements":
            in_elements = True
            remaining = None
            output.append(line)
            continue
        if in_elements and remaining is None:
            if not line.strip():
                output.append(line)
                continue
            remaining = int(line.strip())
            output.append(line)
            continue
        if in_elements and remaining and line.strip():
            fields = line.split()
            fields[-1], fields[-2] = fields[-2], fields[-1]
            output.append(" ".join(fields))
            remaining -= 1
            if remaining == 0:
                in_elements = False
            continue
        output.append(line)
    Path(target).write_text("\n".join(output) + "\n", encoding="utf-8")


class _LabelMesh:
    ne = 1

    def __init__(self, materials, boundaries, bboundaries=(), bbboundaries=()):
        self._labels = (materials, boundaries, bboundaries, bbboundaries)

    def GetMaterials(self):
        return self._labels[0]

    def GetBoundaries(self):
        return self._labels[1]

    def GetBBoundaries(self):
        return self._labels[2]

    def GetBBBoundaries(self):
        return self._labels[3]


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


def test_consistently_negative_element_orientation_is_valid(tmp_path):
    mesh, vol_path = _save_box(tmp_path)
    reversed_path = tmp_path / "box_negative_orientation.vol"
    _reverse_volume_orientations(vol_path, reversed_path)

    result = check_mesh_quality(reversed_path)

    assert result["passed"]
    assert result["negative_orientation_element_count"] == mesh.ne
    assert result["positive_orientation_element_count"] == 0
    assert result["negative_jacobian_sample_count"] > 0
    assert result["invalid_jacobian_sample_count"] == 0
    assert result["mapping_sample_count"] == mesh.ne
    assert result["minimum_absolute_jacobian"] > 0.0


def test_conductor_face_gate_separates_sibc_surface_from_loop_bridges():
    result = check_mesh_quality(
        _conductor_box(),
        conductors=("cond",),
        sibc_boundaries=("sibc",),
        require_all_sibc_labeled=True,
        tet_only=True,
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


def test_standalone_check_passes_without_optional_cad_sidecar(tmp_path):
    _, vol_path = _save_box(tmp_path)

    result = check_consistency(vol_path)

    assert result["schema"] == REPORT_SCHEMA
    assert result["passed"]
    assert result["cad_reference"]["available"] is False
    assert result["mesh"]["n_elements"] > 0
    assert result["labels"]["materials"] == ["body"]
    assert result["labels"]["boundaries"] == ["outer"]
    assert result["quality"]["passed"]


def test_auto_discovered_cad_sidecar_checks_geometry_and_metadata(tmp_path):
    mesh, vol_path = _save_box(tmp_path)
    sidecar = Path(str(vol_path) + ".json")
    sidecar.write_text(
        json.dumps({
            "materials": {"body": 1.0},
            "boundaries": {"outer": 6.0},
            "edges": {},
            "n_elements": mesh.ne,
            "n_points": mesh.nv,
            "order": mesh.GetCurveOrder(),
        }),
        encoding="utf-8",
    )

    result = check_consistency(vol_path)

    assert result["passed"]
    assert result["cad_reference"]["available"]
    assert result["cad_reference"]["auto_discovered"]
    assert all(
        check["passed"]
        for check in result["cad_reference"]["metadata_checks"].values()
    )


def test_strict_label_gate_reports_generated_names_collisions_and_pairs():
    mesh = _LabelMesh(
        ("Copper", "copper", "volume_2"),
        ("source", "Surface_7", "sym_bn=0_x", "sym_ht=0_x"),
        ("default",),
        ("gnd",),
    )

    result = check_label_contract(mesh, strict=True)

    assert not result["passed"]
    assert result["autogenerated"]["materials"] == ["volume_2"]
    assert result["autogenerated"]["boundaries"] == ["Surface_7"]
    assert result["casefold_collisions"]["materials"] == [["Copper", "copper"]]
    assert "source and sink boundaries must appear as a pair" in result[
        "relational_issues"
    ]
    assert "symmetry axis x cannot carry both bn and ht labels" in result[
        "relational_issues"
    ]


def test_versioned_label_contract_enforces_required_and_allowed_labels():
    mesh = _LabelMesh(
        ("coil", "shield"),
        ("source", "sink", "sibc"),
        ("curve_1",),
        (),
    )
    contract = {
        "schema": LABEL_CONTRACT_SCHEMA,
        "application": "radia-ih",
        "strict_labels": True,
        "required": {
            "materials": ["coil", "air"],
            "boundaries": ["source", "sink"],
        },
        "allowed": {
            "materials": ["coil", "air"],
            "boundaries": ["source", "sink", "sibc"],
            "bboundaries": [],
        },
    }

    result = check_label_contract(mesh, contract=contract)

    assert not result["passed"]
    assert result["application"] == "radia-ih"
    assert result["missing"]["materials"] == ["air"]
    assert result["unexpected"]["materials"] == ["shield"]
    assert result["unexpected"]["bboundaries"] == ["curve_1"]


def test_label_contract_rejects_required_label_omitted_from_allowlist():
    mesh = _LabelMesh(("coil",), ("source", "sink"))
    contract = {
        "schema": LABEL_CONTRACT_SCHEMA,
        "required": {"materials": ["coil", "air"]},
        "allowed": {"materials": ["coil"]},
    }

    with pytest.raises(ValueError, match="absent from allowed materials: air"):
        check_label_contract(mesh, contract=contract)

    with pytest.raises(ValueError, match="must declare schema"):
        check_label_contract(mesh, contract={"required_materials": ["coil"]})


def test_cli_writes_versioned_json_report_and_requires_explicit_sidecar(
    tmp_path, capsys,
):
    _, vol_path = _save_box(tmp_path)
    report_path = tmp_path / "vol_check.json"
    contract_path = tmp_path / "labels.json"
    contract_path.write_text(
        json.dumps({
            "schema": LABEL_CONTRACT_SCHEMA,
            "application": "test-box",
            "strict_labels": True,
            "required": {
                "materials": ["body"],
                "boundaries": ["outer"],
            },
            "allowed": {
                "materials": ["body"],
                "boundaries": ["outer"],
            },
        }),
        encoding="utf-8",
    )

    exit_code = check_main([
        str(vol_path),
        "--contract", str(contract_path),
        "--report-json", str(report_path),
    ])

    assert exit_code == 0
    assert json.loads(report_path.read_text(encoding="utf-8"))["schema"] == REPORT_SCHEMA
    assert "All checks PASSED" in capsys.readouterr().out

    error_report = tmp_path / "vol_check_error.json"
    exit_code = check_main([
        str(vol_path),
        "--json", str(tmp_path / "missing.vol.json"),
        "--report-json", str(error_report),
    ])
    assert exit_code == 2
    assert "CAD reference not found" in capsys.readouterr().err
    assert json.loads(error_report.read_text(encoding="utf-8"))["passed"] is False


def test_cubit_toolbar_verifier_delegates_to_canonical_checker(tmp_path):
    mesh, vol_path = _save_box(tmp_path)
    sidecar = Path(str(vol_path) + ".json")
    sidecar.write_text(
        json.dumps({
            "materials": {"body": 1.0},
            "boundaries": {"outer": 6.0},
            "edges": {},
            "n_elements": mesh.ne,
            "n_points": mesh.nv,
            "order": mesh.GetCurveOrder(),
        }),
        encoding="utf-8",
    )

    result = VERIFY_MODULE.verify_vol(str(vol_path))
    updated = json.loads(sidecar.read_text(encoding="utf-8"))

    assert result["schema"] == REPORT_SCHEMA
    assert result["passed"]
    assert updated["ng_materials"]["body"] == pytest.approx(1.0)
    assert updated["ng_boundaries"]["outer"] == pytest.approx(6.0)
    assert updated["warnings"] == []
