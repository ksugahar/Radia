import hashlib
import json
from pathlib import Path

from radia_mcp.fem.uninstall_safety import (
    _validate_production_replacement_proof,
)
from radia_mcp.motor.triple_check_knowledge import (
    validate_motor_triple_check_artifact,
)


ARTIFACT_DIR = (
    Path(__file__).resolve().parents[3]
    / "validation_test"
    / "radia_mcp"
    / "artifacts"
    / "annular_motor_dual_lane_v1"
)


def _load(name: str) -> dict:
    path = ARTIFACT_DIR / name
    data = json.loads(path.read_text(encoding="utf-8"))
    assert next(iter(data)) == "radia_version"
    assert data["radia_version"]
    return data


def _text_sha256(path: Path) -> str:
    data = path.read_bytes()
    canonical = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(canonical).hexdigest()


def test_annular_motor_artifacts_keep_one_honest_harmonic_fixture_identity():
    age = _load("ngsolve_age.json")
    hdiv = _load("hdiv_mmm_hcurl_eddy_bubble.json")
    manifest = _load("manifest.json")

    assert age["shared_mesh_material_identity"] == hdiv[
        "shared_mesh_material_identity"
    ]
    geometry = manifest["physical_identity"]["geometry"]
    assert manifest["execution_environment"]["host_role"] == "compute"
    assert manifest["execution_environment"]["hostname"]
    assert geometry["family"] == "concentric_annular_harmonic_fixture"
    assert geometry["slot_geometry_resolved"] is False
    assert manifest["physical_identity"]["periodicity"]["role"] == (
        "nominal_symmetry_contract_not_resolved_slot_geometry"
    )
    assert age["execution"]["mesh_reused_all_angles"] is True
    assert age["execution"]["operator_reused_all_angles"] is True
    assert age["execution"]["factorization_reused_all_angles"] is True
    assert hdiv["execution"]["checks"]["passed"] is True
    assert hdiv["execution"]["checks"]["motor_angle_grid_complete"] is True


def test_annular_motor_dual_lane_artifact_remains_accepted_for_learning():
    packet = _load("dual_lane.json")
    comparison = packet["dual_lane_comparison"]
    result = validate_motor_triple_check_artifact(packet)

    assert comparison["torque_sign_mapping_age_from_hdiv"] == -1.0
    assert comparison["centered_normalized_correlation"] >= 0.95
    assert comparison["centered_normalized_rms"] <= 0.12
    assert comparison["checks"][
        "same_smooth_annular_harmonic_fixture_identity"
    ] is True
    assert all(comparison["checks"].values())
    assert result["validated_dual_solver_check"] is True
    assert result["shared_model_identity_matches"] is True
    assert result["accepted_for_primary_dual_learning"] is True
    assert result["accepted_for_mcp_learning"] is True


def test_annular_motor_saved_gate_matches_live_gate():
    packet = _load("dual_lane.json")
    saved = _load("gate_result.json")
    current = validate_motor_triple_check_artifact(packet)

    for key in (
        "status",
        "validated_dual_solver_check",
        "shared_model_identity_matches",
        "accepted_for_primary_dual_learning",
        "accepted_for_mcp_learning",
        "errors",
        "warnings",
    ):
        assert saved[key] == current[key]


def test_native_motor_angle_family_artifact_records_live_matlab_evidence():
    native = _load("native_motor_angle_family.json")
    manifest = _load("manifest.json")
    artifact_text = (ARTIFACT_DIR / "native_motor_angle_family.json").read_text(
        encoding="utf-8"
    )

    assert native["schema"] == "radia.validation.motor-angle-family-mex.v1"
    assert native["status"] == "pass"
    assert native["execution_mode"] == "standalone_matlab_batch"
    assert native["execution_environment"]["host_role"] in {
        "compute",
        "developer-smoke",
    }
    assert native["execution_environment"]["hostname"]
    assert native["matlab_release"] == "2026a"
    assert native["test_count"] == native["passed_count"]
    assert native["test_count"] >= 82
    assert native["failed_count"] == native["incomplete_count"] == 0
    assert isinstance(native["optimization_toolbox_available"], bool)
    assert native["foreign_openmp_runtime_dirs_remaining_on_path_count"] == 0
    assert len(native["mex_sha256"]) == 64
    assert len(native["source_sha256"]) == 64
    assert len(native["setup_sha256"]) == 64
    assert len(native["generator_sha256"]) == 64
    assert native["text_sha256_normalization"] == "newline-lf"
    assert "periodic_angle_family_native_interpolation" in native[
        "validated_capabilities"
    ]
    assert "simulink_s_function_compile" in native["validated_capabilities"]
    assert "split_output_update_lifecycle" in native["validated_capabilities"]
    assert "custom_sim_state_roundtrip" in native["validated_capabilities"]
    root = ARTIFACT_DIR.parents[3]
    for path_key, sha_key in (
        ("source_relative_path", "source_sha256"),
        ("setup_relative_path", "setup_sha256"),
        ("generator_relative_path", "generator_sha256"),
    ):
        path = root / native[path_key]
        assert native[sha_key] == _text_sha256(path)
    assert ":\\" not in artifact_text
    assert manifest["artifact_files"]["native_motor_angle_family"] == (
        "native_motor_angle_family.json"
    )


def test_production_replacement_proof_binds_motor_and_native_artifacts():
    path = ARTIFACT_DIR / "production_replacement_proof.json"
    raw = path.read_text(encoding="utf-8")
    proof = json.loads(raw)

    assert next(iter(proof)) == "radia_version"
    assert proof["status"] == "pass"
    assert proof["pass"] is True
    assert set(proof["proofs"]) == {
        "motor_dual_lane",
        "native_motor_angle_family",
    }
    assert all(proof["proofs"]["motor_dual_lane"]["checks"].values())
    assert all(proof["proofs"]["native_motor_angle_family"]["checks"].values())
    assert _validate_production_replacement_proof(
        {
            "artifact_json": raw,
            "artifact_sha256": hashlib.sha256(raw.encode("utf-8")).hexdigest(),
        }
    )
