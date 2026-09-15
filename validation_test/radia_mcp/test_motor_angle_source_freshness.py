"""Native acceptance: historical evidence must bind the candidate source.

This is deliberately separate from MCP package/schema tests. A failure means
new MATLAB/native validation evidence is required; never re-stamp old hashes.
"""

import hashlib
import json
from pathlib import Path

ARTIFACT_DIR = Path(__file__).resolve().parent / "artifacts/annular_motor_dual_lane_v1"


def _load(name):
    data = json.loads((ARTIFACT_DIR / name).read_text(encoding="utf-8"))
    assert next(iter(data)) == "radia_version"
    assert data["radia_version"]
    return data


def _text_sha256(path):
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


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
    provenance = native.get("native_build_provenance")
    assert isinstance(provenance, dict), (
        "Native evidence predates clean-build provenance; rerun MATLAB validation, "
        "do not restamp the historical artifact."
    )
    assert provenance["schema"] == "radia.native-build-provenance.v1"
    assert provenance["source_dirty"] is False
    assert len(provenance["source_commit"]) == 40
    assert provenance["binary_sha256"] == native["mex_sha256"]
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
