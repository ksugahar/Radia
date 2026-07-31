"""Bind production motor and native Simulink artifacts into one proof."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT_DIR = (
    ROOT / "validation_test" / "radia_mcp" / "artifacts" / "annular_motor_dual_lane_v1"
)
OUTPUT = ARTIFACT_DIR / "production_replacement_proof.json"
REQUIRED_NATIVE_CAPABILITIES = {
    "periodic_angle_family_native_interpolation",
    "quadratic_torque_output",
    "persistent_native_state",
    "split_output_update_lifecycle",
    "custom_sim_state_roundtrip",
    "simulink_s_function_compile",
    "foreign_openmp_runtime_isolation",
}


def _load(name: str) -> dict:
    value = json.loads((ARTIFACT_DIR / name).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise RuntimeError(f"{name} must encode an object")
    return value


def _sha_file(name: str) -> str:
    return hashlib.sha256((ARTIFACT_DIR / name).read_bytes()).hexdigest()


def _sha_json(value: object) -> str:
    payload = json.dumps(
        value, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def main() -> int:
    manifest = _load("manifest.json")
    gate = _load("gate_result.json")
    age = _load("ngsolve_age.json")
    hdiv = _load("hdiv_mmm_hcurl_eddy_bubble.json")
    dual = _load("dual_lane.json")
    native = _load("native_motor_angle_family.json")

    identity = manifest.get("shared_mesh_material_identity")
    motor_checks = {
        "both_lanes_solved": (
            age.get("status") == "pass"
            and age.get("execution", {}).get("solver_status") == "solved"
            and hdiv.get("status") == "pass"
            and hdiv.get("execution", {}).get("checks", {}).get("passed") is True
        ),
        "shared_identity_matches": bool(identity)
        and age.get("shared_mesh_material_identity") == identity
        and hdiv.get("shared_mesh_material_identity") == identity
        and gate.get("shared_model_identity_matches") is True,
        "angle_grids_complete": (
            age.get("case_count") == 8
            and hdiv.get("case_count") == 8
            and len(age.get("execution", {}).get("torque_rows", [])) == 8
            and len(hdiv.get("execution", {}).get("motor_angle_rows", [])) == 8
        ),
        "torque_waveforms_nonconstant": dual.get("dual_lane_comparison", {})
        .get("checks", {})
        .get("both_torque_waveforms_nonconstant")
        is True,
        "correlation_gate_pass": (
            float(
                dual.get("dual_lane_comparison", {}).get(
                    "centered_normalized_correlation", 0.0
                )
            )
            >= 0.95
            and float(
                dual.get("dual_lane_comparison", {}).get(
                    "centered_normalized_rms", 1.0
                )
            )
            <= 0.12
        ),
        "dual_mcp_gate_pass": (
            gate.get("status") == "pass"
            and gate.get("validated_dual_solver_check") is True
            and gate.get("accepted_for_primary_dual_learning") is True
            and gate.get("accepted_for_mcp_learning") is True
        ),
    }
    test_names = {str(row.get("name", "")) for row in native.get("test_results", [])}
    native_hashes = [
        str(native.get("mex_sha256", "")),
        str(native.get("source_sha256", "")),
        str(native.get("setup_sha256", "")),
        str(native.get("generator_sha256", "")),
    ]
    native_checks = {
        "standalone_batch": native.get("execution_mode") == "standalone_matlab_batch",
        "all_tests_passed": (
            native.get("status") == "pass"
            and int(native.get("test_count", 0)) > 0
            and native.get("passed_count") == native.get("test_count")
            and native.get("failed_count") == 0
            and native.get("incomplete_count") == 0
        ),
        "mex_and_sources_hashed": all(
            len(value) == 64 and all(char in "0123456789abcdef" for char in value)
            for value in native_hashes
        ),
        "periodic_motor_handle_tested": any(
            name.endswith("/testNativePeriodicMotorAngleFamilyHandle")
            for name in test_names
        ),
        "simulink_compile_tested": any(
            name.endswith("/testNativePeriodicMotorAngleFamilyBlock")
            for name in test_names
        ),
        "foreign_openmp_runtime_isolated": (
            int(
                native.get(
                    "foreign_openmp_runtime_dirs_remaining_on_path_count", -1
                )
            )
            == 0
            and "foreign_openmp_runtime_isolation"
            in set(native.get("validated_capabilities", []))
        ),
    }
    if not all(motor_checks.values()):
        raise RuntimeError(f"motor dual-lane proof failed: {motor_checks}")
    if not all(native_checks.values()):
        raise RuntimeError(f"native motor proof failed: {native_checks}")
    if not REQUIRED_NATIVE_CAPABILITIES.issubset(native.get("validated_capabilities", [])):
        raise RuntimeError("native motor artifact is missing required capabilities")

    artifact = {
        "radia_version": str(manifest["radia_version"]),
        "schema": "radia.validation.production-replacement-proof.v1",
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_version": {
            "radia": str(manifest["radia_version"]),
            "matlab": str(native["matlab_release"]),
            "simulink": str(native["simulink_release"]),
        },
        "proofs": {
            "motor_dual_lane": {
                "status": "pass",
                "shared_model_identity_sha256": identity["aggregate_sha256"],
                "artifact_sha256_by_role": {
                    "manifest": _sha_file("manifest.json"),
                    "gate": _sha_file("gate_result.json"),
                    "ngsolve_age": _sha_file("ngsolve_age.json"),
                    "hdiv_mmm_hcurl_eddy_bubble": _sha_file(
                        "hdiv_mmm_hcurl_eddy_bubble.json"
                    ),
                    "dual_lane": _sha_file("dual_lane.json"),
                },
                "checks": motor_checks,
            },
            "native_motor_angle_family": {
                "status": "pass",
                "test_count": int(native["test_count"]),
                "artifact_sha256_by_role": {
                    "matlab": _sha_file("native_motor_angle_family.json")
                },
                "checks": native_checks,
            },
        },
        "status": "pass",
        "pass": True,
    }
    artifact["proof_payload_sha256"] = _sha_json(artifact)
    with OUTPUT.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(artifact, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": "pass",
                "output": str(OUTPUT.relative_to(ROOT)).replace("\\", "/"),
                "proof_payload_sha256": artifact["proof_payload_sha256"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
