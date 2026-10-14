import copy
import json

from radia_mcp.build123d.server import (
    build123d_heat_exchanger_source_recovery_gate,
    build123d_reflection_rotation_handoff_gate,
)


def handoff() -> dict:
    half_before = 100.0
    half_after = 101.0
    full = 202.0
    rows = [
        ("half_before_fillet", half_before, half_before * (1.0 - 1.0e-6), 1),
        ("half_after_fillet", half_after, half_after * (1.0 - 5.0e-7), 1),
        ("full_after_mirror", full, full * 0.75, 1),
        ("mirrored_half_alone", half_after, half_after * 1.02, 1),
        ("rotated_half_alone", half_after, half_after * (1.0 + 5.0e-7), 1),
        (
            "two_body_rotation_compound",
            full,
            half_after * (1.0 - 5.0e-7) + half_after * (1.0 + 5.0e-7),
            2,
        ),
    ]
    return {
        "volume_contract": {
            "unfilleted_volume_mm3": 200.0,
            "feature_delta_mm3": 2.0,
            "formula_final_volume_mm3": full,
            "build123d_final_volume_mm3": full,
        },
        "checkpoints": [
            {
                "name": name,
                "build123d_volume_mm3": native,
                "step_roundtrip_volume_mm3": native * (1.0 + 1.0e-13),
                "external_volume_mm3": external,
                "body_count": count,
                "volume_count": count,
                "step_sha256": str(index + 1) * 64,
            }
            for index, (name, native, external, count) in enumerate(rows)
        ],
        "disposition": "reflection_translation_bias_proper_rotation_compound_handoff_ready",
    }


def combined() -> dict:
    build_checks = {"source": True, "formula": True, "roundtrip": True}
    checkpoint_checks = {"source": True, "checkpoints": True, "roundtrip": True}
    external_rows = [{"name": row["name"]} for row in handoff()["checkpoints"]]
    return {
        "build": {
            "source_kind": "upstream_native_build123d_example_with_viewer_stub_only",
            "source_example": "examples/heat_exchanger.py",
            "upstream_tag": "v0.10.0",
            "upstream_git_blob_sha1": "a" * 40,
            "source_sha256_before": "b" * 64,
            "source_sha256_after": "b" * 64,
            "parameters_mm": {"tube_count": 148, "tube_location_count": 148},
            "topology": {
                "runtime_tube_count": 148,
                "source_comment_tube_count": 149,
                "source_comment_runtime_mismatch": True,
            },
            "timing": {
                "source_example_execute_s": 8.0,
                "formula_and_topology_inventory_s": 0.2,
                "step_export_s": 0.5,
                "step_roundtrip_s": 1.0,
            },
            "checks": build_checks,
        },
        "checkpoint_run": {
            "instrumentation": "runtime copy immediately before and after the upstream fillet call",
            "checks": checkpoint_checks,
        },
        "external_run": {
            "execution_mode": "python_api_headless_synchronous_commands",
            "headless_flags": ["-nographics", "-batch"],
            "gui_daemon_enabled": False,
            "rows": external_rows,
            "process": {
                "exit_code": 3,
                "unexpected_error_lines": [],
                "result_artifact_fresh": True,
                "owned_processes_remaining": 0,
                "process_exit_policy": "artifact_evidence_over_known_headless_diagnostics",
                "acceptable": True,
            },
        },
        "handoff": handoff(),
    }


def test_public_gate_accepts_failed_reflection_controls_and_rotation_recovery():
    result = json.loads(build123d_reflection_rotation_handoff_gate(json.dumps(handoff())))
    assert result["status"] == "ok"
    assert result["recommended_handoff"] == "two_body_rotation_compound"
    assert result["external_volume_relative_errors"]["full_after_mirror"] >= 0.2


def test_public_gate_rejects_biased_rotation_recovery():
    row = handoff()
    checkpoint = next(
        item for item in row["checkpoints"] if item["name"] == "two_body_rotation_compound"
    )
    checkpoint["external_volume_mm3"] *= 1.01
    result = json.loads(build123d_reflection_rotation_handoff_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["proper_rotation_compound_restores_volume"] is False


def test_public_gate_rejects_missing_reflection_negative_control():
    row = handoff()
    mirrored = next(
        item for item in row["checkpoints"] if item["name"] == "mirrored_half_alone"
    )
    mirrored["external_volume_mm3"] = mirrored["build123d_volume_mm3"]
    result = json.loads(build123d_reflection_rotation_handoff_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["reflection_negative_controls_expose_bias"] is False


def test_source_gate_accepts_bound_upstream_and_classified_headless_replay():
    result = json.loads(build123d_heat_exchanger_source_recovery_gate(json.dumps(combined())))
    assert result["status"] == "ok"
    assert result["runtime_tube_count"] == 148
    assert result["source_comment_tube_count"] == 149


def test_source_gate_rejects_comment_as_topology_and_stale_process():
    row = combined()
    row["build"]["parameters_mm"]["tube_count"] = 149
    row["external_run"]["process"]["result_artifact_fresh"] = False
    row["external_run"]["process"]["acceptable"] = False
    result = json.loads(build123d_heat_exchanger_source_recovery_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["runtime_topology_overrides_stale_comment"] is False
    assert result["checks"]["fresh_result_and_no_owned_process_leak"] is False


def test_source_gate_rejects_unexpected_error_and_failed_public_gate():
    row = combined()
    row["external_run"]["process"]["unexpected_error_lines"] = ["ERROR: unknown"]
    recovery = next(
        item
        for item in row["handoff"]["checkpoints"]
        if item["name"] == "rotated_half_alone"
    )
    recovery["external_volume_mm3"] *= 1.01
    result = json.loads(build123d_heat_exchanger_source_recovery_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["headless_external_cad_process_is_classified"] is False
    assert result["checks"]["independent_reflection_rotation_gate_passed"] is False


def test_server_returns_invalid_input_for_missing_checkpoint():
    row = handoff()
    row["checkpoints"] = row["checkpoints"][:-1]
    result = json.loads(build123d_reflection_rotation_handoff_gate(json.dumps(row)))
    assert result["status"] == "invalid_input"
    assert "six named" in result["error"]
