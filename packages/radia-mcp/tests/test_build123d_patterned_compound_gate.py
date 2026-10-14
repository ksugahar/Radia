import copy
import json

from radia_mcp.build123d.server import (
    build123d_patterned_compound_translation_gate,
    build123d_wrap_faces_rotational_source_replay_gate,
)


def summary() -> dict:
    group = {
        "copy_count": 180,
        "valid_copy_count": 180,
        "copy_to_prototype_max_relative_error": 2.0e-15,
        "radial_center_relative_spread": 3.0e-16,
        "x_center_spread": 0.0,
    }
    external_row = {
        "mode": "noheal",
        "body_count": 1081,
        "volume_count": 1081,
        "surface_count": 6877,
        "curve_count": 14114,
        "vertex_count": 9397,
        "total_volume_relative_error": 0.0081238835,
        "tire_volume_relative_error": 0.0087624337,
        "tread_total_volume_relative_error": 0.0003440066,
        "maximum_sorted_body_volume_relative_error": 0.0087624337,
    }
    return {
        "authoring": {
            "pass": True,
            "version": "0.10.0",
            "upstream_tag": "v0.10.0",
            "upstream_commit": "a" * 40,
            "source_kind": "upstream_native_build123d_example_with_viewer_stub_only",
            "source_example": "examples/bicycle_tire.py",
            "source_sha256_before": "b" * 64,
            "source_sha256_after": "b" * 64,
            "source_operations": ["Bezier", "revolve", "wrap_faces", "thicken", "Rot"],
            "viewer_stubbed": True,
            "prototype_count": 6,
            "rotation_count_per_prototype": 180,
            "rotation_step_degrees": 2,
            "tread_solid_count": 1080,
            "authoring_solid_count": 1081,
            "prototype_groups": [{**group, "prototype_index": index} for index in range(6)],
            "step": {
                "sha256": "c" * 64,
                "roundtrip_solid_count": 1081,
                "roundtrip_volume_relative_error": 3.72e-8,
            },
        },
        "external": {
            "execution_mode": "python_api_headless_synchronous_commands",
            "headless_flags": ["-nographics", "-batch"],
            "gui_daemon_enabled": False,
            "step_sha256": "c" * 64,
            "external_solver_ready": False,
            "disposition": "translation_bias_detected_not_solver_ready",
            "rows": [external_row, {**external_row, "mode": "heal"}],
        },
        "external_process": {
            "exit_code": 3,
            "startup_diagnostics": [
                "ERROR: Could not open file: C:/x/plugins",
                "ERROR: Could not open file: -commandplugindir",
                "ERROR: Could not open file: -nojournal",
            ],
            "script_error_lines": [],
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
        "timing_breakdown_s": {
            "source_execution": 15.1,
            "pattern_inventory_and_step": 17.0,
            "external_heal_noheal": 33.3,
            "mcp_verification": 7.0,
        },
    }


def call(tool, row):
    return json.loads(tool(json.dumps(row)))


def test_translation_gate_accepts_diagnosis_but_rejects_solver_ready_claim():
    result = call(build123d_patterned_compound_translation_gate, summary())
    assert result["status"] == "ok"
    assert result["diagnosis"] == "dominant_curved_body_translation_bias"
    assert result["solver_ready"] is False
    assert result["checks"]["all_patterned_bodies_preserved"] is True


def test_translation_gate_rejects_lost_body_and_false_solver_ready_claim():
    row = summary()
    row["external"]["rows"][0]["body_count"] = 1080
    row["external"]["external_solver_ready"] = True
    result = call(build123d_patterned_compound_translation_gate, row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_patterned_bodies_preserved"] is False
    assert result["checks"]["external_solver_handoff_is_rejected"] is False


def test_translation_gate_rejects_unlocalized_or_hidden_external_bias():
    row = summary()
    for external in row["external"]["rows"]:
        external["total_volume_relative_error"] = 1.0e-7
        external["tread_total_volume_relative_error"] = 0.002
    result = call(build123d_patterned_compound_translation_gate, row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["external_translation_bias_is_material"] is False
    assert result["checks"]["bias_is_localized_to_dominant_curved_body"] is False


def test_source_gate_accepts_exact_six_by_180_upstream_replay():
    result = call(build123d_wrap_faces_rotational_source_replay_gate, summary())
    assert result["status"] == "ok"
    assert result["prototype_count"] == 6
    assert result["copy_count"] == 1080
    assert result["solver_ready"] is False


def test_source_gate_rejects_source_mutation_and_pattern_drift():
    row = summary()
    row["authoring"]["source_sha256_after"] = "d" * 64
    row["authoring"]["source_operations"].remove("wrap_faces")
    row["authoring"]["prototype_groups"][0]["valid_copy_count"] = 179
    result = call(build123d_wrap_faces_rotational_source_replay_gate, row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_preserved_with_viewer_stub_only"] is False
    assert result["checks"]["wrap_thicken_rotation_sequence_recorded"] is False
    assert result["checks"]["every_rotational_copy_is_valid_and_invariant"] is False


def test_source_gate_rejects_exit_code_only_and_stale_output():
    row = summary()
    row["external_process"]["startup_diagnostics"] = []
    row["external_process"]["result_artifact_fresh"] = False
    row["external_process"]["owned_processes_remaining"] = 1
    result = call(build123d_wrap_faces_rotational_source_replay_gate, row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["startup_only_nonzero_exit_explained"] is False
    assert result["checks"]["fresh_artifact_and_no_owned_process_leak"] is False


def test_server_wrappers_report_invalid_input():
    result = json.loads(build123d_patterned_compound_translation_gate("{}"))
    assert result["status"] == "invalid_input"
