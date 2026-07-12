import copy
import json

from radia_mcp.build123d.server import (
    build123d_jointed_assembly_source_replay_gate,
    build123d_jointed_assembly_step_closure_gate,
)


def summary():
    return {
        "components": [
            {
                "name": "body",
                "native_valid": True,
                "native_solid_count": 1,
                "native_volume_mm3": 6000.0,
                "tessellated_volume_mm3": 5999.7,
                "joint_names": ["body_top"],
                "self_roundtrip": {"solid_count": 0, "total_volume_mm3": 0.0},
                "expected_disposition": "reject_solid_closure_loss",
            },
            {
                "name": "cab",
                "native_valid": True,
                "native_solid_count": 1,
                "native_volume_mm3": 500.0,
                "tessellated_volume_mm3": 499.95,
                "joint_names": ["cab_base"],
                "self_roundtrip": {"solid_count": 1, "total_volume_mm3": 500.0001},
                "expected_disposition": "portable_control",
            },
        ],
        "external_rows": [
            {"name": "body.step", "volume_count": 1, "total_volume_mm3": 0.0},
            {"name": "cab.step", "volume_count": 1, "total_volume_mm3": 500.005},
            {"name": "assembly.step", "volume_count": 2, "total_volume_mm3": 500.005},
        ],
        "assembly": {
            "step_name": "assembly.step",
            "native_total_volume_mm3": 6500.0,
            "self_roundtrip": {"solid_count": 1, "total_volume_mm3": 500.0001},
        },
    }


def replay():
    row = summary()
    return {
        "source_kind": "upstream_source_native_example_with_display_stub_only",
        "source_sha256": "a" * 64,
        "source_url": "https://example.invalid/project/blob/v0.10.0/examples/model.py",
        "source_preserved": True,
        "display_stubbed_only": True,
        "components": row["components"],
        "joint_connections": [{"from": "body_top", "to": "cab_base", "kind": "rigid"}],
        "external_execution": {
            "mode": "python_api_headless_synchronous_commands",
            "headless_flags": ["-nographics", "-batch"],
            "gui_daemon_enabled": False,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
        "diagnosis_gate_status": "ok",
        "diagnosis": "component_solid_closure_loss",
        "solver_ready": False,
        "timing_breakdown_s": {
            "step_roundtrips": 2.1,
            "source_model": 1.1,
            "external_import": 0.6,
            "evidence_merge": 0.01,
        },
    }


def test_component_gate_locates_missing_body_with_portable_cab_control():
    result = json.loads(build123d_jointed_assembly_step_closure_gate(json.dumps(summary())))
    assert result["status"] == "ok"
    assert result["diagnosis"] == "component_solid_closure_loss"
    assert result["solver_ready"] is False
    assert result["components"][0]["observed_disposition"] == "reject_solid_closure_loss"
    assert result["components"][1]["observed_disposition"] == "portable_control"


def test_component_gate_rejects_missing_tessellated_support_or_wrong_disposition():
    row = summary()
    row["components"][0]["tessellated_volume_mm3"] = 4000.0
    row["components"][0]["expected_disposition"] = "portable_control"
    result = json.loads(build123d_jointed_assembly_step_closure_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["native_brep_supported_by_tessellated_volume"] is False
    assert result["checks"]["solid_closure_loss_component_present"] is False


def test_source_replay_accepts_joint_graph_and_headless_diagnosis():
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(replay())))
    assert result["status"] == "ok"
    assert result["checks"]["joint_connection_endpoints_resolve"] is True


def test_source_replay_rejects_unresolved_joint_and_gui_or_process_leak():
    row = replay()
    row["joint_connections"][0]["to"] = "missing_joint"
    row["external_execution"]["gui_daemon_enabled"] = True
    row["external_execution"]["owned_processes_remaining"] = 1
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["joint_connection_endpoints_resolve"] is False
    assert result["checks"]["headless_external_cad_replay"] is False
    assert result["checks"]["fresh_result_and_owned_process_cleanup"] is False
