import copy
import json

from radia_mcp.cubit.server import (
    cubit_helical_partition_mesh_gate,
    cubit_source_journal_replay_gate,
)


def good():
    return {
        "source_sha256": "a" * 64,
        "source_kind": "source_native_local_journal_with_filtered_webcut_replay",
        "execution_mode": "python_api_headless_synchronous_commands",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "playback_async_probe": {
            "unsafe_zero_entity_artifact_rejected": True,
            "replacement": "synchronous cubit.cmd sequence",
        },
        "operations": {"webcut_intersection_selection_counts": [253] * 15},
        "volume_count": 4048,
        "element_counts": {"hex": 65706, "pyramid": 0, "wedge": 0, "tet": 0},
        "quality": {"hex": {"scaled_jacobian": {"min": 0.73}}},
        "webcut_volume_relative_drift": 4.0e-6,
        "analytic_volume_relative_error": 8.0e-6,
        "shared_surface_count": 3795,
        "shared_meshed_surface_count": 3795,
        "export_inventory": {
            "volume_elements": 65706,
            "volume_kind_counts": {"hex": 65706},
            "routing_hint": "cubit_hex_or_mixed_path",
            "points": 167393,
        },
        "timing_breakdown_s": {"geometry_boolean": 402.0, "filtered_webcuts": 108.0, "local_merge": 53.0, "mesh_and_export": 251.0},
        "process_exit_code": 3,
        "startup_diagnostics": [
            "ERROR: Could not open file: C:/x/plugins",
            "ERROR: Could not open file: -commandplugindir",
            "ERROR: Could not open file: -nojournal",
        ],
        "script_error_lines": [],
        "result_artifact_fresh": True,
        "owned_processes_remaining": 0,
        "public_gate_status": "ok",
        "expected_public_gate_status": "ok",
    }


def test_accepts_solver_ready_helical_partition_and_safe_replay():
    row = good()
    assert json.loads(cubit_helical_partition_mesh_gate(row))["status"] == "ok"
    assert json.loads(cubit_source_journal_replay_gate(row))["status"] == "ok"


def test_rejects_live_failure_triad_without_overclaiming():
    row = good()
    row["quality"]["hex"]["scaled_jacobian"]["min"] = -0.6978
    row["webcut_volume_relative_drift"] = 5.2e-5
    row["shared_surface_count"] = 0
    row["shared_meshed_surface_count"] = 0
    row["export_inventory"]["volume_elements"] = 0
    row["export_inventory"]["volume_kind_counts"] = {}
    row["export_inventory"]["routing_hint"] = "inspect_before_solver_import"
    result = json.loads(cubit_helical_partition_mesh_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["scaled_jacobian_acceptable"] is False
    assert result["checks"]["export_volume_elements_present"] is False


def test_replay_gate_accepts_expected_rejection_but_rejects_early_artifact():
    row = good()
    row["public_gate_status"] = "needs_attention"
    row["expected_public_gate_status"] = "needs_attention"
    assert json.loads(cubit_source_journal_replay_gate(row))["status"] == "ok"
    bad = copy.deepcopy(row)
    bad["playback_async_probe"]["unsafe_zero_entity_artifact_rejected"] = False
    result = json.loads(cubit_source_journal_replay_gate(bad))
    assert result["status"] == "needs_attention"
    assert result["checks"]["early_playback_artifact_rejected"] is False
