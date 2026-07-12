import copy
import json

from radia_mcp.cubit.server import (
    cubit_boundary_layer_candidate_gate,
    cubit_boundary_layer_journal_recovery_gate,
)


def candidate(height: float, minimum_scaled_jacobian: float, count: int) -> dict:
    return {
        "recovery": {
            "boundary_layer_parameters": {
                "height": height,
                "growth": 1.2,
                "layers": 4,
            },
            "one_surface_volume_pair_per_command": True,
        },
        "element_counts": {"hex": count, "tet": 0, "wedge": 0, "pyramid": 0},
        "quality": {"hex": {"scaled_jacobian": {"min": minimum_scaled_jacobian}}},
        "minimum_scaled_jacobian": minimum_scaled_jacobian,
        "shared_surface_count": 20,
        "shared_meshed_surface_count": 20,
        "volume_relative_error": 5.4e-6,
        "quarter_volume_relative_error": 5.4e-6,
        "rotational_copy_volume_relative_error": 2.0e-16,
        "export_inventory": {
            "volume_element_count": count,
            "volume_node_count_histogram": {"8": count},
        },
        "timing_breakdown_s": {
            "geometry_and_webcut": 0.13,
            "boundary_layer_recovery": 0.03,
            "sweep_mesh_and_rotation": 0.37,
            "netgen_export_and_parse": 0.28,
        },
    }


def ledger() -> list[dict]:
    return [candidate(2.0, -0.60, 7296), candidate(0.5, 0.836, 7980)]


def recovery_summary() -> dict:
    rows = ledger()
    rows[0]["candidate_status"] = "rejected"
    rows[1]["candidate_status"] = "accepted"
    return {
        "source_sha256": "a" * 64,
        "source_kind": "source_native_failed_journal_with_minimal_command_recovery",
        "source_failure": "missing height value plus combined surface-volume pairs",
        "execution_mode": "python_api_headless_synchronous_commands",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "candidates": rows,
        "process_runs": [
            {
                "process_exit_code": code,
                "startup_diagnostics": ["ERROR: Could not open file: C:/x/plugins"],
                "script_error_lines": ["Bad Quality Shear Hex Generated"] if code == 4 else [],
                "unexpected_script_error_lines": [],
                "expected_quality_rejection": code == 4,
                "result_artifact_fresh": True,
                "owned_processes_remaining": 0,
            }
            for code in (4, 3)
        ],
        "block_registered_netgen_export": True,
        "public_gate_status": "ok",
    }


def test_candidate_gate_selects_quality_not_element_count_or_geometry_only():
    result = json.loads(cubit_boundary_layer_candidate_gate(ledger()))
    assert result["status"] == "ok"
    assert result["selected_height"] == 0.5
    assert result["candidates"][0]["status"] == "rejected"
    assert result["candidates"][1]["status"] == "accepted"
    assert result["candidates"][0]["checks"]["analytic_volume_conserved"] is True
    assert result["candidates"][0]["checks"]["scaled_jacobian_above_threshold"] is False


def test_candidate_gate_rejects_a_ledger_without_a_discriminating_failure():
    rows = ledger()
    rows[0]["minimum_scaled_jacobian"] = 0.4
    rows[0]["quality"]["hex"]["scaled_jacobian"]["min"] = 0.4
    result = json.loads(cubit_boundary_layer_candidate_gate(rows))
    assert result["status"] == "needs_attention"
    assert result["checks"]["rejected_candidate_present"] is False


def test_journal_recovery_gate_accepts_three_parameter_pairwise_headless_replay():
    result = json.loads(cubit_boundary_layer_journal_recovery_gate(recovery_summary()))
    assert result["status"] == "ok"
    assert result["checks"]["exactly_three_boundary_layer_parameters"] is True
    assert result["checks"]["all_live_processes_closed_with_fresh_artifacts"] is True


def test_journal_recovery_gate_rejects_ambiguous_pairing_and_leaked_process():
    summary = recovery_summary()
    summary["candidates"][0]["recovery"]["one_surface_volume_pair_per_command"] = False
    summary["process_runs"][0]["owned_processes_remaining"] = 1
    result = json.loads(cubit_boundary_layer_journal_recovery_gate(summary))
    assert result["status"] == "needs_attention"
    assert result["checks"]["one_surface_volume_pair_per_command"] is False
    assert result["checks"]["all_live_processes_closed_with_fresh_artifacts"] is False
