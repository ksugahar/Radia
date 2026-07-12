from __future__ import annotations

import copy
import json

from radia_mcp.cubit.server import (
	cubit_structured_hex_lattice_gate,
	cubit_structured_hex_source_replay_gate,
)


def lattice_summary() -> dict:
	return {
		"brick_dimensions": [4.0, 3.0, 2.0],
		"intervals": [4, 3, 2],
		"element_counts": {"hex": 24, "tet": 0, "pyramid": 0, "wedge": 0},
		"node_count": 60,
		"connectivity_sizes": [8],
		"quality": {
			"scaled_jacobian": {"count": 24, "min": 0.99},
			"shape": {"count": 24, "min": 0.98},
		},
		"geometry": {"cad_volume": 24.0},
		"gmsh_header": {"status": "ok", "mesh_format": "4.1", "binary": False},
		"gmsh_independent": {
			"hex_count": 24,
			"node_count": 60,
			"non_hex_volume_count": 0,
			"integrated_volume": 24.0,
		},
	}


def call_lattice(summary: dict) -> dict:
	return json.loads(cubit_structured_hex_lattice_gate(summary))


def source_summary(public_gate: dict) -> dict:
	return {
		"source_native_journal": True,
		"source_sha256": "a" * 64,
		"source_commands": [
			"reset",
			"brick x 10 y 10 z 10",
			"volume 1 interval 40",
			"mesh vol all",
		],
		"execution_mode": "headless_combined_journal_then_python_inventory",
		"headless_flags": ["-nographics", "-batch"],
		"gui_daemon_enabled": False,
		"license_probe": {"status": "valid", "mesh_command_completed": True},
		"process": {
			"exit_code": 3,
			"error_categories": ["plugin_argument_diagnostic", "session_error_summary"],
			"unexpected_error_lines": [],
			"acceptable": True,
			"result_artifact_fresh": True,
			"owned_processes_remaining": 0,
		},
		"timing_breakdown_s": {"geometry": 0.1, "mesh": 0.2, "export": 0.3, "review": 0.4},
		"public_gate": public_gate,
	}


def test_accepts_interval_lattice_quality_and_independent_volume():
	result = call_lattice(lattice_summary())
	assert result["status"] == "ok"
	assert result["metrics"]["expected_hex_count"] == 24
	assert result["metrics"]["expected_node_count"] == 60


def test_rejects_large_count_without_lattice_or_export_closure():
	summary = copy.deepcopy(lattice_summary())
	summary["element_counts"]["hex"] = 25
	summary["gmsh_independent"]["integrated_volume"] = 23.0
	result = call_lattice(summary)
	assert result["status"] == "needs_attention"
	assert set(result["issues"]) >= {
		"hex_count_is_interval_product",
		"gmsh_counts_match_live_inventory",
		"gmsh_integrated_volume_matches_cad",
	}


def test_source_gate_accepts_classified_nonzero_headless_exit():
	result = json.loads(cubit_structured_hex_source_replay_gate(source_summary(call_lattice(lattice_summary()))))
	assert result["status"] == "ok"


def test_source_gate_rejects_license_only_claim_and_stale_artifact():
	summary = source_summary(call_lattice(lattice_summary()))
	summary["license_probe"]["mesh_command_completed"] = False
	summary["process"]["result_artifact_fresh"] = False
	result = json.loads(cubit_structured_hex_source_replay_gate(summary))
	assert result["status"] == "needs_attention"
	assert set(result["issues"]) >= {
		"license_status_and_operational_evidence_separated",
		"fresh_artifact_and_no_owned_process_leak",
	}
