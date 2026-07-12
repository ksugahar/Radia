"""Gates for structured all-hex meshes and source-journal replay evidence."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _relative(left: float, right: float) -> float:
	return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def cubit_structured_hex_lattice_gate(
	summary: Mapping[str, object],
	*,
	min_scaled_jacobian: float = 0.2,
	max_volume_relative_error: float = 1.0e-10,
) -> dict[str, object]:
	"""Gate interval-product counts, quality, and independent Gmsh closure."""

	if not isinstance(summary, Mapping):
		raise TypeError("summary must be a mapping")
	threshold = float(min_scaled_jacobian)
	volume_tolerance = float(max_volume_relative_error)
	if not math.isfinite(threshold) or threshold <= 0.0:
		raise ValueError("min_scaled_jacobian must be finite and positive")
	if not math.isfinite(volume_tolerance) or volume_tolerance < 0.0:
		raise ValueError("max_volume_relative_error must be finite and nonnegative")

	dimensions = [float(value) for value in summary.get("brick_dimensions") or []]
	intervals = [int(value) for value in summary.get("intervals") or []]
	if len(dimensions) != 3 or not all(math.isfinite(value) and value > 0.0 for value in dimensions):
		raise ValueError("brick_dimensions must contain three finite positive values")
	if len(intervals) != 3 or not all(value > 0 for value in intervals):
		raise ValueError("intervals must contain three positive integers")

	counts = summary.get("element_counts")
	quality = summary.get("quality")
	geometry = summary.get("geometry")
	header = summary.get("gmsh_header")
	independent = summary.get("gmsh_independent")
	if not all(isinstance(value, Mapping) for value in (counts, quality, geometry, header, independent)):
		raise ValueError("element_counts, quality, geometry, gmsh_header, and gmsh_independent must be mappings")

	expected_hexes = math.prod(intervals)
	expected_nodes = math.prod(value + 1 for value in intervals)
	hex_count = int(counts.get("hex", 0))
	non_hex_count = sum(int(counts.get(kind, 0)) for kind in ("tet", "pyramid", "wedge"))
	node_count = int(summary.get("node_count", 0))
	scaled = quality.get("scaled_jacobian") or {}
	shape = quality.get("shape") or {}
	if not isinstance(scaled, Mapping) or not isinstance(shape, Mapping):
		raise ValueError("quality rows must be mappings")
	cad_volume = float(geometry.get("cad_volume", math.nan))
	analytic_volume = math.prod(dimensions)
	independent_volume = float(independent.get("integrated_volume", math.nan))
	checks = {
		"hex_count_is_interval_product": hex_count == expected_hexes,
		"node_count_is_structured_lattice": node_count == expected_nodes,
		"all_volume_elements_are_hex": hex_count > 0 and non_hex_count == 0,
		"all_hexes_have_eight_nodes": list(summary.get("connectivity_sizes") or []) == [8],
		"scaled_jacobian_covers_all_hexes": int(scaled.get("count", 0)) == hex_count
		and float(scaled.get("min", -math.inf)) >= threshold,
		"shape_covers_all_hexes": int(shape.get("count", 0)) == hex_count
		and float(shape.get("min", -math.inf)) > 0.0,
		"cad_volume_matches_brick": math.isfinite(cad_volume)
		and _relative(cad_volume, analytic_volume) <= volume_tolerance,
		"gmsh_ascii_v41_inventory_valid": header.get("status") == "ok"
		and header.get("mesh_format") == "4.1"
		and header.get("binary") is False,
		"gmsh_counts_match_live_inventory": int(independent.get("hex_count", -1)) == hex_count
		and int(independent.get("node_count", -1)) == node_count
		and int(independent.get("non_hex_volume_count", -1)) == 0,
		"gmsh_integrated_volume_matches_cad": math.isfinite(independent_volume)
		and _relative(independent_volume, cad_volume) <= volume_tolerance,
	}
	issues = [name for name, ok in checks.items() if not ok]
	return {
		"policy": "cubit_structured_hex_lattice_gate_v1",
		"status": "ok" if not issues else "needs_attention",
		"checks": checks,
		"issues": issues,
		"metrics": {
			"expected_hex_count": expected_hexes,
			"expected_node_count": expected_nodes,
			"observed_hex_count": hex_count,
			"observed_node_count": node_count,
			"minimum_scaled_jacobian": float(scaled.get("min", math.nan)),
			"cad_volume_relative_error": _relative(cad_volume, analytic_volume),
			"gmsh_volume_relative_error": _relative(independent_volume, cad_volume),
		},
		"lesson": (
			"A structured all-hex mesh is reproducible when element and node counts close "
			"against the interval lattice, quality covers every cell, and an independently "
			"parsed export closes geometry volume. A large element count alone is insufficient."
		),
	}


def cubit_structured_hex_source_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
	"""Gate source-journal identity, headless execution, and semantic exit policy."""

	if not isinstance(summary, Mapping):
		raise TypeError("summary must be a mapping")
	commands = [" ".join(str(value).strip().lower().split()) for value in summary.get("source_commands") or []]
	required = ["reset", "brick x 10 y 10 z 10", "volume 1 interval 40", "mesh vol all"]
	process = summary.get("process") or {}
	license_probe = summary.get("license_probe") or {}
	if not isinstance(process, Mapping) or not isinstance(license_probe, Mapping):
		raise ValueError("process and license_probe must be mappings")
	error_categories = set(str(value) for value in process.get("error_categories") or [])
	allowed_categories = {"plugin_argument_diagnostic", "session_error_summary"}
	timing = summary.get("timing_breakdown_s") or {}
	public_gate = summary.get("public_gate") or {}
	checks = {
		"source_native_journal_recorded": summary.get("source_native_journal") is True,
		"source_digest_is_sha256": len(str(summary.get("source_sha256", ""))) == 64,
		"source_command_contract_preserved": commands[:4] == required,
		"headless_combined_replay_recorded": summary.get("execution_mode")
		== "headless_combined_journal_then_python_inventory"
		and {"-nographics", "-batch"}.issubset(set(summary.get("headless_flags") or []))
		and summary.get("gui_daemon_enabled") is False,
		"license_status_and_operational_evidence_separated": bool(str(license_probe.get("status", "")).strip())
		and license_probe.get("mesh_command_completed") is True,
		"process_exit_semantically_classified": int(process.get("exit_code", -1)) in {0, 3, 4}
		and error_categories <= allowed_categories
		and process.get("unexpected_error_lines") == []
		and process.get("acceptable") is True,
		"fresh_artifact_and_no_owned_process_leak": process.get("result_artifact_fresh") is True
		and int(process.get("owned_processes_remaining", -1)) == 0,
		"exactly_four_timing_stages_recorded": isinstance(timing, Mapping)
		and len(timing) == 4
		and all(math.isfinite(float(value)) and float(value) >= 0.0 for value in timing.values()),
		"structured_lattice_gate_passed": isinstance(public_gate, Mapping)
		and public_gate.get("policy") == "cubit_structured_hex_lattice_gate_v1"
		and public_gate.get("status") == "ok",
	}
	issues = [name for name, ok in checks.items() if not ok]
	return {
		"policy": "cubit_structured_hex_source_replay_gate_v1",
		"status": "ok" if not issues else "needs_attention",
		"checks": checks,
		"issues": issues,
		"license_status": license_probe.get("status"),
		"process_exit_code": process.get("exit_code"),
		"error_categories": sorted(error_categories),
		"notes": [
			"Do not infer mesh success or failure from a license line alone; require completed commands and fresh numerical artifacts.",
			"A nonzero launcher exit is acceptable only for named diagnostics with no unexpected errors, a fresh passing result, and no owned process leak.",
		],
	}
