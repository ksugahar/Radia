"""Validation gates for Cubit mesh-carrying straight sweeps."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence


_STRAIGHT_MODES = {"vector", "perpendicular", "direction"}
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _normalized_command(value: object) -> str:
    return " ".join(str(value).strip().lower().split())


def _counts(value: object) -> dict[str, int]:
    row = _mapping(value, "element_counts")
    return {str(key).strip().lower(): int(count) for key, count in row.items()}


def _float_sequence(value: object, name: str) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    return [float(item) for item in value]


def cubit_mesh_carrying_straight_sweep_gate(
    summary: Mapping[str, object],
    *,
    min_scaled_jacobian: float = 0.2,
    volume_relative_tolerance: float = 1.0e-12,
) -> dict[str, object]:
    """Gate a mapped surface mesh carried through a straight all-hex sweep."""

    summary = _mapping(summary, "summary")
    if min_scaled_jacobian <= 0.0 or volume_relative_tolerance < 0.0:
        raise ValueError("quality threshold must be positive and tolerance nonnegative")

    mode = str(summary.get("sweep_mode", "")).strip().lower()
    command = _normalized_command(summary.get("command", ""))
    counts = _counts(summary.get("element_counts", {}))
    quality = _mapping(summary.get("quality", {}), "quality")
    scaled = _mapping(quality.get("scaled_jacobian", {}), "scaled_jacobian")
    gmsh = _mapping(summary.get("gmsh_export", {}), "gmsh_export")
    z_levels = _float_sequence(summary.get("z_levels", []), "z_levels")
    column_depths = [
        int(value)
        for value in _float_sequence(
            summary.get("xy_column_depths", []), "xy_column_depths"
        )
    ]

    source_quads = int(summary.get("source_quad_count", 0))
    source_nodes = int(summary.get("source_node_count", 0))
    source_xy = int(summary.get("source_xy_count", 0))
    intervals = int(summary.get("sweep_interval_count", 0))
    node_count = int(summary.get("node_count", 0))
    xy_columns = int(summary.get("xy_column_count", 0))
    complete_columns = int(summary.get("complete_xy_column_count", 0))
    total_volume = float(summary.get("cad_total_volume", math.nan))
    expected_volume = float(summary.get("expected_volume", math.nan))
    scaled_min = float(scaled.get("min", math.nan))
    expected_hexes = source_quads * intervals
    expected_nodes = source_nodes * (intervals + 1)
    volume_error = abs(total_volume - expected_volume) / max(
        abs(expected_volume), 1.0
    )
    z_strictly_increasing = len(z_levels) == intervals + 1 and all(
        right > left for left, right in zip(z_levels, z_levels[1:])
    )

    checks = {
        "straight_sweep_mode_supported": mode in _STRAIGHT_MODES,
        "command_is_surface_sweep": command.startswith("sweep surface "),
        "command_matches_mode": bool(mode) and mode in command.split(),
        "command_carries_existing_mesh": "include_mesh" in command.split(),
        "source_surface_mesh_present": source_quads > 0 and source_nodes > 0,
        "hex_only_volume_mesh": counts.get("hex", 0) > 0
        and all(counts.get(kind, 0) == 0 for kind in ("tet", "pyramid", "wedge")),
        "hex_count_is_source_quads_times_layers": counts.get("hex", 0)
        == expected_hexes,
        "strict_sweep_coordinate_layers": z_strictly_increasing,
        "node_lattice_count_conserved": node_count == expected_nodes,
        "source_xy_columns_conserved": source_xy == source_nodes
        and xy_columns == source_xy
        and complete_columns == source_xy,
        "uniform_column_depth": column_depths == [intervals + 1],
        "scaled_jacobian_above_threshold": scaled_min >= min_scaled_jacobian,
        "cad_volume_matches_expected": volume_error <= volume_relative_tolerance,
        "gmsh_ascii_v41": str(gmsh.get("mesh_format", "")) == "4.1"
        and gmsh.get("binary") is False,
        "gmsh_counts_match_live": int(gmsh.get("node_count", -1)) == node_count
        and int(gmsh.get("hex_count", -1)) == counts.get("hex", 0)
        and int(gmsh.get("other_volume_count", -1)) == 0,
    }
    return {
        "policy": "cubit_mesh_carrying_straight_sweep_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "source_quad_count": source_quads,
            "sweep_interval_count": intervals,
            "hex_count": counts.get("hex", 0),
            "expected_hex_count": expected_hexes,
            "node_count": node_count,
            "expected_node_count": expected_nodes,
            "scaled_jacobian_min": scaled_min,
            "volume_relative_error": volume_error,
        },
        "notes": [
            "Cubit reports surface mesh elements through the face entity family; identify quads by four-node face connectivity.",
            "For a straight include_mesh sweep, each mapped source quad becomes one hex per sweep interval.",
            "Omitting include_mesh can still create valid CAD volume while leaving the volume unmeshed.",
        ],
    }


def cubit_mesh_carrying_straight_sweep_source_replay_gate(
    summary: Mapping[str, object],
) -> dict[str, object]:
    """Gate official-command provenance and the headless no-mesh control."""

    summary = _mapping(summary, "summary")
    process = _mapping(summary.get("process", {}), "process")
    api = _mapping(summary.get("api_entity_contract", {}), "api_entity_contract")
    public_gate = _mapping(summary.get("public_gate", {}), "public_gate")
    negative = _mapping(summary.get("negative_control", {}), "negative_control")
    negative_counts = _counts(negative.get("element_counts", {}))
    timing = _mapping(summary.get("timing_breakdown_s", {}), "timing_breakdown_s")
    flags_raw = summary.get("headless_flags", [])
    if not isinstance(flags_raw, Sequence) or isinstance(flags_raw, (str, bytes)):
        raise ValueError("headless_flags must be a sequence")
    flags = {str(value).strip().lower() for value in flags_raw}
    source_contract = str(summary.get("source_contract", "")).strip().lower()
    source_sha256 = str(summary.get("source_sha256", "")).strip().lower()
    negative_command = _normalized_command(negative.get("command", ""))

    checks = {
        "official_help_command_seed": str(summary.get("source_kind", "")).strip()
        == "installed-official-help-command-with-synthetic-replay",
        "source_digest_recorded": bool(_SHA256_RE.fullmatch(source_sha256)),
        "source_contract_binds_existing_mesh": "include_mesh" in source_contract
        and "already meshed" in source_contract,
        "combined_journal_headless": str(summary.get("execution_mode", "")).strip()
        == "combined_journal_headless"
        and {"-nographics", "-batch"}.issubset(flags),
        "persistent_gui_not_started": summary.get("gui_daemon_enabled") is False,
        "fresh_classified_process": process.get("acceptable") is True
        and process.get("result_artifact_fresh") is True
        and not list(process.get("unexpected_error_lines", [])),
        "surface_mesh_uses_face_entities": str(api.get("surface_mesh_entity", ""))
        == "face"
        and int(api.get("quad_connectivity_size", 0)) == 4
        and int(api.get("quad_alias_count", -1)) == 0,
        "public_topology_gate_passed": public_gate.get("status") == "ok",
        "negative_command_omits_include_mesh": "include_mesh"
        not in negative_command.split(),
        "negative_still_creates_cad_volume": len(negative.get("volume_ids", []))
        == 1
        and float(negative.get("cad_total_volume", 0.0)) > 0.0,
        "negative_has_no_volume_elements": all(
            negative_counts.get(kind, 0) == 0
            for kind in ("hex", "tet", "pyramid", "wedge")
        ),
        "four_stage_timing_recorded": len(timing) == 4
        and all(float(value) >= 0.0 for value in timing.values()),
    }
    return {
        "policy": "cubit_mesh_carrying_straight_sweep_source_replay_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "notes": [
            "Use an already meshed planar surface before requesting include_mesh.",
            "Run Cubit through a combined journal in headless mode; do not rely on a persistent GUI daemon.",
            "Keep a no-include_mesh control because CAD creation alone does not prove that volume mesh was carried.",
        ],
    }
