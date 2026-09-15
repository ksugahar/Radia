"""Neutral v46 Coreform/Cubit mesh and headless replay identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _same(row: Mapping[str, object], *names: str) -> bool:
    return all(row.get(f"result_{name}") == row.get(name) for name in names)


def _mixed_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation", "")).strip()
    elements = row.get("element_types")
    return (
        bool(generation)
        and row.get("orientation_generation") == generation
        and row.get("jacobian_generation") == generation
        and row.get("partial_export_generation") == generation
        and row.get("result_generation") == generation
        and isinstance(elements, list)
        and elements == row.get("result_element_types")
        and set(elements) >= {"hex", "tet"}
        and row.get("orientation") == row.get("result_orientation") == "positive"
        and row.get("degenerate_jacobian_count") == row.get("result_degenerate_jacobian_count") == 0
        and row.get("partial_export_state") == row.get("result_partial_export_state") == "complete"
        and bool(str(row.get("owner") or ""))
        and row.get("accepted_owner") == row.get("owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _quality_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation", "")).strip()
    scale = row.get("unit_scale_to_si")
    order = row.get("node_order")
    minimum = row.get("minimum_quality")
    return (
        bool(generation)
        and row.get("unit_scale_generation") == generation
        and row.get("coordinate_transform_generation") == generation
        and row.get("node_order_generation") == generation
        and row.get("quality_generation") == generation
        and row.get("result_generation") == generation
        and isinstance(scale, (int, float))
        and math.isfinite(float(scale))
        and float(scale) > 0.0
        and row.get("result_unit_scale_to_si") == scale
        and row.get("coordinate_transform") == row.get("result_coordinate_transform") == "global_cartesian"
        and isinstance(order, list)
        and order
        and order == row.get("result_node_order")
        and row.get("finite_quality_status") == row.get("result_finite_quality_status") == "finite"
        and isinstance(minimum, (int, float))
        and math.isfinite(float(minimum))
        and float(minimum) >= 0.0
        and row.get("result_minimum_quality") == minimum
        and _same(row, "unit_name", "coordinate_transform")
        and bool(str(row.get("owner") or ""))
        and row.get("accepted_owner") == row.get("owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _journal_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation", "")).strip()
    commands = row.get("commands")
    statuses = row.get("command_status")
    return (
        bool(generation)
        and row.get("restart_generation") == generation
        and row.get("command_status_generation") == generation
        and row.get("partial_database_generation") == generation
        and row.get("result_generation") == generation
        and isinstance(commands, list)
        and commands == row.get("result_commands")
        and isinstance(statuses, list)
        and statuses
        and all(item == "success" for item in statuses)
        and statuses == row.get("result_command_status")
        and row.get("restart_state") == row.get("result_restart_state") == "resumed_clean"
        and row.get("partial_database_state") == row.get("result_partial_database_state") == "complete"
        and bool(str(row.get("owner") or ""))
        and row.get("accepted_owner") == row.get("owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _export_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation", "")).strip()
    checksum = str(row.get("checksum_sha256") or "")
    return (
        bool(generation)
        and row.get("stream_generation") == generation
        and row.get("checksum_generation") == generation
        and row.get("process_generation") == generation
        and row.get("result_generation") == generation
        and row.get("stream_truncated") is row.get("result_stream_truncated") is False
        and _digest(checksum)
        and row.get("result_checksum_sha256") == checksum
        and row.get("process_exit_code") == row.get("result_process_exit_code") == 0
        and row.get("export_complete") == row.get("result_export_complete") is True
        and bool(str(row.get("owner") or ""))
        and row.get("accepted_owner") == row.get("owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    records = {
        "v46_mixed_hex_tet_orientation_partial_export": payload.get("mixed_hex_tet_orientation_degenerate_jacobian_partial_export_identity"),
        "v46_unit_transform_node_quality": payload.get("unit_scale_coordinate_transform_node_order_nonfinite_quality_identity"),
    }
    checks = {name: (isinstance(row, Mapping) and (_mixed_ok(row) if "mixed" in name else _quality_ok(row))) for name, row in records.items() if row is not None}
    if not checks:
        return {}
    return {"policy": "cubit_v46_public_identity_v1", "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, ok in checks.items() if not ok]}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    records = {
        "v46_headless_journal_restart_database": payload.get("headless_journal_restart_command_status_partial_database_identity"),
        "v46_export_stream_checksum_exit": payload.get("mesh_export_stream_truncation_checksum_process_exit_identity"),
    }
    checks = {name: (isinstance(row, Mapping) and (_journal_ok(row) if "journal" in name else _export_ok(row))) for name, row in records.items() if row is not None}
    if not checks:
        return {}
    return {"policy": "cubit_v46_source_identity_v1", "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, ok in checks.items() if not ok]}
