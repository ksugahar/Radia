from __future__ import annotations

from radia_mcp.cubit.cubit_v46_identity import validate_public_identity, validate_source_identity


PROMOTED_CASE_IDS = {
    "v46_public_mixed_hex_tet_orientation_degenerate_jacobian_partial_export_mismatch",
    "v46_public_unit_scale_coordinate_transform_node_order_nonfinite_quality_mismatch",
    "v46_source_tool_headless_journal_restart_command_status_partial_database_mismatch",
    "v46_source_tool_mesh_export_stream_truncation_checksum_process_exit_mismatch",
}


def _public() -> dict[str, object]:
    generation = "mixed-hex-tet-v46-901"
    quality_generation = "mesh-quality-v46-901"
    return {
        "mixed_hex_tet_orientation_degenerate_jacobian_partial_export_identity": {
            "generation": generation, "orientation_generation": generation,
            "jacobian_generation": generation, "partial_export_generation": generation,
            "result_generation": generation, "element_types": ["hex", "tet"],
            "result_element_types": ["hex", "tet"], "orientation": "positive",
            "result_orientation": "positive", "degenerate_jacobian_count": 0,
            "result_degenerate_jacobian_count": 0, "partial_export_state": "complete",
            "result_partial_export_state": "complete", "owner": "headless:mixed-v46",
            "accepted_owner": "headless:mixed-v46", "result_sha256": "5" * 64,
            "accepted_result_sha256": "5" * 64,
        },
        "unit_scale_coordinate_transform_node_order_nonfinite_quality_identity": {
            "generation": quality_generation, "unit_scale_generation": quality_generation,
            "coordinate_transform_generation": quality_generation, "node_order_generation": quality_generation,
            "quality_generation": quality_generation, "result_generation": quality_generation,
            "unit_name": "m", "result_unit_name": "m", "unit_scale_to_si": 1.0,
            "result_unit_scale_to_si": 1.0, "coordinate_transform": "global_cartesian",
            "result_coordinate_transform": "global_cartesian", "node_order": [1, 2, 3],
            "result_node_order": [1, 2, 3], "finite_quality_status": "finite",
            "result_finite_quality_status": "finite", "minimum_quality": 0.31,
            "result_minimum_quality": 0.31, "owner": "headless:quality-v46",
            "accepted_owner": "headless:quality-v46", "result_sha256": "6" * 64,
            "accepted_result_sha256": "6" * 64,
        },
    }


def _source() -> dict[str, object]:
    generation = "journal-restart-v46-901"
    export_generation = "stream-export-v46-901"
    commands = ["reset", "set units mm", "mesh volume 1"]
    statuses = ["success"] * 3
    checksum = "8" * 64
    return {
        "headless_journal_restart_command_status_partial_database_identity": {
            "generation": generation, "restart_generation": generation,
            "command_status_generation": generation, "partial_database_generation": generation,
            "result_generation": generation, "commands": commands, "result_commands": commands,
            "command_status": statuses, "result_command_status": statuses,
            "restart_state": "resumed_clean", "result_restart_state": "resumed_clean",
            "partial_database_state": "complete", "result_partial_database_state": "complete",
            "owner": "headless:journal-v46", "accepted_owner": "headless:journal-v46",
            "result_sha256": "7" * 64, "accepted_result_sha256": "7" * 64,
        },
        "mesh_export_stream_truncation_checksum_process_exit_identity": {
            "generation": export_generation, "stream_generation": export_generation,
            "checksum_generation": export_generation, "process_generation": export_generation,
            "result_generation": export_generation, "stream_truncated": False,
            "result_stream_truncated": False, "checksum_sha256": checksum,
            "result_checksum_sha256": checksum, "process_exit_code": 0,
            "result_process_exit_code": 0, "export_complete": True,
            "result_export_complete": True, "owner": "headless:export-v46",
            "accepted_owner": "headless:export-v46", "result_sha256": "9" * 64,
            "accepted_result_sha256": "9" * 64,
        },
    }


def test_v46_public_positive_is_accepted() -> None:
    result = validate_public_identity(_public())
    assert result["status"] == "ok"


def test_v46_public_orientation_and_quality_mutations_are_rejected() -> None:
    payload = _public()
    payload["mixed_hex_tet_orientation_degenerate_jacobian_partial_export_identity"]["result_orientation"] = "negative"
    payload["unit_scale_coordinate_transform_node_order_nonfinite_quality_identity"]["result_unit_scale_to_si"] = 1000.0
    result = validate_public_identity(payload)
    assert result["status"] == "needs_attention"


def test_v46_source_positive_is_accepted() -> None:
    assert validate_source_identity(_source())["status"] == "ok"


def test_v46_source_restart_and_export_mutations_are_rejected() -> None:
    payload = _source()
    payload["headless_journal_restart_command_status_partial_database_identity"]["result_restart_state"] = "fresh"
    payload["mesh_export_stream_truncation_checksum_process_exit_identity"]["result_stream_truncated"] = True
    assert validate_source_identity(payload)["status"] == "needs_attention"
