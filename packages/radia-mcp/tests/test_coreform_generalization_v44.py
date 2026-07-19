from __future__ import annotations

from radia_mcp.cubit.mixed_transition_gate import (
    cubit_conformal_hex_pyramid_tet_interface_gate,
    cubit_mixed_transition_source_gate,
)


_PUBLIC_CASE = (
    "v44_public_hex_sweep_periodic_interface_quality_jacobian_block_sideset_export_mismatch"
)
_MIXED_CASE = (
    "v44_public_tet_hex_transition_pyramid_orientation_volume_area_quality_export_mismatch"
)
_JOURNAL_CASE = (
    "v44_source_journal_replay_command_order_session_units_geometry_generation_database_owner_mismatch"
)
_QUALITY_CASE = (
    "v44_source_mesh_quality_metric_reference_element_dimension_block_export_owner_mismatch"
)


def _summary() -> dict:
    generation = "periodic-hex-v44-731"
    return {
        "hex_sweep_periodic_interface_quality_jacobian_block_sideset_export_generation_identity": {
            "sweep_generation": generation,
            "periodic_generation": generation,
            "interface_generation": generation,
            "quality_generation": generation,
            "jacobian_generation": generation,
            "block_generation": generation,
            "sideset_generation": generation,
            "export_generation": generation,
            "result_generation": generation,
            "source_surface_id": 11,
            "result_source_surface_id": 11,
            "target_surface_id": 21,
            "result_target_surface_id": 21,
            "source_node_count": 25,
            "result_source_node_count": 25,
            "target_node_count": 25,
            "result_target_node_count": 25,
            "paired_source_node_ids": [101, 102, 103, 104],
            "result_paired_source_node_ids": [101, 102, 103, 104],
            "paired_target_node_ids": [201, 202, 203, 204],
            "result_paired_target_node_ids": [201, 202, 203, 204],
            "coordinate_frame": "global_cartesian",
            "result_coordinate_frame": "global_cartesian",
            "periodic_transform_matrix": [[1.0, 0.0], [0.0, 1.0]],
            "result_periodic_transform_matrix": [[1.0, 0.0], [0.0, 1.0]],
            "interface_normal_dot": -1.0,
            "result_interface_normal_dot": -1.0,
            "minimum_scaled_jacobian": 0.36,
            "result_minimum_scaled_jacobian": 0.36,
            "minimum_allowed_scaled_jacobian": 0.20,
            "block_membership": {"block:periodic_hex": [1]},
            "result_block_membership": {"block:periodic_hex": [1]},
            "sideset_membership": {"sideset:periodic": [11, 21]},
            "result_sideset_membership": {"sideset:periodic": [11, 21]},
            "mesh_owner": "headless:periodic-hex-731",
            "result_mesh_owner": "headless:periodic-hex-731",
            "mesh_export_sha256": "a" * 64,
            "accepted_mesh_export_sha256": "a" * 64,
        },
        "journal_replay_command_order_session_units_geometry_generation_database_owner_identity": {
            "journal_generation": "journal-replay-v44-731",
            "command_order_generation": "journal-replay-v44-731",
            "session_units_generation": "journal-replay-v44-731",
            "geometry_generation": "journal-replay-v44-731",
            "status_generation": "journal-replay-v44-731",
            "database_generation": "journal-replay-v44-731",
            "mesh_generation": "journal-replay-v44-731",
            "result_generation": "journal-replay-v44-731",
            "journal_commands": ["reset", "set units mm", "create brick 1", "mesh volume 1"],
            "replay_journal_commands": ["reset", "set units mm", "create brick 1", "mesh volume 1"],
            "session_units": "mm",
            "replay_session_units": "mm",
            "command_status": ["success", "success", "success", "success"],
            "replay_command_status": ["success", "success", "success", "success"],
            "geometry_generation_id": 731,
            "replay_geometry_generation_id": 731,
            "database_owner": "headless:journal-731",
            "replay_database_owner": "headless:journal-731",
            "mesh_export_sha256": "b" * 64,
            "replay_mesh_export_sha256": "b" * 64,
            "result_sha256": "c" * 64,
            "replay_result_sha256": "c" * 64,
            "accepted_result_sha256": "c" * 64,
        },
        "mesh_quality_metric_reference_element_dimension_block_export_owner_identity": {
            "quality_generation": "quality-report-v44-731",
            "reference_element_generation": "quality-report-v44-731",
            "dimension_generation": "quality-report-v44-731",
            "metric_generation": "quality-report-v44-731",
            "block_generation": "quality-report-v44-731",
            "export_generation": "quality-report-v44-731",
            "owner_generation": "quality-report-v44-731",
            "result_generation": "quality-report-v44-731",
            "reference_element": "hex8",
            "replay_reference_element": "hex8",
            "dimension": 3,
            "replay_dimension": 3,
            "metric_definition": "scaled_jacobian",
            "replay_metric_definition": "scaled_jacobian",
            "minimum_scaled_jacobian": 0.36,
            "replay_minimum_scaled_jacobian": 0.36,
            "block_membership": {"block:periodic_hex": [1]},
            "replay_block_membership": {"block:periodic_hex": [1]},
            "export_generation_id": 731,
            "replay_export_generation_id": 731,
            "database_owner": "headless:quality-731",
            "replay_database_owner": "headless:quality-731",
            "quality_export_sha256": "d" * 64,
            "accepted_quality_export_sha256": "d" * 64,
        },
    }


def test_v44_positive_public_and_source_contracts() -> None:
    row = _summary()
    assert cubit_conformal_hex_pyramid_tet_interface_gate(row)["status"] == "ok"
    assert cubit_mixed_transition_source_gate(row)["status"] == "ok"


def test_v44_rejects_periodic_hex_interface_mismatch() -> None:
    row = _summary()
    identity = row["hex_sweep_periodic_interface_quality_jacobian_block_sideset_export_generation_identity"]
    identity["result_paired_target_node_ids"] = [201, 202, 203, 999]
    identity["result_interface_normal_dot"] = 1.0
    assert cubit_conformal_hex_pyramid_tet_interface_gate(row)["status"] == "needs_attention"


def test_v44_rejects_mixed_transition_quality_mismatch() -> None:
    row = _summary()
    identity = row["hex_sweep_periodic_interface_quality_jacobian_block_sideset_export_generation_identity"]
    identity["result_minimum_scaled_jacobian"] = -0.2
    identity["result_mesh_owner"] = "gui:old"
    assert cubit_conformal_hex_pyramid_tet_interface_gate(row)["status"] == "needs_attention"


def test_v44_rejects_journal_replay_mismatch() -> None:
    row = _summary()
    identity = row["journal_replay_command_order_session_units_geometry_generation_database_owner_identity"]
    identity["replay_journal_commands"] = ["reset", "create brick 1", "set units mm", "mesh volume 1"]
    identity["replay_session_units"] = "in"
    assert cubit_mixed_transition_source_gate(row)["status"] == "needs_attention"


def test_v44_rejects_quality_report_mismatch() -> None:
    row = _summary()
    identity = row["mesh_quality_metric_reference_element_dimension_block_export_owner_identity"]
    identity["replay_reference_element"] = "tet4"
    identity["replay_dimension"] = 2
    assert cubit_mixed_transition_source_gate(row)["status"] == "needs_attention"


def test_v44_case_ids_are_promoted_to_owner_regression() -> None:
    assert all(case_id.startswith("v44_") for case_id in (_PUBLIC_CASE, _MIXED_CASE, _JOURNAL_CASE, _QUALITY_CASE))


def test_v44_combined_source_gate_does_not_mask_journal_lineage_failure() -> None:
    row = _summary()
    row[
        "journal_replay_command_order_session_units_geometry_generation_database_owner_identity"
    ]["command_order_generation"] = "journal-replay-v44-stale"
    result = cubit_mixed_transition_source_gate(row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["journal.generation_lineage"] is False
    assert result["checks"]["quality.generation_lineage"] is True


def test_v44_rejects_missing_result_digest() -> None:
    row = _summary()
    identity = row[
        "journal_replay_command_order_session_units_geometry_generation_database_owner_identity"
    ]
    identity.pop("result_sha256")
    assert cubit_mixed_transition_source_gate(row)["status"] == "needs_attention"


def test_v44_rejects_malformed_numeric_values_without_raising() -> None:
    row = _summary()
    identity = row[
        "hex_sweep_periodic_interface_quality_jacobian_block_sideset_export_generation_identity"
    ]
    identity["minimum_scaled_jacobian"] = {"bad": "value"}
    assert cubit_conformal_hex_pyramid_tet_interface_gate(row)["status"] == "needs_attention"

    row = _summary()
    identity = row[
        "mesh_quality_metric_reference_element_dimension_block_export_owner_identity"
    ]
    identity["export_generation_id"] = {"bad": "value"}
    assert cubit_mixed_transition_source_gate(row)["status"] == "needs_attention"
