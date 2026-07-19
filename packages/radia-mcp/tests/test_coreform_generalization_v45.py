from __future__ import annotations

from radia_mcp.cubit.mixed_transition_gate import cubit_conformal_hex_pyramid_tet_interface_gate, cubit_mixed_transition_source_gate
from test_coreform_generalization_v44 import _summary


_PROMOTED_CASE_IDS = (
    "v45_public_hex_sweep_periodic_pairing_curvature_jacobian_block_export_mismatch",
    "v45_public_mixed_tet_hex_pyramid_transition_orientation_volume_area_quality_owner_mismatch",
    "v45_source_headless_journal_units_command_status_geometry_mesh_database_owner_mismatch",
    "v45_source_quality_metric_reference_element_dimension_threshold_block_export_owner_mismatch",
)


def _with_v45(summary: dict) -> dict:
    summary["hex_sweep_periodic_pairing_curvature_jacobian_block_export_identity"] = {
        "generation": "periodic-hex-v45-812", "pairing_generation": "periodic-hex-v45-812",
        "curvature_generation": "periodic-hex-v45-812", "jacobian_generation": "periodic-hex-v45-812",
        "block_generation": "periodic-hex-v45-812", "export_generation": "periodic-hex-v45-812",
        "source_surface_id": 11, "result_source_surface_id": 11, "target_surface_id": 21, "result_target_surface_id": 21,
        "paired_source_node_ids": [101, 102, 103, 104], "result_paired_source_node_ids": [101, 102, 103, 104],
        "paired_target_node_ids": [201, 202, 203, 204], "result_paired_target_node_ids": [201, 202, 203, 204],
        "curvature_order": 2, "result_curvature_order": 2, "minimum_scaled_jacobian": 0.36, "result_minimum_scaled_jacobian": 0.36,
        "minimum_allowed_scaled_jacobian": 0.20, "block_membership": {"block:periodic_hex": [1]}, "result_block_membership": {"block:periodic_hex": [1]},
        "mesh_owner": "headless:periodic-hex-v45-812", "result_mesh_owner": "headless:periodic-hex-v45-812", "mesh_export_sha256": "1" * 64, "accepted_mesh_export_sha256": "1" * 64,
    }
    summary["headless_journal_units_command_status_geometry_mesh_database_owner_identity"] = {
        "generation": "journal-v45-812", "command_order_generation": "journal-v45-812", "session_units_generation": "journal-v45-812", "status_generation": "journal-v45-812", "geometry_generation": "journal-v45-812", "mesh_generation": "journal-v45-812", "database_generation": "journal-v45-812",
        "commands": ["reset", "set units mm", "create brick 1", "mesh volume 1"], "replay_commands": ["reset", "set units mm", "create brick 1", "mesh volume 1"], "session_units": "mm", "replay_session_units": "mm", "command_status": ["success"] * 4, "replay_command_status": ["success"] * 4, "geometry_generation_id": 812, "replay_geometry_generation_id": 812, "mesh_generation_id": 812, "replay_mesh_generation_id": 812, "database_owner": "headless:journal-v45-812", "replay_database_owner": "headless:journal-v45-812", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    summary["quality_metric_reference_element_dimension_threshold_block_export_owner_identity"] = {
        "generation": "quality-v45-812", "reference_element_generation": "quality-v45-812", "dimension_generation": "quality-v45-812", "threshold_generation": "quality-v45-812", "block_generation": "quality-v45-812", "export_generation": "quality-v45-812", "reference_element": "hex8", "replay_reference_element": "hex8", "dimension": 3, "replay_dimension": 3, "quality_threshold": 0.20, "replay_quality_threshold": 0.20, "block_name": "block:periodic_hex", "replay_block_name": "block:periodic_hex", "export_generation_id": 812, "replay_export_generation_id": 812, "database_owner": "headless:quality-v45-812", "replay_database_owner": "headless:quality-v45-812", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
    }
    return summary


def test_v45_positive_identity_contracts() -> None:
    summary = _with_v45(_summary())
    assert cubit_conformal_hex_pyramid_tet_interface_gate(summary)["status"] == "ok"
    assert cubit_mixed_transition_source_gate(summary)["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 4


def test_v45_rejects_periodic_pairing_or_curvature_mismatch() -> None:
    summary = _with_v45(_summary())
    summary["hex_sweep_periodic_pairing_curvature_jacobian_block_export_identity"]["result_curvature_order"] = 1
    assert cubit_conformal_hex_pyramid_tet_interface_gate(summary)["status"] == "needs_attention"


def test_v45_rejects_journal_replay_mismatch() -> None:
    summary = _with_v45(_summary())
    summary["headless_journal_units_command_status_geometry_mesh_database_owner_identity"]["replay_session_units"] = "in"
    assert cubit_mixed_transition_source_gate(summary)["status"] == "needs_attention"


def test_v45_rejects_quality_threshold_mismatch() -> None:
    summary = _with_v45(_summary())
    summary["quality_metric_reference_element_dimension_threshold_block_export_owner_identity"]["replay_quality_threshold"] = -0.1
    assert cubit_mixed_transition_source_gate(summary)["status"] == "needs_attention"
