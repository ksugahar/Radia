from __future__ import annotations

import json

from test_coreform_generalization_v31 import (
    _with_v31_high_order_periodic_lineage_exodus_identity,
    cubit_conformal_hex_pyramid_tet_interface_gate,
    cubit_mixed_transition_source_gate,
    summary,
)


_PROMOTED_CASE_IDS = (
    "v32_public_hex_boundary_layer_thickness_growth_collision_jacobian_block_sideset_mismatch",
    "v32_public_pyramid_transition_tet_hex_interface_diagonal_orientation_region_export_mismatch",
    "v32_source_journal_undo_transaction_idempotency_entity_allocation_save_restore_mismatch",
    "v32_source_exodus_assembly_coordinate_frame_qa_record_time_metadata_checksum_mismatch",
)


def _with_v32_boundary_transition_journal_exodus_identity(row: dict) -> dict:
    row = _with_v31_high_order_periodic_lineage_exodus_identity(row)
    generation = "hex-boundary-layer-closure-191"
    row[
        "hex_boundary_layer_thickness_growth_count_collision_jacobian_block_sideset_mesh_generation_identity"
    ] = {
        "boundary_layer_generation": generation,
        "thickness_generation": generation,
        "growth_generation": generation,
        "collision_generation": generation,
        "jacobian_generation": generation,
        "block_generation": generation,
        "sideset_generation": generation,
        "mesh_generation": generation,
        "result_generation": generation,
        "first_layer_thickness_m": 1.0e-4,
        "result_first_layer_thickness_m": 1.0e-4,
        "growth_law": "geometric",
        "result_growth_law": "geometric",
        "growth_ratio": 1.25,
        "result_growth_ratio": 1.25,
        "layer_count": 4,
        "result_layer_count": 4,
        "layer_thicknesses_m": [1.0e-4, 1.25e-4, 1.5625e-4, 1.953125e-4],
        "result_layer_thicknesses_m": [1.0e-4, 1.25e-4, 1.5625e-4, 1.953125e-4],
        "collision_handling": "truncate_and_rebalance",
        "result_collision_handling": "truncate_and_rebalance",
        "minimum_scaled_jacobian": 0.31,
        "result_minimum_scaled_jacobian": 0.31,
        "block_id": 101,
        "result_block_id": 101,
        "wall_sideset_id": 201,
        "result_wall_sideset_id": 201,
        "boundary_layer_mesh_sha256": "1" * 64,
        "accepted_boundary_layer_mesh_sha256": "1" * 64,
    }
    generation = "pyramid-transition-closure-191"
    row[
        "pyramid_transition_diagonal_orientation_shared_face_region_jacobian_export_generation_identity"
    ] = {
        "transition_generation": generation,
        "diagonal_generation": generation,
        "orientation_generation": generation,
        "face_generation": generation,
        "region_generation": generation,
        "jacobian_generation": generation,
        "export_generation": generation,
        "result_generation": generation,
        "interface_diagonal": [2, 4],
        "result_interface_diagonal": [2, 4],
        "pyramid_orientation_signs": [1, 1, 1, 1],
        "result_pyramid_orientation_signs": [1, 1, 1, 1],
        "shared_face_connectivity": [[1, 2, 3, 4], [1, 2, 5], [2, 3, 5]],
        "result_shared_face_connectivity": [[1, 2, 3, 4], [1, 2, 5], [2, 3, 5]],
        "region_owners": ["hex_region", "transition_region", "tet_region"],
        "result_region_owners": ["hex_region", "transition_region", "tet_region"],
        "scaled_jacobians": [0.44, 0.28, 0.35],
        "result_scaled_jacobians": [0.44, 0.28, 0.35],
        "transition_export_sha256": "2" * 64,
        "accepted_transition_export_sha256": "2" * 64,
    }
    generation = "journal-transaction-closure-191"
    row[
        "journal_undo_transaction_idempotency_entity_allocation_save_restore_session_result_generation_identity"
    ] = {
        "transaction_generation": generation,
        "undo_generation": generation,
        "allocation_generation": generation,
        "idempotency_generation": generation,
        "save_generation": generation,
        "restore_generation": generation,
        "session_generation": generation,
        "result_generation": generation,
        "transaction_id": "txn-191",
        "replayed_transaction_id": "txn-191",
        "undo_count": 1,
        "replayed_undo_count": 1,
        "redo_count": 1,
        "replayed_redo_count": 1,
        "allocated_entity_ids": [301, 302, 303],
        "replayed_allocated_entity_ids": [301, 302, 303],
        "idempotency_sha256": "3" * 64,
        "replayed_idempotency_sha256": "3" * 64,
        "save_point_id": "save-191",
        "restored_save_point_id": "save-191",
        "saved_model_sha256": "4" * 64,
        "restored_model_sha256": "4" * 64,
        "session_owner": "headless-batch-session-191",
        "restored_session_owner": "headless-batch-session-191",
        "journal_result_sha256": "5" * 64,
        "accepted_journal_result_sha256": "5" * 64,
    }
    generation = "exodus-assembly-closure-191"
    row[
        "exodus_assembly_membership_frame_qa_time_block_sideset_file_result_generation_identity"
    ] = {
        "exodus_generation": generation,
        "assembly_generation": generation,
        "frame_generation": generation,
        "qa_generation": generation,
        "time_generation": generation,
        "block_generation": generation,
        "sideset_generation": generation,
        "file_generation": generation,
        "result_generation": generation,
        "assembly_ids": [1001, 1002],
        "decoded_assembly_ids": [1001, 1002],
        "assembly_members": [[101, 102], [201]],
        "decoded_assembly_members": [[101, 102], [201]],
        "coordinate_frame": "global_cartesian_m",
        "decoded_coordinate_frame": "global_cartesian_m",
        "qa_record": ["CAE-AI Lab", "2026.6", "2026-07-18", "14:00:00"],
        "decoded_qa_record": ["CAE-AI Lab", "2026.6", "2026-07-18", "14:00:00"],
        "time_values_s": [0.0, 0.1, 0.2],
        "decoded_time_values_s": [0.0, 0.1, 0.2],
        "block_owners": [[101, 1001], [102, 1001], [201, 1002]],
        "decoded_block_owners": [[101, 1001], [102, 1001], [201, 1002]],
        "sideset_owners": [[301, 1001], [302, 1002]],
        "decoded_sideset_owners": [[301, 1001], [302, 1002]],
        "exodus_file_sha256": "6" * 64,
        "decoded_exodus_file_sha256": "6" * 64,
        "assembly_table_sha256": "7" * 64,
        "decoded_assembly_table_sha256": "7" * 64,
    }
    return row


def test_v32_positive_boundary_transition_journal_and_exodus() -> None:
    row = _with_v32_boundary_transition_journal_exodus_identity(summary())
    assert json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))["status"] == "ok"
    assert json.loads(cubit_mixed_transition_source_gate(row))["status"] == "ok"


def test_v32_public_hex_boundary_layer_thickness_growth_collision_jacobian_block_sideset_mismatch() -> None:
    row = _with_v32_boundary_transition_journal_exodus_identity(summary())
    row[
        "hex_boundary_layer_thickness_growth_count_collision_jacobian_block_sideset_mesh_generation_identity"
    ].update(
        {
            "thickness_generation": "hex-boundary-layer-190",
            "collision_generation": "hex-boundary-layer-189",
            "result_generation": "hex-boundary-layer-188",
            "result_first_layer_thickness_m": 2.0e-4,
            "result_growth_law": "linear",
            "result_growth_ratio": 1.0,
            "result_layer_count": 3,
            "result_layer_thicknesses_m": [2.0e-4, 2.0e-4, 2.0e-4],
            "result_collision_handling": "overlap",
            "result_minimum_scaled_jacobian": -0.03,
            "result_block_id": 102,
            "result_wall_sideset_id": 202,
            "accepted_boundary_layer_mesh_sha256": "8" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hex_boundary_layers_use_current_thickness_growth_count_collision_jacobian_block_sideset_and_mesh"
    ]


def test_v32_public_pyramid_transition_tet_hex_interface_diagonal_orientation_region_export_mismatch() -> None:
    row = _with_v32_boundary_transition_journal_exodus_identity(summary())
    row[
        "pyramid_transition_diagonal_orientation_shared_face_region_jacobian_export_generation_identity"
    ].update(
        {
            "diagonal_generation": "pyramid-transition-190",
            "orientation_generation": "pyramid-transition-189",
            "result_generation": "pyramid-transition-188",
            "result_interface_diagonal": [1, 3],
            "result_pyramid_orientation_signs": [1, -1, 1, -1],
            "result_shared_face_connectivity": [[1, 4, 3, 2], [1, 5, 2]],
            "result_region_owners": ["tet_region", "transition_region", "hex_region"],
            "result_scaled_jacobians": [0.44, -0.08, 0.35],
            "accepted_transition_export_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pyramid_transitions_use_current_diagonal_orientation_faces_regions_jacobians_and_export"
    ]


def test_v32_source_journal_undo_transaction_idempotency_entity_allocation_save_restore_mismatch() -> None:
    row = _with_v32_boundary_transition_journal_exodus_identity(summary())
    row[
        "journal_undo_transaction_idempotency_entity_allocation_save_restore_session_result_generation_identity"
    ].update(
        {
            "undo_generation": "journal-transaction-190",
            "restore_generation": "journal-transaction-189",
            "result_generation": "journal-transaction-188",
            "replayed_transaction_id": "txn-190",
            "replayed_undo_count": 2,
            "replayed_redo_count": 0,
            "replayed_allocated_entity_ids": [301, 303, 304],
            "replayed_idempotency_sha256": "a" * 64,
            "restored_save_point_id": "save-190",
            "restored_model_sha256": "b" * 64,
            "restored_session_owner": "gui-session-old",
            "accepted_journal_result_sha256": "c" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "journal_replays_use_current_transaction_undo_allocation_idempotency_save_restore_session_and_result"
    ]


def test_v32_source_exodus_assembly_coordinate_frame_qa_record_time_metadata_checksum_mismatch() -> None:
    row = _with_v32_boundary_transition_journal_exodus_identity(summary())
    row[
        "exodus_assembly_membership_frame_qa_time_block_sideset_file_result_generation_identity"
    ].update(
        {
            "assembly_generation": "exodus-assembly-190",
            "time_generation": "exodus-assembly-189",
            "result_generation": "exodus-assembly-188",
            "decoded_assembly_ids": [1002, 1001],
            "decoded_assembly_members": [[201], [101]],
            "decoded_coordinate_frame": "part_local_mm",
            "decoded_qa_record": ["unknown", "old", "", ""],
            "decoded_time_values_s": [0.0, 0.2, 0.1],
            "decoded_block_owners": [[101, 1002], [201, 1001]],
            "decoded_sideset_owners": [[301, 1002]],
            "decoded_exodus_file_sha256": "d" * 64,
            "decoded_assembly_table_sha256": "e" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "exodus_assemblies_use_current_membership_frame_qa_time_block_sideset_file_and_result"
    ]
