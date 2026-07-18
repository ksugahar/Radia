from __future__ import annotations

import json

from test_coreform_generalization_v32 import (
    _with_v32_boundary_transition_journal_exodus_identity,
    cubit_conformal_hex_pyramid_tet_interface_gate,
    cubit_mixed_transition_source_gate,
    summary,
)


_PROMOTED_CASE_IDS = (
    "v33_public_hex_sheet_pillow_cell_face_incidence_euler_boundary_shell_jacobian_mismatch",
    "v33_public_multiblock_conformal_interface_merge_tolerance_face_owner_duplicate_cell_mismatch",
    "v33_source_command_failure_atomic_rollback_error_index_entity_allocator_session_mismatch",
    "v33_source_cub_save_open_roundtrip_kernel_entity_name_attribute_group_mesh_checksum_mismatch",
)


def _with_v33_topology_transaction_roundtrip_identity(row: dict) -> dict:
    row = _with_v32_boundary_transition_journal_exodus_identity(row)
    generation = "hex-sheet-pillow-201"
    row[
        "hex_sheet_pillow_incidence_euler_shell_orientation_block_jacobian_mesh_result_generation_identity"
    ] = {
        "sheet_generation": generation,
        "incidence_generation": generation,
        "euler_generation": generation,
        "shell_generation": generation,
        "orientation_generation": generation,
        "block_generation": generation,
        "jacobian_generation": generation,
        "mesh_generation": generation,
        "result_generation": generation,
        "operation": "pillow",
        "result_operation": "pillow",
        "cell_ids": [1, 2, 3, 4],
        "result_cell_ids": [1, 2, 3, 4],
        "cell_face_incidence": [
            [1, 101],
            [1, 102],
            [2, 102],
            [2, 103],
            [3, 103],
            [3, 104],
            [4, 104],
        ],
        "result_cell_face_incidence": [
            [1, 101],
            [1, 102],
            [2, 102],
            [2, 103],
            [3, 103],
            [3, 104],
            [4, 104],
        ],
        "vertex_count": 18,
        "result_vertex_count": 18,
        "edge_count": 33,
        "result_edge_count": 33,
        "face_count": 20,
        "result_face_count": 20,
        "cell_count": 4,
        "result_cell_count": 4,
        "euler_characteristic": 1,
        "result_euler_characteristic": 1,
        "boundary_shell_count": 1,
        "result_boundary_shell_count": 1,
        "cell_orientation_signs": [1, 1, 1, 1],
        "result_cell_orientation_signs": [1, 1, 1, 1],
        "block_id": 101,
        "result_block_id": 101,
        "scaled_jacobians": [0.42, 0.38, 0.51, 0.46],
        "result_scaled_jacobians": [0.42, 0.38, 0.51, 0.46],
        "sheet_mesh_sha256": "1" * 64,
        "accepted_sheet_mesh_sha256": "1" * 64,
    }
    generation = "multiblock-interface-201"
    row[
        "multiblock_interface_merge_face_owner_block_sideset_duplicate_jacobian_export_generation_identity"
    ] = {
        "interface_generation": generation,
        "merge_generation": generation,
        "node_generation": generation,
        "face_generation": generation,
        "owner_generation": generation,
        "block_generation": generation,
        "sideset_generation": generation,
        "duplicate_generation": generation,
        "jacobian_generation": generation,
        "export_generation": generation,
        "result_generation": generation,
        "merge_tolerance_m": 1.0e-8,
        "result_merge_tolerance_m": 1.0e-8,
        "interface_node_pairs": [[1, 101], [2, 102], [3, 103], [4, 104]],
        "result_interface_node_pairs": [[1, 101], [2, 102], [3, 103], [4, 104]],
        "coincident_face_pairs": [[10, 20]],
        "result_coincident_face_pairs": [[10, 20]],
        "face_owners": [[10, 101], [20, 102]],
        "result_face_owners": [[10, 101], [20, 102]],
        "block_ids": [101, 102],
        "result_block_ids": [101, 102],
        "sideset_ids": [201, 202],
        "result_sideset_ids": [201, 202],
        "duplicate_cell_count": 0,
        "result_duplicate_cell_count": 0,
        "minimum_scaled_jacobian": 0.36,
        "result_minimum_scaled_jacobian": 0.36,
        "multiblock_export_sha256": "2" * 64,
        "accepted_multiblock_export_sha256": "2" * 64,
    }
    generation = "command-failure-transaction-201"
    row[
        "command_failure_atomic_rollback_error_entity_allocator_undo_session_result_generation_identity"
    ] = {
        "failure_generation": generation,
        "command_generation": generation,
        "error_generation": generation,
        "rollback_generation": generation,
        "allocator_generation": generation,
        "undo_generation": generation,
        "session_generation": generation,
        "result_generation": generation,
        "failed_command_index": 7,
        "reported_failed_command_index": 7,
        "error_code": 1201,
        "reported_error_code": 1201,
        "failure_category": "invalid_entity_reference",
        "reported_failure_category": "invalid_entity_reference",
        "pre_transaction_model_sha256": "3" * 64,
        "rolled_back_model_sha256": "3" * 64,
        "next_entity_id": 41,
        "rolled_back_next_entity_id": 41,
        "undo_depth": 2,
        "rolled_back_undo_depth": 2,
        "session_owner": "headless-batch-session-201",
        "rolled_back_session_owner": "headless-batch-session-201",
        "failure_result_sha256": "4" * 64,
        "accepted_failure_result_sha256": "4" * 64,
    }
    generation = "cub-roundtrip-201"
    row[
        "cub_roundtrip_kernel_entity_name_attribute_group_mesh_model_file_result_generation_identity"
    ] = {
        "roundtrip_generation": generation,
        "kernel_generation": generation,
        "entity_generation": generation,
        "attribute_generation": generation,
        "group_generation": generation,
        "mesh_generation": generation,
        "model_generation": generation,
        "file_generation": generation,
        "result_generation": generation,
        "kernel_version": "2026.6",
        "reopened_kernel_version": "2026.6",
        "entity_names": [["volume", 1, "core"], ["surface", 7, "wall"]],
        "reopened_entity_names": [["volume", 1, "core"], ["surface", 7, "wall"]],
        "entity_attributes": [["volume", 1, "material", "steel"]],
        "reopened_entity_attributes": [["volume", 1, "material", "steel"]],
        "group_memberships": [["moving", "volume", 1], ["walls", "surface", 7]],
        "reopened_group_memberships": [
            ["moving", "volume", 1],
            ["walls", "surface", 7],
        ],
        "mesh_state": {
            "hex": 64,
            "quad": 96,
            "block_ids": [101],
            "sideset_ids": [201],
        },
        "reopened_mesh_state": {
            "hex": 64,
            "quad": 96,
            "block_ids": [101],
            "sideset_ids": [201],
        },
        "model_generation": "model-201",
        "reopened_model_generation": "model-201",
        "cub_file_sha256": "5" * 64,
        "reopened_cub_file_sha256": "5" * 64,
        "roundtrip_result_sha256": "6" * 64,
        "accepted_roundtrip_result_sha256": "6" * 64,
    }
    return row


def test_v33_positive_topology_transaction_and_roundtrip() -> None:
    row = _with_v33_topology_transaction_roundtrip_identity(summary())
    assert json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))["status"] == "ok"
    assert json.loads(cubit_mixed_transition_source_gate(row))["status"] == "ok"


def test_v33_public_hex_sheet_pillow_cell_face_incidence_euler_boundary_shell_jacobian_mismatch() -> None:
    row = _with_v33_topology_transaction_roundtrip_identity(summary())
    row[
        "hex_sheet_pillow_incidence_euler_shell_orientation_block_jacobian_mesh_result_generation_identity"
    ].update(
        {
            "incidence_generation": "hex-sheet-pillow-200",
            "shell_generation": "hex-sheet-pillow-199",
            "result_generation": "hex-sheet-pillow-198",
            "result_operation": "sheet_extract",
            "result_cell_ids": [1, 2, 4],
            "result_cell_face_incidence": [[1, 101], [2, 999], [4, 104]],
            "result_edge_count": 31,
            "result_face_count": 18,
            "result_cell_count": 3,
            "result_euler_characteristic": 2,
            "result_boundary_shell_count": 2,
            "result_cell_orientation_signs": [1, -1, 1],
            "result_block_id": 102,
            "result_scaled_jacobians": [0.42, -0.05, 0.46],
            "accepted_sheet_mesh_sha256": "7" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hex_sheet_pillows_use_current_incidence_euler_shell_orientation_block_jacobian_and_mesh"
    ]


def test_v33_public_multiblock_conformal_interface_merge_tolerance_face_owner_duplicate_cell_mismatch() -> None:
    row = _with_v33_topology_transaction_roundtrip_identity(summary())
    row[
        "multiblock_interface_merge_face_owner_block_sideset_duplicate_jacobian_export_generation_identity"
    ].update(
        {
            "merge_generation": "multiblock-interface-200",
            "owner_generation": "multiblock-interface-199",
            "result_generation": "multiblock-interface-198",
            "result_merge_tolerance_m": 1.0e-3,
            "result_interface_node_pairs": [[1, 104], [2, 103]],
            "result_coincident_face_pairs": [[10, 21]],
            "result_face_owners": [[10, 102], [20, 101]],
            "result_block_ids": [101, 101],
            "result_sideset_ids": [201, 201],
            "result_duplicate_cell_count": 2,
            "result_minimum_scaled_jacobian": -0.02,
            "accepted_multiblock_export_sha256": "8" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "multiblock_interfaces_use_current_merge_nodes_faces_owners_sets_duplicates_jacobian_and_export"
    ]


def test_v33_source_command_failure_atomic_rollback_error_index_entity_allocator_session_mismatch() -> None:
    row = _with_v33_topology_transaction_roundtrip_identity(summary())
    row[
        "command_failure_atomic_rollback_error_entity_allocator_undo_session_result_generation_identity"
    ].update(
        {
            "command_generation": "command-failure-transaction-200",
            "rollback_generation": "command-failure-transaction-199",
            "result_generation": "command-failure-transaction-198",
            "reported_failed_command_index": 6,
            "reported_error_code": 0,
            "reported_failure_category": "success",
            "rolled_back_model_sha256": "9" * 64,
            "rolled_back_next_entity_id": 44,
            "rolled_back_undo_depth": 3,
            "rolled_back_session_owner": "gui-session-old",
            "accepted_failure_result_sha256": "a" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "failed_commands_rollback_model_allocator_undo_session_and_result_atomically"
    ]


def test_v33_source_cub_save_open_roundtrip_kernel_entity_name_attribute_group_mesh_checksum_mismatch() -> None:
    row = _with_v33_topology_transaction_roundtrip_identity(summary())
    row[
        "cub_roundtrip_kernel_entity_name_attribute_group_mesh_model_file_result_generation_identity"
    ].update(
        {
            "kernel_generation": "cub-roundtrip-200",
            "mesh_generation": "cub-roundtrip-199",
            "result_generation": "cub-roundtrip-198",
            "reopened_kernel_version": "2025.8",
            "reopened_entity_names": [["volume", 1, "old_core"]],
            "reopened_entity_attributes": [["volume", 1, "material", "air"]],
            "reopened_group_memberships": [["walls", "surface", 8]],
            "reopened_mesh_state": {
                "hex": 32,
                "quad": 48,
                "block_ids": [102],
                "sideset_ids": [],
            },
            "reopened_model_generation": "model-200",
            "reopened_cub_file_sha256": "b" * 64,
            "accepted_roundtrip_result_sha256": "c" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cub_roundtrips_preserve_kernel_entities_attributes_groups_mesh_model_file_and_result"
    ]
