from __future__ import annotations

import json

from test_coreform_generalization_v29 import (
    _with_v29_map_morph_exodus_cad_identity,
    cubit_conformal_hex_pyramid_tet_interface_gate,
    cubit_mixed_transition_source_gate,
    summary,
)


_PROMOTED_CASE_IDS = (
    "v30_public_sweep_hex_twist_source_target_interval_edge_map_scaled_jacobian_mismatch",
    "v30_public_hex_tet_pyramid_transition_face_conformity_orientation_quality_export_mismatch",
    "v30_source_journal_undo_transaction_command_order_entity_id_session_digest_mismatch",
    "v30_source_exodus_sideset_distribution_factor_topology_orientation_block_time_digest_mismatch",
)


def _with_v30_sweep_transition_journal_sideset_identity(row: dict) -> dict:
    row = _with_v29_map_morph_exodus_cad_identity(row)
    generation = "sweep-hex-twist-171"
    row[
        "sweep_hex_twist_source_target_interval_edge_map_orientation_jacobian_block_mesh_generation_identity"
    ] = {
        "sweep_generation": generation,
        "twist_sweep_generation": generation,
        "pairing_sweep_generation": generation,
        "interval_sweep_generation": generation,
        "edge_map_sweep_generation": generation,
        "orientation_sweep_generation": generation,
        "quality_sweep_generation": generation,
        "block_sweep_generation": generation,
        "result_sweep_generation": generation,
        "twist_angle_deg": 90.0,
        "result_twist_angle_deg": 90.0,
        "source_face_id": 101,
        "result_source_face_id": 101,
        "target_face_id": 201,
        "result_target_face_id": 201,
        "axial_interval_count": 12,
        "result_axial_interval_count": 12,
        "source_edge_ids": [1, 2, 3, 4],
        "result_source_edge_ids": [1, 2, 3, 4],
        "target_edge_ids": [5, 6, 7, 8],
        "result_target_edge_ids": [5, 6, 7, 8],
        "edge_map": [[1, 5], [2, 6], [3, 7], [4, 8]],
        "result_edge_map": [[1, 5], [2, 6], [3, 7], [4, 8]],
        "sweep_orientation": "right_handed_source_to_target",
        "result_sweep_orientation": "right_handed_source_to_target",
        "scaled_jacobians": [0.62, 0.55, 0.48],
        "result_scaled_jacobians": [0.62, 0.55, 0.48],
        "block_id": 30,
        "result_block_id": 30,
        "sweep_mesh_sha256": "1" * 64,
        "result_sweep_mesh_sha256": "1" * 64,
    }
    generation = "mixed-transition-171"
    row[
        "hex_tet_pyramid_transition_face_nodes_orientation_conformity_quality_block_sideset_export_generation_identity"
    ] = {
        "transition_generation": generation,
        "face_transition_generation": generation,
        "node_transition_generation": generation,
        "orientation_transition_generation": generation,
        "conformity_transition_generation": generation,
        "quality_transition_generation": generation,
        "block_transition_generation": generation,
        "sideset_transition_generation": generation,
        "export_transition_generation": generation,
        "result_transition_generation": generation,
        "hex_element_ids": [10],
        "result_hex_element_ids": [10],
        "pyramid_element_ids": [301, 302],
        "result_pyramid_element_ids": [301, 302],
        "tet_element_ids": [401, 402],
        "result_tet_element_ids": [401, 402],
        "transition_face_node_ids": [[1, 2, 3, 4], [2, 5, 6], [3, 6, 7]],
        "result_transition_face_node_ids": [[1, 2, 3, 4], [2, 5, 6], [3, 6, 7]],
        "face_owner_pairs": [["hex:10", "pyramid:301"], ["pyramid:301", "tet:401"], ["pyramid:302", "tet:402"]],
        "result_face_owner_pairs": [["hex:10", "pyramid:301"], ["pyramid:301", "tet:401"], ["pyramid:302", "tet:402"]],
        "opposed_face_orientation_signs": [[1, -1], [1, -1], [1, -1]],
        "result_opposed_face_orientation_signs": [[1, -1], [1, -1], [1, -1]],
        "unmatched_transition_face_count": 0,
        "result_unmatched_transition_face_count": 0,
        "minimum_scaled_jacobian": 0.31,
        "result_minimum_scaled_jacobian": 0.31,
        "block_ids": [30, 40, 50],
        "result_block_ids": [30, 40, 50],
        "sideset_ids": [100, 200],
        "result_sideset_ids": [100, 200],
        "transition_mesh_sha256": "2" * 64,
        "result_transition_mesh_sha256": "2" * 64,
        "transition_export_sha256": "3" * 64,
        "accepted_transition_export_sha256": "3" * 64,
    }
    generation = "journal-transaction-171"
    row[
        "journal_undo_transaction_command_order_entity_session_model_digest_generation_identity"
    ] = {
        "journal_generation": generation,
        "undo_journal_generation": generation,
        "command_journal_generation": generation,
        "entity_journal_generation": generation,
        "session_journal_generation": generation,
        "model_journal_generation": generation,
        "digest_journal_generation": generation,
        "result_journal_generation": generation,
        "undo_transaction_id": "undo-tx-171",
        "replayed_undo_transaction_id": "undo-tx-171",
        "command_sequence": ["reset", "brick x 1 y 1 z 1", "mesh volume 1"],
        "replayed_command_sequence": ["reset", "brick x 1 y 1 z 1", "mesh volume 1"],
        "command_ordinals": [1, 2, 3],
        "replayed_command_ordinals": [1, 2, 3],
        "entity_ids": [1],
        "replayed_entity_ids": [1],
        "active_session_id": "headless-session-171",
        "replayed_active_session_id": "headless-session-171",
        "model_generation_id": "model-generation-171",
        "replayed_model_generation_id": "model-generation-171",
        "journal_sha256": "4" * 64,
        "loaded_journal_sha256": "4" * 64,
        "journal_result_sha256": "5" * 64,
        "accepted_journal_result_sha256": "5" * 64,
    }
    generation = "exodus-sideset-171"
    row[
        "exodus_sideset_distribution_factor_topology_orientation_block_time_file_generation_identity"
    ] = {
        "sideset_generation": generation,
        "topology_sideset_generation": generation,
        "orientation_sideset_generation": generation,
        "factor_sideset_generation": generation,
        "block_sideset_generation": generation,
        "time_sideset_generation": generation,
        "file_sideset_generation": generation,
        "result_sideset_generation": generation,
        "sideset_ids": [100, 200],
        "decoded_sideset_ids": [100, 200],
        "face_topologies": ["quad4", "tri3"],
        "decoded_face_topologies": ["quad4", "tri3"],
        "face_orientation_signs": [1, -1],
        "decoded_face_orientation_signs": [1, -1],
        "distribution_factors": [[1.0, 1.0, 1.0, 1.0], [0.5, 0.5, 0.5]],
        "decoded_distribution_factors": [[1.0, 1.0, 1.0, 1.0], [0.5, 0.5, 0.5]],
        "block_owner_ids": [30, 40],
        "decoded_block_owner_ids": [30, 40],
        "time_step_index": 2,
        "decoded_time_step_index": 2,
        "exodus_file_sha256": "6" * 64,
        "decoded_exodus_file_sha256": "6" * 64,
        "sideset_table_sha256": "7" * 64,
        "decoded_sideset_table_sha256": "7" * 64,
    }
    return row


def test_v30_positive_sweep_transition_journal_and_sideset() -> None:
    row = _with_v30_sweep_transition_journal_sideset_identity(summary())
    assert json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))["status"] == "ok"
    assert json.loads(cubit_mixed_transition_source_gate(row))["status"] == "ok"


def test_v30_public_sweep_hex_twist_source_target_interval_edge_map_scaled_jacobian_mismatch() -> None:
    row = _with_v30_sweep_transition_journal_sideset_identity(summary())
    row[
        "sweep_hex_twist_source_target_interval_edge_map_orientation_jacobian_block_mesh_generation_identity"
    ].update(
        {
            "twist_sweep_generation": "sweep-hex-twist-170",
            "edge_map_sweep_generation": "sweep-hex-twist-169",
            "result_sweep_generation": "sweep-hex-twist-168",
            "result_twist_angle_deg": -90.0,
            "result_source_face_id": 201,
            "result_target_face_id": 101,
            "result_axial_interval_count": 11,
            "result_source_edge_ids": [1, 2, 3],
            "result_target_edge_ids": [8, 7, 6, 5],
            "result_edge_map": [[1, 8], [2, 7]],
            "result_sweep_orientation": "left_handed_target_to_source",
            "result_scaled_jacobians": [0.62, -0.2],
            "result_block_id": 40,
            "result_sweep_mesh_sha256": "8" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "swept_hexes_use_current_twist_pairing_intervals_edge_map_orientation_jacobians_block_and_mesh"
    ]


def test_v30_public_hex_tet_pyramid_transition_face_conformity_orientation_quality_export_mismatch() -> None:
    row = _with_v30_sweep_transition_journal_sideset_identity(summary())
    row[
        "hex_tet_pyramid_transition_face_nodes_orientation_conformity_quality_block_sideset_export_generation_identity"
    ].update(
        {
            "face_transition_generation": "mixed-transition-170",
            "conformity_transition_generation": "mixed-transition-169",
            "result_transition_generation": "mixed-transition-168",
            "result_hex_element_ids": [11],
            "result_pyramid_element_ids": [302, 301],
            "result_tet_element_ids": [402],
            "result_transition_face_node_ids": [[1, 2, 3], [2, 6, 5]],
            "result_face_owner_pairs": [["hex:10", "tet:401"]],
            "result_opposed_face_orientation_signs": [[1, 1]],
            "result_unmatched_transition_face_count": 2,
            "result_minimum_scaled_jacobian": -0.1,
            "result_block_ids": [30, 50],
            "result_sideset_ids": [100, 300],
            "result_transition_mesh_sha256": "9" * 64,
            "accepted_transition_export_sha256": "a" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mixed_transitions_use_current_elements_faces_nodes_orientation_conformity_quality_blocks_sidesets_and_export"
    ]


def test_v30_source_journal_undo_transaction_command_order_entity_id_session_digest_mismatch() -> None:
    row = _with_v30_sweep_transition_journal_sideset_identity(summary())
    row[
        "journal_undo_transaction_command_order_entity_session_model_digest_generation_identity"
    ].update(
        {
            "undo_journal_generation": "journal-transaction-170",
            "session_journal_generation": "journal-transaction-169",
            "result_journal_generation": "journal-transaction-168",
            "replayed_undo_transaction_id": "undo-tx-old",
            "replayed_command_sequence": ["reset", "mesh volume 1", "brick x 1 y 1 z 1"],
            "replayed_command_ordinals": [1, 3, 2],
            "replayed_entity_ids": [2],
            "replayed_active_session_id": "gui-session-old",
            "replayed_model_generation_id": "model-generation-170",
            "loaded_journal_sha256": "b" * 64,
            "accepted_journal_result_sha256": "c" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "journals_use_current_undo_transaction_command_order_entities_headless_session_model_and_digests"
    ]


def test_v30_source_exodus_sideset_distribution_factor_topology_orientation_block_time_digest_mismatch() -> None:
    row = _with_v30_sweep_transition_journal_sideset_identity(summary())
    row[
        "exodus_sideset_distribution_factor_topology_orientation_block_time_file_generation_identity"
    ].update(
        {
            "topology_sideset_generation": "exodus-sideset-170",
            "factor_sideset_generation": "exodus-sideset-169",
            "result_sideset_generation": "exodus-sideset-168",
            "decoded_sideset_ids": [200, 100],
            "decoded_face_topologies": ["tri3", "quad8"],
            "decoded_face_orientation_signs": [-1, -1],
            "decoded_distribution_factors": [[0.5, 0.5], [1.0]],
            "decoded_block_owner_ids": [40, 30],
            "decoded_time_step_index": 1,
            "decoded_exodus_file_sha256": "d" * 64,
            "decoded_sideset_table_sha256": "e" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "exodus_sidesets_use_current_ids_topology_orientation_distribution_blocks_time_and_digests"
    ]
