from __future__ import annotations

from test_coreform_generalization_v35 import (
    _public_result,
    _source_result,
    _with_v35_metric_curved_sideset_sculpt_identity,
    summary,
)


_PROMOTED_CASE_IDS = (
    "v36_public_hex_sweep_source_target_layers_orientation_jacobian_boundary_owner_mismatch",
    "v36_public_mixed_hex_tet_pyramid_transition_conformity_face_region_owner_mismatch",
    "v36_source_journal_undo_group_idempotence_entity_reset_checkpoint_replay_mismatch",
    "v36_source_exodus_block_sideset_nodeset_int64_topology_map_export_digest_mismatch",
)


def _with_v36_sweep_transition_journal_exodus_identity(row: dict) -> dict:
    row = _with_v35_metric_curved_sideset_sculpt_identity(row)
    generation = "hex-sweep-contract-231"
    row[
        "hex_sweep_source_target_topology_layer_interval_orientation_jacobian_boundary_volume_mesh_result_generation_identity"
    ] = {
        "sweep_generation": generation,
        **{
            key: generation
            for key in (
                "source_generation",
                "target_generation",
                "layer_generation",
                "orientation_generation",
                "quality_generation",
                "boundary_generation",
                "volume_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "source_surface_topology": "structured_quad_4x4",
        "result_source_surface_topology": "structured_quad_4x4",
        "target_surface_topology": "structured_quad_4x4",
        "result_target_surface_topology": "structured_quad_4x4",
        "source_node_count": 25,
        "result_source_node_count": 25,
        "target_node_count": 25,
        "result_target_node_count": 25,
        "layer_count": 4,
        "result_layer_count": 4,
        "interval_bias": 1.0,
        "result_interval_bias": 1.0,
        "hex_count": 64,
        "result_hex_count": 64,
        "minimum_scaled_jacobian": 0.42,
        "result_minimum_scaled_jacobian": 0.42,
        "minimum_allowed_scaled_jacobian": 0.2,
        "result_minimum_allowed_scaled_jacobian": 0.2,
        "orientation": "source_to_target_positive",
        "result_orientation": "source_to_target_positive",
        "boundary_owners": ["surface:source", "surface:target", "surface:lateral"],
        "result_boundary_owners": ["surface:source", "surface:target", "surface:lateral"],
        "cad_volume_m3": 1.0,
        "mesh_volume_m3": 1.0,
        "result_mesh_volume_m3": 1.0,
        "volume_tolerance_m3": 1.0e-10,
        "result_volume_tolerance_m3": 1.0e-10,
        "sweep_mesh_sha256": "1" * 64,
        "accepted_sweep_mesh_sha256": "1" * 64,
    }

    generation = "mixed-transition-contract-231"
    row[
        "mixed_hex_tet_pyramid_count_interface_conformity_face_orientation_node_region_volume_mesh_result_generation_identity"
    ] = {
        "transition_generation": generation,
        **{
            key: generation
            for key in (
                "hex_generation",
                "tet_generation",
                "pyramid_generation",
                "interface_generation",
                "orientation_generation",
                "region_generation",
                "volume_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "hex_count": 8,
        "result_hex_count": 8,
        "tet_count": 24,
        "result_tet_count": 24,
        "pyramid_count": 6,
        "result_pyramid_count": 6,
        "interface_face_ids": [101, 102, 103, 104, 105, 106],
        "result_interface_face_ids": [101, 102, 103, 104, 105, 106],
        "shared_interface_node_ids": [11, 12, 13, 14, 15, 16, 17, 18],
        "result_shared_interface_node_ids": [11, 12, 13, 14, 15, 16, 17, 18],
        "interface_face_orientation": "hex_outward_pyramid_inward",
        "result_interface_face_orientation": "hex_outward_pyramid_inward",
        "region_labels": ["hex_core", "transition", "tet_shell"],
        "result_region_labels": ["hex_core", "transition", "tet_shell"],
        "signed_region_volumes_m3": [0.4, 0.1, 0.5],
        "result_signed_region_volumes_m3": [0.4, 0.1, 0.5],
        "cad_volume_m3": 1.0,
        "mesh_volume_m3": 1.0,
        "result_mesh_volume_m3": 1.0,
        "volume_tolerance_m3": 1.0e-10,
        "result_volume_tolerance_m3": 1.0e-10,
        "mixed_mesh_owner": "headless:mixed-transition-31",
        "result_mixed_mesh_owner": "headless:mixed-transition-31",
        "mixed_mesh_sha256": "2" * 64,
        "accepted_mixed_mesh_sha256": "2" * 64,
    }

    generation = "journal-replay-contract-231"
    row[
        "journal_undo_idempotence_entity_reset_checkpoint_replay_order_database_result_generation_identity"
    ] = {
        "journal_generation": generation,
        **{
            key: generation
            for key in (
                "undo_generation",
                "idempotence_generation",
                "entity_generation",
                "reset_generation",
                "checkpoint_generation",
                "replay_generation",
                "database_generation",
                "result_generation",
            )
        },
        "command_count": 12,
        "replay_command_count": 12,
        "undo_groups": [[0, 3], [4, 8], [9, 11]],
        "replay_undo_groups": [[0, 3], [4, 8], [9, 11]],
        "entity_ids_after_first_run": [1, 2, 3, 11, 12, 21],
        "entity_ids_after_replay": [1, 2, 3, 11, 12, 21],
        "reset_sequence": ["reset", "create", "mesh", "export"],
        "replay_reset_sequence": ["reset", "create", "mesh", "export"],
        "checkpoint_generation_id": 31,
        "replay_checkpoint_generation_id": 31,
        "command_order_sha256": "3" * 64,
        "replay_command_order_sha256": "3" * 64,
        "first_database_sha256": "4" * 64,
        "replay_database_sha256": "4" * 64,
        "journal_owner": "headless:batch31",
        "replay_journal_owner": "headless:batch31",
        "journal_result_sha256": "5" * 64,
        "accepted_journal_result_sha256": "5" * 64,
    }

    generation = "exodus-semantics-contract-231"
    row[
        "exodus_block_sideset_nodeset_int64_topology_element_map_owner_mesh_export_generation_identity"
    ] = {
        "exodus_generation": generation,
        **{
            key: generation
            for key in (
                "block_generation",
                "sideset_generation",
                "nodeset_generation",
                "integer_generation",
                "topology_generation",
                "map_generation",
                "owner_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "block_ids": [101, 102, 103],
        "exported_block_ids": [101, 102, 103],
        "block_topologies": ["HEX8", "PYRAMID5", "TET4"],
        "exported_block_topologies": ["HEX8", "PYRAMID5", "TET4"],
        "sideset_ids": [201, 202],
        "exported_sideset_ids": [201, 202],
        "sideset_orientation": [1, -1],
        "exported_sideset_orientation": [1, -1],
        "nodeset_ids": [301],
        "exported_nodeset_ids": [301],
        "int64_ids": True,
        "exported_int64_ids": True,
        "maximum_entity_id": 3_000_000_001,
        "exported_maximum_entity_id": 3_000_000_001,
        "element_map": [1001, 1002, 1003, 1004],
        "exported_element_map": [1001, 1002, 1003, 1004],
        "assembly_owner": "assembly:mixed31",
        "exported_assembly_owner": "assembly:mixed31",
        "mesh_generation_id": 31,
        "exported_mesh_generation_id": 31,
        "exodus_mesh_sha256": "6" * 64,
        "exported_exodus_mesh_sha256": "6" * 64,
        "exodus_file_sha256": "7" * 64,
        "accepted_exodus_file_sha256": "7" * 64,
    }
    return row


def test_v36_positive_public_and_source_contracts() -> None:
    row = _with_v36_sweep_transition_journal_exodus_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"


def test_v36_public_hex_sweep_source_target_layers_orientation_jacobian_boundary_owner_mismatch() -> None:
    row = _with_v36_sweep_transition_journal_exodus_identity(summary())
    row[
        "hex_sweep_source_target_topology_layer_interval_orientation_jacobian_boundary_volume_mesh_result_generation_identity"
    ].update({"target_generation": "hex-sweep-contract-230", "result_target_surface_topology": "unstructured_tri", "result_source_node_count": 24, "result_target_node_count": 17, "result_layer_count": 3, "result_interval_bias": -1.0, "result_hex_count": 47, "result_minimum_scaled_jacobian": -0.1, "result_orientation": "target_to_source_negative", "result_boundary_owners": ["surface:stale"], "result_mesh_volume_m3": 0.8, "accepted_sweep_mesh_sha256": "8" * 64})
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["hex_sweeps_use_current_source_target_layers_bias_orientation_jacobian_boundaries_volume_and_mesh"]


def test_v36_public_mixed_hex_tet_pyramid_transition_conformity_face_region_owner_mismatch() -> None:
    row = _with_v36_sweep_transition_journal_exodus_identity(summary())
    row[
        "mixed_hex_tet_pyramid_count_interface_conformity_face_orientation_node_region_volume_mesh_result_generation_identity"
    ].update({"pyramid_generation": "mixed-transition-contract-230", "result_hex_count": 7, "result_tet_count": 25, "result_pyramid_count": 0, "result_interface_face_ids": [101, 101, 103], "result_shared_interface_node_ids": [11, 12, 99], "result_interface_face_orientation": "same_normal_nonconformal", "result_region_labels": ["region0"], "result_signed_region_volumes_m3": [0.4, -0.1, 0.5], "result_mesh_volume_m3": 0.8, "result_mixed_mesh_owner": "gui:stale", "accepted_mixed_mesh_sha256": "9" * 64})
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["mixed_transitions_use_current_counts_interface_nodes_orientation_regions_volume_owner_and_mesh"]


def test_v36_source_journal_undo_group_idempotence_entity_reset_checkpoint_replay_mismatch() -> None:
    row = _with_v36_sweep_transition_journal_exodus_identity(summary())
    row[
        "journal_undo_idempotence_entity_reset_checkpoint_replay_order_database_result_generation_identity"
    ].update({"undo_generation": "journal-replay-contract-230", "replay_command_count": 11, "replay_undo_groups": [[0, 11]], "entity_ids_after_replay": [7, 8, 9], "replay_reset_sequence": ["create", "reset", "export"], "replay_checkpoint_generation_id": 30, "replay_command_order_sha256": "a" * 64, "replay_database_sha256": "b" * 64, "replay_journal_owner": "gui:interactive", "accepted_journal_result_sha256": "c" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["journal_replays_use_current_undo_groups_entities_reset_checkpoint_order_database_owner_and_result"]


def test_v36_source_exodus_block_sideset_nodeset_int64_topology_map_export_digest_mismatch() -> None:
    row = _with_v36_sweep_transition_journal_exodus_identity(summary())
    row[
        "exodus_block_sideset_nodeset_int64_topology_element_map_owner_mesh_export_generation_identity"
    ].update({"block_generation": "exodus-semantics-contract-230", "exported_block_ids": [103, 101], "exported_block_topologies": ["TET10", "HEX20"], "exported_sideset_ids": [202], "exported_sideset_orientation": [1, 1], "exported_nodeset_ids": [], "exported_int64_ids": False, "exported_maximum_entity_id": 2_147_483_647, "exported_element_map": [1004, 1003, 1003], "exported_assembly_owner": "assembly:stale", "exported_mesh_generation_id": 30, "exported_exodus_mesh_sha256": "d" * 64, "accepted_exodus_file_sha256": "e" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["exodus_exports_use_current_blocks_sets_int64_topology_map_owner_mesh_and_file"]


def test_v36_rejects_self_consistent_sweep_volume_loss() -> None:
    row = _with_v36_sweep_transition_journal_exodus_identity(summary())
    identity = row["hex_sweep_source_target_topology_layer_interval_orientation_jacobian_boundary_volume_mesh_result_generation_identity"]
    identity["mesh_volume_m3"] = 0.8
    identity["result_mesh_volume_m3"] = 0.8
    assert _public_result(row)["status"] == "needs_attention"


def test_v36_rejects_self_consistent_duplicate_exodus_map() -> None:
    row = _with_v36_sweep_transition_journal_exodus_identity(summary())
    identity = row["exodus_block_sideset_nodeset_int64_topology_element_map_owner_mesh_export_generation_identity"]
    identity["element_map"] = [1001, 1002, 1002, 1004]
    identity["exported_element_map"] = [1001, 1002, 1002, 1004]
    assert _source_result(row)["status"] == "needs_attention"
