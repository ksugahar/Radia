from copy import deepcopy

from radia_mcp.cubit.topology_replay_identity_v49 import validate_public_identity, validate_source_identity


PROMOTED_CASE_IDS = {
    "v49_public_hex_sweep_source_target_topology_layer_interval_bias_periodic_owner_mismatch",
    "v49_public_pyramid_transition_apex_orientation_jacobian_boundary_block_owner_mismatch",
    "v49_source_tool_journal_undo_checkpoint_entity_allocator_replay_cursor_revision_owner_mismatch",
    "v49_source_tool_netgen_export_order_curved_node_boundary_name_index_digest_owner_mismatch",
}


def _records() -> dict[str, object]:
    generation = "v49-coreform"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    return {
        "hex_sweep_source_target_topology_layer_interval_bias_periodic_owner_identity": {
            "generation": generation, "topology_generation": generation, "layer_generation": generation,
            "interval_generation": generation, "periodic_generation": generation, "result_generation": generation,
            "source_target_topology": {"source": "surface:11", "target": "surface:19", "side_count": 8},
            "result_source_target_topology": {"source": "surface:11", "target": "surface:19", "side_count": 8},
            "layer_count": 12, "result_layer_count": 12, "interval_bias": 1.15, "result_interval_bias": 1.15,
            "periodic_node_pairs": [[101, 201], [102, 202]], "result_periodic_node_pairs": [[101, 201], [102, 202]],
            "mesh_owner": "headless:hex-v49", "result_mesh_owner": "headless:hex-v49", **result,
        },
        "pyramid_transition_apex_orientation_jacobian_boundary_block_owner_identity": {
            "generation": generation, "apex_generation": generation, "orientation_generation": generation,
            "jacobian_generation": generation, "boundary_generation": generation, "result_generation": generation,
            "apex_node_id": 501, "result_apex_node_id": 501,
            "face_orientations": [1, 1, -1, 1], "result_face_orientations": [1, 1, -1, 1],
            "minimum_scaled_jacobian": 0.42, "result_minimum_scaled_jacobian": 0.42,
            "boundary_blocks": {"hex_side": 21, "tet_side": 22}, "result_boundary_blocks": {"hex_side": 21, "tet_side": 22},
            "mesh_owner": "headless:pyramid-v49", "result_mesh_owner": "headless:pyramid-v49", **result,
        },
        "journal_undo_checkpoint_entity_allocator_replay_cursor_revision_owner_identity": {
            "generation": generation, "checkpoint_generation": generation, "allocator_generation": generation,
            "cursor_generation": generation, "revision_generation": generation, "result_generation": generation,
            "undo_checkpoint": "checkpoint:before-webcut", "result_undo_checkpoint": "checkpoint:before-webcut",
            "entity_id_allocator": {"next_vertex": 301, "next_surface": 41},
            "result_entity_id_allocator": {"next_vertex": 301, "next_surface": 41},
            "replay_cursor": 184, "result_replay_cursor": 184, "model_revision": "model-r12", "result_model_revision": "model-r12",
            "journal_owner": "headless:journal-v49", "result_journal_owner": "headless:journal-v49", **result,
        },
        "netgen_export_order_curved_node_boundary_name_index_digest_owner_identity": {
            "generation": generation, "order_generation": generation, "curved_node_generation": generation,
            "boundary_generation": generation, "index_generation": generation, "file_generation": generation, "result_generation": generation,
            "element_order": 3, "result_element_order": 3,
            "curved_node_map": {"edge:11": [1001, 1002]}, "result_curved_node_map": {"edge:11": [1001, 1002]},
            "boundary_names": {"1": "air", "2": "conductor"}, "result_boundary_names": {"1": "air", "2": "conductor"},
            "index_base": 1, "result_index_base": 1, "export_sha256": "e" * 64, "result_export_sha256": "e" * 64,
            "export_owner": "headless:netgen-v49", "result_export_owner": "headless:netgen-v49", **result,
        },
    }


def test_v49_positive_public_and_source_replays_are_accepted() -> None:
    records = _records()
    assert validate_public_identity(records)["status"] == "ok"
    assert validate_source_identity(records)["status"] == "ok"


def test_v49_public_topology_mutations_are_rejected() -> None:
    records = deepcopy(_records())
    records["hex_sweep_source_target_topology_layer_interval_bias_periodic_owner_identity"]["result_layer_count"] = 10
    records["pyramid_transition_apex_orientation_jacobian_boundary_block_owner_identity"]["result_minimum_scaled_jacobian"] = -0.05
    result = validate_public_identity(records)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v49_source_replay_mutations_are_rejected() -> None:
    records = deepcopy(_records())
    records["journal_undo_checkpoint_entity_allocator_replay_cursor_revision_owner_identity"]["result_replay_cursor"] = 177
    records["netgen_export_order_curved_node_boundary_name_index_digest_owner_identity"]["result_index_base"] = 0
    result = validate_source_identity(records)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2

