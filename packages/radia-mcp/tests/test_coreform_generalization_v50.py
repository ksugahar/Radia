from copy import deepcopy

from radia_mcp.cubit.topology_replay_identity_v50 import validate_public_identity, validate_source_identity


PROMOTED_CASE_IDS = {
    "v50_public_hex_node_order_face_adjacency_orientation_jacobian_mesh_owner_mismatch",
    "v50_public_block_sideset_entity_dimension_overlap_duplicate_membership_owner_mismatch",
    "v50_source_tool_journal_include_relative_path_aprepro_scope_expansion_digest_owner_mismatch",
    "v50_source_tool_save_restore_session_entity_id_block_sideset_mesh_revision_owner_mismatch",
}


def _records() -> dict[str, object]:
    generation = "v50-coreform"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    dimensions = {"block:steel": 3, "block:air": 3, "sideset:wall": 2, "sideset:symmetry": 2}
    memberships = {
        "block:steel": [11, 12],
        "block:air": [21, 22],
        "sideset:wall": [101, 102],
        "sideset:symmetry": [103, 104],
    }
    entity_ids = {"volume:body": 1, "surface:wall": 11, "surface:symmetry": 12}
    return {
        "hex_node_order_face_adjacency_orientation_jacobian_mesh_owner_identity": {
            "generation": generation,
            "node_order_generation": generation,
            "face_adjacency_generation": generation,
            "orientation_generation": generation,
            "jacobian_generation": generation,
            "owner_generation": generation,
            "result_generation": generation,
            "hex_node_order": [101, 102, 103, 104, 105, 106, 107, 108],
            "result_hex_node_order": [101, 102, 103, 104, 105, 106, 107, 108],
            "face_adjacency": {"face:bottom": "boundary:source", "face:top": "hex:2002"},
            "result_face_adjacency": {"face:bottom": "boundary:source", "face:top": "hex:2002"},
            "face_orientation_signs": [1, 1, 1, 1, 1, 1],
            "result_face_orientation_signs": [1, 1, 1, 1, 1, 1],
            "minimum_scaled_jacobian": 0.38,
            "result_minimum_scaled_jacobian": 0.38,
            "mesh_owner": "headless:hex-v50",
            "result_mesh_owner": "headless:hex-v50",
            **result,
        },
        "block_sideset_entity_dimension_overlap_duplicate_membership_owner_identity": {
            "generation": generation,
            "dimension_generation": generation,
            "membership_generation": generation,
            "overlap_generation": generation,
            "duplicate_generation": generation,
            "owner_generation": generation,
            "result_generation": generation,
            "entity_dimensions": dimensions,
            "result_entity_dimensions": dimensions,
            "group_memberships": memberships,
            "result_group_memberships": memberships,
            "allowed_overlaps": [],
            "result_allowed_overlaps": [],
            "duplicate_memberships": [],
            "result_duplicate_memberships": [],
            "group_owner": "headless:groups-v50",
            "result_group_owner": "headless:groups-v50",
            **result,
        },
        "journal_include_relative_path_aprepro_scope_expansion_digest_owner_identity": {
            "generation": generation,
            "include_generation": generation,
            "scope_generation": generation,
            "expansion_generation": generation,
            "digest_generation": generation,
            "owner_generation": generation,
            "result_generation": generation,
            "relative_include_paths": ["include/materials.jou", "include/mesh_controls.jou"],
            "result_relative_include_paths": ["include/materials.jou", "include/mesh_controls.jou"],
            "aprepro_scope": {"mesh_size": "0.25", "layers": "8"},
            "result_aprepro_scope": {"mesh_size": "0.25", "layers": "8"},
            "expanded_commands": ["volume 1 size 0.25", "volume 1 interval 8"],
            "result_expanded_commands": ["volume 1 size 0.25", "volume 1 interval 8"],
            "expanded_journal_sha256": "d" * 64,
            "result_expanded_journal_sha256": "d" * 64,
            "journal_owner": "headless:journal-v50",
            "result_journal_owner": "headless:journal-v50",
            **result,
        },
        "save_restore_session_entity_id_block_sideset_mesh_revision_owner_identity": {
            "generation": generation,
            "session_generation": generation,
            "entity_generation": generation,
            "block_generation": generation,
            "sideset_generation": generation,
            "mesh_generation": generation,
            "owner_generation": generation,
            "result_generation": generation,
            "session_id": "session:headless-v50",
            "result_session_id": "session:headless-v50",
            "entity_ids": entity_ids,
            "result_entity_ids": entity_ids,
            "blocks": {"block:steel": [1]},
            "result_blocks": {"block:steel": [1]},
            "sidesets": {"sideset:wall": [11], "sideset:symmetry": [12]},
            "result_sidesets": {"sideset:wall": [11], "sideset:symmetry": [12]},
            "mesh_revision": "mesh-revision:v50-r7",
            "result_mesh_revision": "mesh-revision:v50-r7",
            "model_owner": "headless:model-v50",
            "result_model_owner": "headless:model-v50",
            **result,
        },
    }


def test_v50_positive_public_and_source_replays_are_accepted() -> None:
    records = _records()
    assert validate_public_identity(records)["status"] == "ok"
    assert validate_source_identity(records)["status"] == "ok"


def test_v50_public_mutations_are_rejected() -> None:
    records = deepcopy(_records())
    records["hex_node_order_face_adjacency_orientation_jacobian_mesh_owner_identity"]["result_minimum_scaled_jacobian"] = -0.08
    records["block_sideset_entity_dimension_overlap_duplicate_membership_owner_identity"]["result_duplicate_memberships"] = [11]
    result = validate_public_identity(records)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v50_source_mutations_are_rejected() -> None:
    records = deepcopy(_records())
    records["journal_include_relative_path_aprepro_scope_expansion_digest_owner_identity"]["result_expanded_journal_sha256"] = "a" * 64
    records["save_restore_session_entity_id_block_sideset_mesh_revision_owner_identity"]["result_mesh_revision"] = "mesh-revision:old"
    result = validate_source_identity(records)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v50_self_consistent_invalid_topology_and_paths_are_rejected() -> None:
    records = deepcopy(_records())
    hex_row = records["hex_node_order_face_adjacency_orientation_jacobian_mesh_owner_identity"]
    hex_row["hex_node_order"] = hex_row["result_hex_node_order"] = [101, 101, 103, 104, 105, 106, 107, 108]
    journal = records["journal_include_relative_path_aprepro_scope_expansion_digest_owner_identity"]
    journal["relative_include_paths"] = journal["result_relative_include_paths"] = ["../private/materials.jou"]
    assert validate_public_identity(records)["status"] == "needs_attention"
    assert validate_source_identity(records)["status"] == "needs_attention"


def test_v50_self_consistent_group_overlap_and_restore_id_drift_are_rejected() -> None:
    records = deepcopy(_records())
    group = records["block_sideset_entity_dimension_overlap_duplicate_membership_owner_identity"]
    group["group_memberships"] = group["result_group_memberships"] = {
        "block:steel": [11, 12],
        "block:air": [12, 21],
        "sideset:wall": [101, 102],
        "sideset:symmetry": [103, 104],
    }
    restore = records["save_restore_session_entity_id_block_sideset_mesh_revision_owner_identity"]
    restore["blocks"] = restore["result_blocks"] = {"block:steel": [99]}
    assert validate_public_identity(records)["status"] == "needs_attention"
    assert validate_source_identity(records)["status"] == "needs_attention"
