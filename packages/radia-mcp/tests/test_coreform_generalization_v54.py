from copy import deepcopy

from radia_mcp.cubit.topology_transaction_identity_v54 import EXODUS, HEX, JOURNAL, TRANSITION, validate_public_identity, validate_source_identity


CASE_IDS = {
    "v54_public_hex_sweep_edgecorrespondence_highorder_jacobian_block_owner_mismatch",
    "v54_public_tethex_pyramid_transition_facediagonal_orientation_quality_owner_mismatch",
    "v54_source_tool_journal_undo_transaction_entityallocator_checkpoint_owner_mismatch",
    "v54_source_tool_exodus_sideset_distributionfactor_topology_qa_owner_mismatch",
}


def _records() -> dict[str, object]:
    generation = "coreform-v54-test"
    generations = lambda names: {name: generation for name in names}
    correspondence = [[1, 101], [2, 102], [3, 103], [4, 104]]
    transitions = [{"element": 201, "base_nodes": [1, 2, 3, 4], "face_diagonal": [1, 3], "orientation": 1, "scaled_jacobian": 0.34}]
    topology = [{"element": 101, "side": 2, "topology": "HEX8"}, {"element": 205, "side": 4, "topology": "TET4"}]
    qa = {"application": "headless-export", "version": "v54", "date": "2026-07-20", "time": "04:50:00"}
    return {
        HEX: {"generation": generation, **generations(("edge_generation", "jacobian_generation", "block_generation", "owner_generation", "result_generation")), "edge_correspondence": correspondence, "result_edge_correspondence": correspondence, "high_order_jacobian_samples": [0.62, 0.71, 0.68, 0.59], "result_high_order_jacobian_samples": [0.62, 0.71, 0.68, 0.59], "element_block": "block:hex-v54", "result_element_block": "block:hex-v54", "mesh_owner": "headless:hex-v54", "result_mesh_owner": "headless:hex-v54", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64},
        TRANSITION: {"generation": generation, **generations(("diagonal_generation", "orientation_generation", "quality_generation", "topology_generation", "owner_generation", "result_generation")), "pyramid_transitions": transitions, "result_pyramid_transitions": transitions, "interface_topology": "tet-pyramid-hex", "result_interface_topology": "tet-pyramid-hex", "mesh_owner": "headless:transition-v54", "result_mesh_owner": "headless:transition-v54", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64},
        JOURNAL: {"generation": generation, **generations(("transaction_generation", "allocator_generation", "checkpoint_generation", "revision_generation", "owner_generation", "result_generation")), "transaction_id": "transaction:42", "replayed_transaction_id": "transaction:42", "undo_depth": 1, "replayed_undo_depth": 1, "entity_allocator_before": 208, "replayed_entity_allocator_before": 208, "entity_allocator_after": 208, "replayed_entity_allocator_after": 208, "command_checkpoint": 73, "replayed_command_checkpoint": 73, "database_revision": "database:v54-r7", "replayed_database_revision": "database:v54-r7", "journal_owner": "headless:journal-v54", "replayed_journal_owner": "headless:journal-v54", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64},
        EXODUS: {"generation": generation, **generations(("sideset_generation", "factor_generation", "topology_generation", "qa_generation", "owner_generation", "result_generation")), "sideset_topology": topology, "replayed_sideset_topology": topology, "distribution_factors": [1.0, 0.75], "replayed_distribution_factors": [1.0, 0.75], "qa_record": qa, "replayed_qa_record": qa, "mesh_revision": "mesh:v54-r2", "replayed_mesh_revision": "mesh:v54-r2", "export_owner": "headless:exodus-v54", "replayed_export_owner": "headless:exodus-v54", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64},
    }


def test_v54_positive_identities_are_accepted() -> None:
    assert validate_public_identity(_records())["status"] == "ok"
    assert validate_source_identity(_records())["status"] == "ok"


def test_v54_frozen_mutations_are_rejected() -> None:
    value = deepcopy(_records())
    value[HEX].update({"result_edge_correspondence": [[1, 104]], "result_high_order_jacobian_samples": [-0.2], "result_element_block": "block:stale", "result_mesh_owner": "headless:stale"})
    value[TRANSITION].update({"result_pyramid_transitions": [{"element": 201, "base_nodes": [4, 3, 2, 1], "face_diagonal": [2, 4], "orientation": -1, "scaled_jacobian": -0.1}], "result_interface_topology": "tet-hex-direct", "result_mesh_owner": "headless:stale"})
    value[JOURNAL].update({"replayed_transaction_id": "transaction:41", "replayed_undo_depth": 0, "replayed_entity_allocator_after": 211, "replayed_command_checkpoint": 70, "replayed_database_revision": "database:stale", "replayed_journal_owner": "headless:stale"})
    value[EXODUS].update({"replayed_sideset_topology": [{"element": 101, "side": 7, "topology": "HEX8"}], "replayed_distribution_factors": [9.0], "replayed_qa_record": {"application": "stale"}, "replayed_mesh_revision": "mesh:stale", "replayed_export_owner": "headless:stale"})
    assert validate_public_identity(value)["status"] == "needs_attention"
    assert validate_source_identity(value)["status"] == "needs_attention"


def test_v54_self_consistent_invalid_topology_and_transaction_are_rejected() -> None:
    value = deepcopy(_records())
    bad_transition = [{"element": 201, "base_nodes": [1, 2, 3, 4], "face_diagonal": [1, 2], "orientation": -1, "scaled_jacobian": -0.1}]
    value[TRANSITION]["pyramid_transitions"] = value[TRANSITION]["result_pyramid_transitions"] = bad_transition
    value[JOURNAL]["entity_allocator_after"] = value[JOURNAL]["replayed_entity_allocator_after"] = 211
    value[EXODUS]["distribution_factors"] = value[EXODUS]["replayed_distribution_factors"] = [1.0]
    assert validate_public_identity(value)["status"] == "needs_attention"
    assert validate_source_identity(value)["status"] == "needs_attention"


def test_v54_malformed_mesh_values_are_rejected_without_raising() -> None:
    value = deepcopy(_records())
    malformed_transition = [{"element": 201, "base_nodes": [[1], 2, 3, 4], "face_diagonal": [[1], 3], "orientation": 1, "scaled_jacobian": 0.4}]
    value[TRANSITION]["pyramid_transitions"] = value[TRANSITION]["result_pyramid_transitions"] = malformed_transition
    malformed_topology = [{"element": [101], "side": 2, "topology": ["HEX8"]}]
    value[EXODUS]["sideset_topology"] = value[EXODUS]["replayed_sideset_topology"] = malformed_topology
    assert validate_public_identity(value)["status"] == "needs_attention"
    assert validate_source_identity(value)["status"] == "needs_attention"
