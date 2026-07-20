from copy import deepcopy

from radia_mcp.cubit.sweep_export_identity_v53 import EXODUS, HEX_SWEEP, JOURNAL, PYRAMID, validate_public_identity, validate_source_identity


CASE_IDS = {
    "v53_public_hex_sweep_source_target_face_layer_orientation_volume_owner_mismatch",
    "v53_public_pyramid_transition_apex_basequad_conformity_jacobian_owner_mismatch",
    "v53_source_tool_exodus_sideset_element_side_ordinal_topology_block_owner_mismatch",
    "v53_source_tool_journal_include_aprepro_variable_expansion_workdir_owner_mismatch",
}


def _records() -> dict[str, object]:
    generation = "coreform-v53-test"
    generations = lambda names: {name: generation for name in names}
    entries = [{"element_id": 101, "side_ordinal": 2, "topology": "HEX8", "block_id": 1}]
    includes = ["units.apr", "geometry.jou", "mesh.jou"]
    hashes = {name: str(index) * 64 for index, name in enumerate(includes, 4)}
    return {
        HEX_SWEEP: {"generation": generation, **generations(("face_generation", "layer_generation", "orientation_generation", "volume_generation", "owner_generation", "result_generation")), "source_surface": "surface:11", "result_source_surface": "surface:11", "target_surface": "surface:12", "result_target_surface": "surface:12", "source_quad_count": 24, "result_source_quad_count": 24, "target_quad_count": 24, "result_target_quad_count": 24, "layer_count": 8, "result_layer_count": 8, "orientation": "source_to_target", "result_orientation": "source_to_target", "volume_id": "volume:1", "result_volume_id": "volume:1", "volume_owner": "headless:sweep-v53", "result_volume_owner": "headless:sweep-v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64},
        PYRAMID: {"generation": generation, **generations(("apex_generation", "base_generation", "conformity_generation", "jacobian_generation", "owner_generation", "result_generation")), "apex_node": 101, "result_apex_node": 101, "base_quad_nodes": [1, 2, 3, 4], "result_base_quad_nodes": [1, 2, 3, 4], "conformal_neighbor_faces": 5, "result_conformal_neighbor_faces": 5, "minimum_scaled_jacobian": 0.31, "result_minimum_scaled_jacobian": 0.31, "transition_owner": "headless:pyramid-v53", "result_transition_owner": "headless:pyramid-v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64},
        EXODUS: {"generation": generation, **generations(("sideset_generation", "element_generation", "ordinal_generation", "topology_generation", "block_generation", "owner_generation", "result_generation")), "sideset_name": "pressure", "replayed_sideset_name": "pressure", "sideset_entries": entries, "replayed_sideset_entries": entries, "export_owner": "headless:exodus-v53", "replayed_export_owner": "headless:exodus-v53", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64},
        JOURNAL: {"generation": generation, **generations(("include_generation", "variable_generation", "expansion_generation", "workdir_generation", "owner_generation", "result_generation")), "include_order": includes, "replayed_include_order": includes, "include_sha256": hashes, "replayed_include_sha256": hashes, "aprepro_variables": {"gap_mm": 0.8}, "replayed_aprepro_variables": {"gap_mm": 0.8}, "expanded_variables": {"gap_m": 0.0008}, "replayed_expanded_variables": {"gap_m": 0.0008}, "working_directory": "workspace:coreform/project", "replayed_working_directory": "workspace:coreform/project", "journal_owner": "headless:journal-v53", "replayed_journal_owner": "headless:journal-v53", "result_sha256": "7" * 64, "accepted_result_sha256": "7" * 64},
    }


def test_v53_positive_identities_are_accepted() -> None:
    assert validate_public_identity(_records())["status"] == "ok"
    assert validate_source_identity(_records())["status"] == "ok"


def test_v53_frozen_mutations_are_rejected() -> None:
    value = deepcopy(_records())
    value[HEX_SWEEP].update({"result_target_quad_count": 20, "result_volume_owner": "headless:stale"})
    value[PYRAMID].update({"result_minimum_scaled_jacobian": -0.1, "result_transition_owner": "headless:stale"})
    value[EXODUS].update({"replayed_sideset_entries": [{"element_id": 101, "side_ordinal": 7, "topology": "HEX8", "block_id": 1}], "replayed_export_owner": "headless:stale"})
    value[JOURNAL].update({"replayed_include_order": list(reversed(value[JOURNAL]["include_order"])), "replayed_working_directory": "workspace:other", "replayed_journal_owner": "headless:stale"})
    assert validate_public_identity(value)["status"] == "needs_attention"
    assert validate_source_identity(value)["status"] == "needs_attention"


def test_v53_self_consistent_invalid_topology_and_paths_are_rejected() -> None:
    value = deepcopy(_records())
    value[PYRAMID]["conformal_neighbor_faces"] = value[PYRAMID]["result_conformal_neighbor_faces"] = 4
    bad = [{"element_id": 101, "side_ordinal": 7, "topology": "HEX8", "block_id": 1}]
    value[EXODUS]["sideset_entries"] = value[EXODUS]["replayed_sideset_entries"] = bad
    value[JOURNAL]["working_directory"] = value[JOURNAL]["replayed_working_directory"] = "workspace:coreform/../private"
    assert validate_public_identity(value)["status"] == "needs_attention"
    assert validate_source_identity(value)["status"] == "needs_attention"
