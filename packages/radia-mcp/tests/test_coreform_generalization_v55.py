from copy import deepcopy

from radia_mcp.cubit.periodic_boundary_identity_v55 import (
    BOUNDARY_LAYER,
    EXODUS,
    MERGE,
    PERIODIC,
    validate_public_identity,
    validate_source_identity,
)


CASE_IDS = {
    "v55_public_periodic_hex_nodeequivalence_transform_highorderjacobian_block_owner_mismatch",
    "v55_public_boundarylayer_hex_thickness_growth_cornercollapse_quality_owner_mismatch",
    "v55_source_tool_merge_tolerance_entityprovenance_idremap_group_owner_mismatch",
    "v55_source_tool_exodus_block_attribute_truth_table_elementorder_owner_mismatch",
}


def _records() -> dict[str, object]:
    generation = "coreform-v55-test"
    generations = lambda names: {name: generation for name in names}
    pairs = [[1, 101], [2, 102], [3, 103], [4, 104]]
    transform = [
        [1.0, 0.0, 0.0, 0.1],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    attributes = {
        "block:1": {"material_id": 7, "section": "solid"},
        "block:2": {"material_id": 8, "section": "solid"},
    }
    return {
        PERIODIC: {
            "generation": generation,
            **generations(("equivalence_generation", "transform_generation", "jacobian_generation", "block_generation", "owner_generation", "result_generation")),
            "node_equivalence": pairs,
            "result_node_equivalence": pairs,
            "periodic_transform_4x4": transform,
            "result_periodic_transform_4x4": transform,
            "high_order_jacobian_samples": [0.58, 0.63, 0.61, 0.55],
            "result_high_order_jacobian_samples": [0.58, 0.63, 0.61, 0.55],
            "element_block": "block:periodic-hex-v55",
            "result_element_block": "block:periodic-hex-v55",
            "mesh_owner": "headless:periodic-hex-v55",
            "result_mesh_owner": "headless:periodic-hex-v55",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        BOUNDARY_LAYER: {
            "generation": generation,
            **generations(("thickness_generation", "growth_generation", "topology_generation", "collapse_generation", "quality_generation", "owner_generation", "result_generation")),
            "first_layer_thickness_m": 1.0e-3,
            "result_first_layer_thickness_m": 1.0e-3,
            "growth_ratio": 1.2,
            "result_growth_ratio": 1.2,
            "layer_count": 5,
            "result_layer_count": 5,
            "corner_topology": ["corner:convex-1", "corner:concave-2"],
            "result_corner_topology": ["corner:convex-1", "corner:concave-2"],
            "collapsed_layer_count": 0,
            "result_collapsed_layer_count": 0,
            "minimum_scaled_jacobian": 0.31,
            "result_minimum_scaled_jacobian": 0.31,
            "mesh_owner": "headless:boundary-layer-v55",
            "result_mesh_owner": "headless:boundary-layer-v55",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
        MERGE: {
            "generation": generation,
            **generations(("tolerance_generation", "provenance_generation", "remap_generation", "group_generation", "revision_generation", "owner_generation", "result_generation")),
            "merge_tolerance_m": 1.0e-6,
            "replayed_merge_tolerance_m": 1.0e-6,
            "source_entity_provenance": {"surface:1": "cad:A", "surface:2": "cad:B"},
            "replayed_source_entity_provenance": {"surface:1": "cad:A", "surface:2": "cad:B"},
            "entity_id_remap": {"surface:2": "surface:1"},
            "replayed_entity_id_remap": {"surface:2": "surface:1"},
            "group_membership": {"group:interface": ["surface:1"]},
            "replayed_group_membership": {"group:interface": ["surface:1"]},
            "database_revision": "database:v55-r4",
            "replayed_database_revision": "database:v55-r4",
            "merge_owner": "headless:merge-v55",
            "replayed_merge_owner": "headless:merge-v55",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
        EXODUS: {
            "generation": generation,
            **generations(("attribute_generation", "truth_table_generation", "order_generation", "qa_generation", "owner_generation", "result_generation")),
            "block_attributes": attributes,
            "replayed_block_attributes": attributes,
            "variable_truth_table": {"stress": [1, 1], "temperature": [1, 0]},
            "replayed_variable_truth_table": {"stress": [1, 1], "temperature": [1, 0]},
            "element_order": {"block:1": "HEX27", "block:2": "TET10"},
            "replayed_element_order": {"block:1": "HEX27", "block:2": "TET10"},
            "qa_revision": "qa:v55-r2",
            "replayed_qa_revision": "qa:v55-r2",
            "file_owner": "headless:exodus-v55",
            "replayed_file_owner": "headless:exodus-v55",
            "result_sha256": "4" * 64,
            "accepted_result_sha256": "4" * 64,
        },
    }


def test_v55_positive_identities_are_accepted() -> None:
    assert validate_public_identity(_records())["status"] == "ok"
    assert validate_source_identity(_records())["status"] == "ok"


def test_v55_frozen_mutations_are_rejected() -> None:
    value = deepcopy(_records())
    value[PERIODIC].update({"result_node_equivalence": [[1, 104]], "result_periodic_transform_4x4": [[1.0, 0.0], [0.0, 1.0]], "result_high_order_jacobian_samples": [-0.2], "result_element_block": "block:stale", "result_mesh_owner": "headless:stale"})
    value[BOUNDARY_LAYER].update({"result_first_layer_thickness_m": 2.0e-3, "result_growth_ratio": 1.8, "result_corner_topology": ["corner:collapsed"], "result_collapsed_layer_count": 2, "result_minimum_scaled_jacobian": -0.1, "result_mesh_owner": "headless:stale"})
    value[MERGE].update({"replayed_merge_tolerance_m": 1.0e-3, "replayed_source_entity_provenance": {"surface:9": "unknown"}, "replayed_entity_id_remap": {"surface:9": "surface:8"}, "replayed_group_membership": {"group:other": ["surface:9"]}, "replayed_database_revision": "database:stale", "replayed_merge_owner": "headless:stale"})
    value[EXODUS].update({"replayed_block_attributes": {"block:9": {"material_id": 99}}, "replayed_variable_truth_table": {"stress": [1, 0, 1]}, "replayed_element_order": {"block:1": "HEX8"}, "replayed_qa_revision": "qa:stale", "replayed_file_owner": "headless:stale"})
    assert validate_public_identity(value)["status"] == "needs_attention"
    assert validate_source_identity(value)["status"] == "needs_attention"


def test_v55_self_consistent_singular_transform_or_collapsed_layer_is_rejected() -> None:
    value = deepcopy(_records())
    singular = [[1.0, 0.0, 0.0, 0.1], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    value[PERIODIC]["periodic_transform_4x4"] = value[PERIODIC]["result_periodic_transform_4x4"] = singular
    value[BOUNDARY_LAYER]["collapsed_layer_count"] = value[BOUNDARY_LAYER]["result_collapsed_layer_count"] = 1
    assert validate_public_identity(value)["status"] == "needs_attention"


def test_v55_self_consistent_dangling_group_or_truth_table_is_rejected() -> None:
    value = deepcopy(_records())
    value[MERGE]["group_membership"] = value[MERGE]["replayed_group_membership"] = {"group:interface": ["surface:2"]}
    value[EXODUS]["variable_truth_table"] = value[EXODUS]["replayed_variable_truth_table"] = {"stress": [1, 0, 1]}
    assert validate_source_identity(value)["status"] == "needs_attention"


def test_v55_numeric_sha256_values_are_rejected() -> None:
    value = _records()
    numeric_digest = int("9" * 64)
    for row in value.values():
        row["result_sha256"] = numeric_digest
        row["accepted_result_sha256"] = numeric_digest
    assert validate_public_identity(value)["status"] == "needs_attention"
    assert validate_source_identity(value)["status"] == "needs_attention"
