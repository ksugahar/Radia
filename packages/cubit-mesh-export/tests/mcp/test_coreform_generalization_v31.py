from __future__ import annotations

import json

from test_coreform_generalization_v30 import (
    _with_v30_sweep_transition_journal_sideset_identity,
    cubit_conformal_hex_pyramid_tet_interface_gate,
    cubit_mixed_transition_source_gate,
    summary,
)


_PROMOTED_CASE_IDS = (
    "v31_public_mixed_hex_tet_pyramid_transition_face_conformity_jacobian_orientation_mismatch",
    "v31_public_periodic_high_order_mesh_affine_node_pair_edge_face_owner_mismatch",
    "v31_source_boolean_imprint_merge_sideset_stable_entity_lineage_measure_mismatch",
    "v31_source_exodus_high_order_block_midnode_qa_record_int64_restart_owner_mismatch",
)


def _with_v31_high_order_periodic_lineage_exodus_identity(row: dict) -> dict:
    row = _with_v30_sweep_transition_journal_sideset_identity(row)
    generation = "curved-transition-181"
    row[
        "mixed_high_order_transition_midnode_permutation_parametric_face_quadrature_jacobian_export_generation_identity"
    ] = {
        "transition_generation": generation,
        "order_transition_generation": generation,
        "midnode_transition_generation": generation,
        "permutation_transition_generation": generation,
        "parametric_transition_generation": generation,
        "jacobian_transition_generation": generation,
        "export_transition_generation": generation,
        "result_transition_generation": generation,
        "element_topologies": ["hex27", "pyramid14", "tet10"],
        "result_element_topologies": ["hex27", "pyramid14", "tet10"],
        "polynomial_orders": [2, 2, 2],
        "result_polynomial_orders": [2, 2, 2],
        "interface_corner_node_ids": [[1, 2, 3, 4], [1, 2, 5]],
        "result_interface_corner_node_ids": [[1, 2, 3, 4], [1, 2, 5]],
        "interface_midnode_ids": [[11, 12, 13, 14], [11, 15, 16]],
        "result_interface_midnode_ids": [[11, 12, 13, 14], [11, 15, 16]],
        "face_node_permutations": [
            [0, 3, 2, 1, 7, 6, 5, 4],
            [0, 2, 1, 5, 4, 3],
        ],
        "result_face_node_permutations": [
            [0, 3, 2, 1, 7, 6, 5, 4],
            [0, 2, 1, 5, 4, 3],
        ],
        "parametric_face_coordinate_sha256": "1" * 64,
        "result_parametric_face_coordinate_sha256": "1" * 64,
        "quadrature_scaled_jacobians": [
            [0.42, 0.38, 0.35],
            [0.31, 0.29, 0.27],
        ],
        "result_quadrature_scaled_jacobians": [
            [0.42, 0.38, 0.35],
            [0.31, 0.29, 0.27],
        ],
        "export_connectivity_sha256": "2" * 64,
        "accepted_export_connectivity_sha256": "2" * 64,
    }
    generation = "periodic-high-order-181"
    transform = [
        [-1.0, 0.0, 0.0, 1.0],
        [0.0, -1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    row[
        "periodic_high_order_affine_corner_edge_face_node_pair_orientation_sideset_residual_generation_identity"
    ] = {
        "periodic_generation": generation,
        "affine_periodic_generation": generation,
        "corner_periodic_generation": generation,
        "edge_periodic_generation": generation,
        "face_periodic_generation": generation,
        "orientation_periodic_generation": generation,
        "sideset_periodic_generation": generation,
        "residual_periodic_generation": generation,
        "result_periodic_generation": generation,
        "polynomial_order": 3,
        "result_polynomial_order": 3,
        "affine_transform_4x4": transform,
        "result_affine_transform_4x4": transform,
        "corner_node_pairs": [[1, 101], [2, 102]],
        "result_corner_node_pairs": [[1, 101], [2, 102]],
        "edge_node_pairs": [[11, 111], [12, 112], [13, 113]],
        "result_edge_node_pairs": [[11, 111], [12, 112], [13, 113]],
        "face_node_pairs": [[21, 121], [22, 122]],
        "result_face_node_pairs": [[21, 121], [22, 122]],
        "orientation_signs": [-1, -1],
        "result_orientation_signs": [-1, -1],
        "source_sideset_id": 100,
        "result_source_sideset_id": 100,
        "target_sideset_id": 200,
        "result_target_sideset_id": 200,
        "maximum_periodic_residual_m": 1.0e-13,
        "result_maximum_periodic_residual_m": 1.0e-13,
        "periodic_residual_tolerance_m": 1.0e-10,
        "periodic_mesh_sha256": "3" * 64,
        "result_periodic_mesh_sha256": "3" * 64,
    }
    generation = "imprint-lineage-181"
    row[
        "boolean_imprint_merge_entity_lineage_surface_orientation_measure_adjacency_journal_digest_generation_identity"
    ] = {
        "lineage_generation": generation,
        "imprint_lineage_generation": generation,
        "merge_lineage_generation": generation,
        "orientation_lineage_generation": generation,
        "measure_lineage_generation": generation,
        "adjacency_lineage_generation": generation,
        "journal_lineage_generation": generation,
        "result_lineage_generation": generation,
        "pre_entity_ids": [10, 20],
        "resolved_pre_entity_ids": [10, 20],
        "post_entity_ids": [101, 102, 201],
        "resolved_post_entity_ids": [101, 102, 201],
        "entity_lineage_pairs": [[10, 101], [10, 102], [20, 201]],
        "resolved_entity_lineage_pairs": [[10, 101], [10, 102], [20, 201]],
        "surface_orientation_signs": [1, 1, -1],
        "resolved_surface_orientation_signs": [1, 1, -1],
        "surface_measures_m2": [0.4, 0.6, 1.0],
        "resolved_surface_measures_m2": [0.4, 0.6, 1.0],
        "block_adjacency": [[101, 30], [102, 30], [201, 40]],
        "resolved_block_adjacency": [[101, 30], [102, 30], [201, 40]],
        "journal_generation_id": "journal-181",
        "resolved_journal_generation_id": "journal-181",
        "lineage_table_sha256": "4" * 64,
        "resolved_lineage_table_sha256": "4" * 64,
        "model_digest_sha256": "5" * 64,
        "accepted_model_digest_sha256": "5" * 64,
    }
    generation = "exodus-high-order-181"
    row[
        "exodus_high_order_block_midnode_order_int64_qa_restart_sideset_digest_generation_identity"
    ] = {
        "exodus_generation": generation,
        "block_exodus_generation": generation,
        "midnode_exodus_generation": generation,
        "id_exodus_generation": generation,
        "qa_exodus_generation": generation,
        "restart_exodus_generation": generation,
        "sideset_exodus_generation": generation,
        "digest_exodus_generation": generation,
        "result_exodus_generation": generation,
        "block_ids": [3000000001, 3000000002],
        "decoded_block_ids": [3000000001, 3000000002],
        "element_topologies": ["HEX27", "PYRAMID14"],
        "decoded_element_topologies": ["HEX27", "PYRAMID14"],
        "midnode_orderings": [[9, 10, 11, 12], [6, 7, 8, 9]],
        "decoded_midnode_orderings": [[9, 10, 11, 12], [6, 7, 8, 9]],
        "integer_word_size_bits": 64,
        "decoded_integer_word_size_bits": 64,
        "qa_records": [["Coreform Cubit", "2026.6", "2026-07-18", "12:00:00"]],
        "decoded_qa_records": [
            ["Coreform Cubit", "2026.6", "2026-07-18", "12:00:00"]
        ],
        "restart_step_index": 4,
        "decoded_restart_step_index": 4,
        "sideset_owner_block_ids": [
            [100, 3000000001],
            [200, 3000000002],
        ],
        "decoded_sideset_owner_block_ids": [
            [100, 3000000001],
            [200, 3000000002],
        ],
        "exodus_file_sha256": "6" * 64,
        "decoded_exodus_file_sha256": "6" * 64,
        "connectivity_table_sha256": "7" * 64,
        "decoded_connectivity_table_sha256": "7" * 64,
    }
    return row


def test_v31_positive_high_order_periodic_lineage_and_exodus() -> None:
    row = _with_v31_high_order_periodic_lineage_exodus_identity(summary())
    assert json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))["status"] == "ok"
    assert json.loads(cubit_mixed_transition_source_gate(row))["status"] == "ok"


def test_v31_public_mixed_hex_tet_pyramid_transition_face_conformity_jacobian_orientation_mismatch() -> None:
    row = _with_v31_high_order_periodic_lineage_exodus_identity(summary())
    row[
        "mixed_high_order_transition_midnode_permutation_parametric_face_quadrature_jacobian_export_generation_identity"
    ].update(
        {
            "midnode_transition_generation": "curved-transition-180",
            "jacobian_transition_generation": "curved-transition-179",
            "result_transition_generation": "curved-transition-178",
            "result_element_topologies": ["hex20", "pyramid13", "tet10"],
            "result_polynomial_orders": [2, 1, 2],
            "result_interface_corner_node_ids": [[1, 2, 4, 3]],
            "result_interface_midnode_ids": [[11, 13, 12, 14], [11, 16, 15]],
            "result_face_node_permutations": [[0, 1, 2, 3, 4, 5, 6, 7]],
            "result_parametric_face_coordinate_sha256": "8" * 64,
            "result_quadrature_scaled_jacobians": [[0.42, -0.02, 0.35]],
            "accepted_export_connectivity_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "curved_mixed_transitions_use_current_orders_midnodes_permutations_parametric_faces_quadrature_jacobians_and_export"
    ]


def test_v31_public_periodic_high_order_mesh_affine_node_pair_edge_face_owner_mismatch() -> None:
    row = _with_v31_high_order_periodic_lineage_exodus_identity(summary())
    row[
        "periodic_high_order_affine_corner_edge_face_node_pair_orientation_sideset_residual_generation_identity"
    ].update(
        {
            "affine_periodic_generation": "periodic-high-order-180",
            "face_periodic_generation": "periodic-high-order-179",
            "result_periodic_generation": "periodic-high-order-178",
            "result_polynomial_order": 2,
            "result_affine_transform_4x4": [[1.0, 0.0, 0.0, 0.0]],
            "result_corner_node_pairs": [[1, 102], [2, 101]],
            "result_edge_node_pairs": [[11, 111], [13, 112]],
            "result_face_node_pairs": [[21, 122]],
            "result_orientation_signs": [1, -1],
            "result_source_sideset_id": 200,
            "result_target_sideset_id": 100,
            "result_maximum_periodic_residual_m": 1.0e-4,
            "result_periodic_mesh_sha256": "a" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "periodic_high_order_meshes_use_current_affine_corner_edge_face_pairs_orientation_sidesets_and_residual"
    ]


def test_v31_source_boolean_imprint_merge_sideset_stable_entity_lineage_measure_mismatch() -> None:
    row = _with_v31_high_order_periodic_lineage_exodus_identity(summary())
    row[
        "boolean_imprint_merge_entity_lineage_surface_orientation_measure_adjacency_journal_digest_generation_identity"
    ].update(
        {
            "imprint_lineage_generation": "imprint-lineage-180",
            "measure_lineage_generation": "imprint-lineage-179",
            "result_lineage_generation": "imprint-lineage-178",
            "resolved_pre_entity_ids": [20, 10],
            "resolved_post_entity_ids": [101, 201],
            "resolved_entity_lineage_pairs": [[10, 101], [20, 102]],
            "resolved_surface_orientation_signs": [1, -1, -1],
            "resolved_surface_measures_m2": [0.4, 0.5, 1.0],
            "resolved_block_adjacency": [[101, 40], [201, 30]],
            "resolved_journal_generation_id": "journal-180",
            "resolved_lineage_table_sha256": "b" * 64,
            "accepted_model_digest_sha256": "c" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "boolean_imprint_merges_use_current_entity_lineage_orientation_measure_adjacency_journal_and_digest"
    ]


def test_v31_source_exodus_high_order_block_midnode_qa_record_int64_restart_owner_mismatch() -> None:
    row = _with_v31_high_order_periodic_lineage_exodus_identity(summary())
    row[
        "exodus_high_order_block_midnode_order_int64_qa_restart_sideset_digest_generation_identity"
    ].update(
        {
            "block_exodus_generation": "exodus-high-order-180",
            "qa_exodus_generation": "exodus-high-order-179",
            "result_exodus_generation": "exodus-high-order-178",
            "decoded_block_ids": [3000000002, 3000000001],
            "decoded_element_topologies": ["HEX20", "PYRAMID13"],
            "decoded_midnode_orderings": [[12, 11, 10, 9], [6, 8, 7]],
            "decoded_integer_word_size_bits": 32,
            "decoded_qa_records": [["unknown", "old", "", ""]],
            "decoded_restart_step_index": 3,
            "decoded_sideset_owner_block_ids": [[100, 3000000002]],
            "decoded_exodus_file_sha256": "d" * 64,
            "decoded_connectivity_table_sha256": "e" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "high_order_exodus_uses_current_blocks_midnodes_int64_qa_restart_sideset_owners_and_digests"
    ]
