from __future__ import annotations

import json

from test_coreform_generalization_v33 import (
    _with_v33_topology_transaction_roundtrip_identity,
    cubit_conformal_hex_pyramid_tet_interface_gate,
    cubit_mixed_transition_source_gate,
    summary,
)


_PROMOTED_CASE_IDS = (
    "v34_public_hex20_hex27_serendipity_lagrange_node_order_jacobian_quadrature_mismatch",
    "v34_public_sheet_midplane_thickness_normal_block_sideset_mass_property_mismatch",
    "v34_source_sweep_source_target_interval_bias_match_periodic_layer_journal_mismatch",
    "v34_source_checkpoint_restore_partition_ghost_owner_entity_id_quality_digest_mismatch",
)


def _with_v34_high_order_sheet_sweep_checkpoint_identity(row: dict) -> dict:
    row = _with_v33_topology_transaction_roundtrip_identity(row)
    generation = "high-order-hex-211"
    row[
        "hex20_hex27_family_node_role_reference_order_jacobian_volume_face_mesh_result_generation_identity"
    ] = {
        "high_order_hex_generation": generation,
        **{
            key: generation
            for key in (
                "family_generation",
                "node_generation",
                "reference_generation",
                "jacobian_generation",
                "volume_generation",
                "face_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "element_ids": [20, 27],
        "result_element_ids": [20, 27],
        "element_families": ["hex20_serendipity", "hex27_lagrange"],
        "result_element_families": ["hex20_serendipity", "hex27_lagrange"],
        "node_counts": [20, 27],
        "result_node_counts": [20, 27],
        "corner_node_order": list(range(8)),
        "result_corner_node_order": list(range(8)),
        "edge_node_order": list(range(8, 20)),
        "result_edge_node_order": list(range(8, 20)),
        "hex27_face_node_order": list(range(20, 26)),
        "result_hex27_face_node_order": list(range(20, 26)),
        "hex27_center_node": 26,
        "result_hex27_center_node": 26,
        "minimum_jacobians": [0.42, 0.39],
        "result_minimum_jacobians": [0.42, 0.39],
        "quadrature_volumes_m3": [1.0, 1.0],
        "result_quadrature_volumes_m3": [1.0, 1.0],
        "geometric_volumes_m3": [1.0, 1.0],
        "result_geometric_volumes_m3": [1.0, 1.0],
        "curved_face_owners": [[20, 101], [27, 102]],
        "result_curved_face_owners": [[20, 101], [27, 102]],
        "high_order_mesh_sha256": "1" * 64,
        "accepted_high_order_mesh_sha256": "1" * 64,
    }

    generation = "sheet-midplane-211"
    row[
        "sheet_midplane_source_offset_thickness_normal_block_sideset_area_mass_geometry_result_generation_identity"
    ] = {
        "sheet_generation": generation,
        **{
            key: generation
            for key in (
                "source_generation",
                "midplane_generation",
                "thickness_generation",
                "normal_generation",
                "block_generation",
                "sideset_generation",
                "mass_generation",
                "geometry_generation",
                "result_generation",
            )
        },
        "source_volume_id": 11,
        "result_source_volume_id": 11,
        "midplane_offset_m": 0.0,
        "result_midplane_offset_m": 0.0,
        "thickness_m": 0.002,
        "result_thickness_m": 0.002,
        "normal_orientation": "outward_source_volume",
        "result_normal_orientation": "outward_source_volume",
        "shell_block_id": 301,
        "result_shell_block_id": 301,
        "top_sideset_id": 401,
        "result_top_sideset_id": 401,
        "bottom_sideset_id": 402,
        "result_bottom_sideset_id": 402,
        "midplane_area_m2": 2.0,
        "result_midplane_area_m2": 2.0,
        "source_volume_m3": 0.004,
        "result_source_volume_m3": 0.004,
        "density_kg_m3": 7800.0,
        "result_density_kg_m3": 7800.0,
        "shell_mass_kg": 31.2,
        "result_shell_mass_kg": 31.2,
        "sheet_geometry_sha256": "2" * 64,
        "accepted_sheet_geometry_sha256": "2" * 64,
    }

    generation = "sweep-replay-211"
    row[
        "sweep_source_target_interval_bias_match_periodic_layer_scheme_journal_result_generation_identity"
    ] = {
        "sweep_generation": generation,
        **{
            key: generation
            for key in (
                "source_generation",
                "target_generation",
                "interval_generation",
                "bias_generation",
                "match_generation",
                "periodic_generation",
                "layer_generation",
                "journal_generation",
                "result_generation",
            )
        },
        "volume_id": 21,
        "replayed_volume_id": 21,
        "source_surface_id": 101,
        "replayed_source_surface_id": 101,
        "target_surface_id": 102,
        "replayed_target_surface_id": 102,
        "interval_count": 8,
        "replayed_interval_count": 8,
        "bias_ratio": 2.0,
        "replayed_bias_ratio": 2.0,
        "bias_direction": "source_to_target",
        "replayed_bias_direction": "source_to_target",
        "match_pairs": [[201, 301], [202, 302]],
        "replayed_match_pairs": [[201, 301], [202, 302]],
        "periodic_layer_map": [[0, 8], [1, 7], [2, 6], [3, 5], [4, 4]],
        "replayed_periodic_layer_map": [[0, 8], [1, 7], [2, 6], [3, 5], [4, 4]],
        "node_layer_count": 9,
        "replayed_node_layer_count": 9,
        "volume_scheme": "sweep",
        "replayed_volume_scheme": "sweep",
        "journal_sha256": "3" * 64,
        "replayed_journal_sha256": "3" * 64,
        "sweep_result_sha256": "4" * 64,
        "accepted_sweep_result_sha256": "4" * 64,
    }

    generation = "checkpoint-partition-211"
    row[
        "checkpoint_partition_owned_ghost_persistent_block_sideset_quality_model_result_generation_identity"
    ] = {
        "checkpoint_generation": generation,
        **{
            key: generation
            for key in (
                "partition_generation",
                "owned_generation",
                "ghost_generation",
                "persistent_generation",
                "block_generation",
                "sideset_generation",
                "quality_generation",
                "model_generation",
                "result_generation",
            )
        },
        "partition_count": 4,
        "restored_partition_count": 4,
        "owned_entity_ids": [[0, [1, 2]], [1, [3, 4]], [2, [5, 6]], [3, [7, 8]]],
        "restored_owned_entity_ids": [[0, [1, 2]], [1, [3, 4]], [2, [5, 6]], [3, [7, 8]]],
        "ghost_entity_ids": [[0, [3]], [1, [2, 5]], [2, [4, 7]], [3, [6]]],
        "restored_ghost_entity_ids": [[0, [3]], [1, [2, 5]], [2, [4, 7]], [3, [6]]],
        "persistent_entity_ids": list(range(1, 9)),
        "restored_persistent_entity_ids": list(range(1, 9)),
        "block_membership": [[101, [1, 2, 3, 4]], [102, [5, 6, 7, 8]]],
        "restored_block_membership": [[101, [1, 2, 3, 4]], [102, [5, 6, 7, 8]]],
        "sideset_membership": [[201, [1, 8]]],
        "restored_sideset_membership": [[201, [1, 8]]],
        "partition_minimum_scaled_jacobian": [0.41, 0.38, 0.44, 0.36],
        "restored_partition_minimum_scaled_jacobian": [0.41, 0.38, 0.44, 0.36],
        "checkpoint_model_sha256": "5" * 64,
        "restored_checkpoint_model_sha256": "5" * 64,
        "checkpoint_result_sha256": "6" * 64,
        "accepted_checkpoint_result_sha256": "6" * 64,
    }
    return row


def _public_result(row: dict) -> dict:
    return json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))


def _source_result(row: dict) -> dict:
    return json.loads(cubit_mixed_transition_source_gate(row))


def test_v34_positive_high_order_sheet_sweep_and_checkpoint() -> None:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"


def test_v34_public_hex20_hex27_serendipity_lagrange_node_order_jacobian_quadrature_mismatch() -> None:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(summary())
    row[
        "hex20_hex27_family_node_role_reference_order_jacobian_volume_face_mesh_result_generation_identity"
    ].update(
        {
            "family_generation": "high-order-hex-210",
            "node_generation": "high-order-hex-209",
            "result_generation": "high-order-hex-208",
            "result_element_ids": [27, 20],
            "result_element_families": ["hex20_lagrange", "hex27_serendipity"],
            "result_node_counts": [27, 20],
            "result_corner_node_order": [0, 2, 1, 3, 4, 5, 7, 6],
            "result_edge_node_order": list(reversed(range(8, 20))),
            "result_hex27_face_node_order": [20, 21, 22, 23, 24, 26],
            "result_hex27_center_node": 25,
            "result_minimum_jacobians": [0.42, -0.1],
            "result_quadrature_volumes_m3": [1.0, 0.8],
            "result_geometric_volumes_m3": [1.0, 1.2],
            "result_curved_face_owners": [[20, 102], [27, 101]],
            "accepted_high_order_mesh_sha256": "7" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "high_order_hexes_use_current_family_node_roles_reference_order_jacobian_volume_faces_and_mesh"
    ]


def test_v34_public_sheet_midplane_thickness_normal_block_sideset_mass_property_mismatch() -> None:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(summary())
    row[
        "sheet_midplane_source_offset_thickness_normal_block_sideset_area_mass_geometry_result_generation_identity"
    ].update(
        {
            "source_generation": "sheet-midplane-210",
            "mass_generation": "sheet-midplane-209",
            "result_generation": "sheet-midplane-208",
            "result_source_volume_id": 12,
            "result_midplane_offset_m": 0.001,
            "result_thickness_m": 0.004,
            "result_normal_orientation": "inward_source_volume",
            "result_shell_block_id": 302,
            "result_top_sideset_id": 402,
            "result_bottom_sideset_id": 401,
            "result_midplane_area_m2": 1.5,
            "result_source_volume_m3": 0.003,
            "result_density_kg_m3": 2700.0,
            "result_shell_mass_kg": 100.0,
            "accepted_sheet_geometry_sha256": "8" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "sheet_midplanes_use_current_source_offset_thickness_normal_sets_area_mass_and_geometry"
    ]


def test_v34_source_sweep_source_target_interval_bias_match_periodic_layer_journal_mismatch() -> None:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(summary())
    row[
        "sweep_source_target_interval_bias_match_periodic_layer_scheme_journal_result_generation_identity"
    ].update(
        {
            "source_generation": "sweep-replay-210",
            "journal_generation": "sweep-replay-209",
            "result_generation": "sweep-replay-208",
            "replayed_volume_id": 22,
            "replayed_source_surface_id": 102,
            "replayed_target_surface_id": 101,
            "replayed_interval_count": 6,
            "replayed_bias_ratio": 0.5,
            "replayed_bias_direction": "target_to_source",
            "replayed_match_pairs": [[201, 302], [202, 301]],
            "replayed_periodic_layer_map": [[0, 7], [1, 6]],
            "replayed_node_layer_count": 8,
            "replayed_volume_scheme": "map",
            "replayed_journal_sha256": "9" * 64,
            "accepted_sweep_result_sha256": "a" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "sweep_replays_use_current_source_target_intervals_bias_matches_periodic_layers_journal_and_result"
    ]


def test_v34_source_checkpoint_restore_partition_ghost_owner_entity_id_quality_digest_mismatch() -> None:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(summary())
    row[
        "checkpoint_partition_owned_ghost_persistent_block_sideset_quality_model_result_generation_identity"
    ].update(
        {
            "partition_generation": "checkpoint-partition-210",
            "model_generation": "checkpoint-partition-209",
            "result_generation": "checkpoint-partition-208",
            "restored_partition_count": 3,
            "restored_owned_entity_ids": [[0, [1, 2]], [1, [2, 4]], [2, [5, 6]], [3, [7, 9]]],
            "restored_ghost_entity_ids": [[0, [99]], [1, [2]], [2, [4]], [3, [6]]],
            "restored_persistent_entity_ids": [1, 2, 3, 4, 5, 6, 7, 9],
            "restored_block_membership": [[101, [1, 2]], [102, [4, 5]]],
            "restored_sideset_membership": [[201, [99]]],
            "restored_partition_minimum_scaled_jacobian": [0.41, -0.02, 0.44],
            "restored_checkpoint_model_sha256": "b" * 64,
            "accepted_checkpoint_result_sha256": "c" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "checkpoint_restores_use_current_partitions_owned_ghost_persistent_sets_quality_model_and_result"
    ]


def test_v34_rejects_self_consistent_hex27_face_center_role_swap() -> None:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(summary())
    identity = row[
        "hex20_hex27_family_node_role_reference_order_jacobian_volume_face_mesh_result_generation_identity"
    ]
    identity["hex27_face_node_order"] = [20, 21, 22, 23, 24, 26]
    identity["result_hex27_face_node_order"] = [20, 21, 22, 23, 24, 26]
    identity["hex27_center_node"] = 25
    identity["result_hex27_center_node"] = 25
    assert _public_result(row)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_shell_mass_without_geometric_closure() -> None:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(summary())
    identity = row[
        "sheet_midplane_source_offset_thickness_normal_block_sideset_area_mass_geometry_result_generation_identity"
    ]
    identity["shell_mass_kg"] = 30.0
    identity["result_shell_mass_kg"] = 30.0
    assert _public_result(row)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_sweep_layer_count_mismatch() -> None:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(summary())
    identity = row[
        "sweep_source_target_interval_bias_match_periodic_layer_scheme_journal_result_generation_identity"
    ]
    identity["node_layer_count"] = 8
    identity["replayed_node_layer_count"] = 8
    assert _source_result(row)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_duplicate_partition_ownership() -> None:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(summary())
    identity = row[
        "checkpoint_partition_owned_ghost_persistent_block_sideset_quality_model_result_generation_identity"
    ]
    duplicate = [[0, [1, 2]], [1, [2, 4]], [2, [5, 6]], [3, [7, 8]]]
    identity["owned_entity_ids"] = duplicate
    identity["restored_owned_entity_ids"] = duplicate
    assert _source_result(row)["status"] == "needs_attention"
