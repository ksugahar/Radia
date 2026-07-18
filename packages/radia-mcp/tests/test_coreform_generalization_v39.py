from __future__ import annotations

from test_coreform_generalization_v37 import _public_result, _source_result, summary
from test_coreform_generalization_v38 import (
    _with_v38_shell_cohesive_virtual_anisotropic_identity,
)


_PERIODIC = (
    "periodic_hex_node_pair_transform_face_orientation_edge_order_jacobian_"
    "block_sideset_database_export_generation_identity"
)
_SWEEP = (
    "thin_sweep_hex_source_target_interval_propagation_layer_thickness_side_"
    "topology_orientation_volume_export_generation_identity"
)
_JOURNAL = (
    "journal_delete_recreate_entity_group_block_sideset_undo_checkpoint_"
    "database_result_generation_identity"
)
_ADAPTIVE = (
    "adaptive_size_field_node_projection_region_boundary_quality_block_"
    "session_export_generation_identity"
)
_PROMOTED_CASE_IDS = (
    "v39_public_periodic_hex_node_pair_transform_face_orientation_jacobian_export_mismatch",
    "v39_public_thin_sweep_hex_source_target_interval_propagation_layer_thickness_mismatch",
    "v39_source_journal_delete_recreate_entity_id_group_block_undo_checkpoint_mismatch",
    "v39_source_adaptive_size_field_node_projection_region_boundary_quality_export_mismatch",
)


def _generations(generation: str, *names: str) -> dict[str, str]:
    return {name: generation for name in names}


def _with_v39_coreform_identity(row: dict) -> dict:
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(row)

    generation = "periodic-hex-271"
    row[_PERIODIC] = {
        "periodic_generation": generation,
        **_generations(
            generation,
            "pair_generation",
            "transform_generation",
            "orientation_generation",
            "edge_order_generation",
            "jacobian_generation",
            "set_generation",
            "database_generation",
            "export_generation",
            "result_generation",
        ),
        "source_surface": "surface:11",
        "result_source_surface": "surface:11",
        "target_surface": "surface:12",
        "result_target_surface": "surface:12",
        "rigid_transform": [
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "result_rigid_transform": [
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "periodic_node_pairs": [[1, 101], [2, 102], [3, 103], [4, 104]],
        "result_periodic_node_pairs": [[1, 101], [2, 102], [3, 103], [4, 104]],
        "source_edge_node_order": [[1, 2], [2, 3], [3, 4], [4, 1]],
        "result_source_edge_node_order": [[1, 2], [2, 3], [3, 4], [4, 1]],
        "target_edge_node_order": [[101, 102], [102, 103], [103, 104], [104, 101]],
        "result_target_edge_node_order": [[101, 102], [102, 103], [103, 104], [104, 101]],
        "face_orientation": "source_outward_target_inward",
        "result_face_orientation": "source_outward_target_inward",
        "minimum_scaled_jacobian": 0.38,
        "result_minimum_scaled_jacobian": 0.38,
        "minimum_allowed_scaled_jacobian": 0.2,
        "result_minimum_allowed_scaled_jacobian": 0.2,
        "block_owner": "block:periodic-21",
        "result_block_owner": "block:periodic-21",
        "sideset_owners": ["sideset:periodic-source", "sideset:periodic-target"],
        "result_sideset_owners": ["sideset:periodic-source", "sideset:periodic-target"],
        "database_generation_id": 271,
        "result_database_generation_id": 271,
        "periodic_export_sha256": "1" * 64,
        "accepted_periodic_export_sha256": "1" * 64,
    }

    generation = "thin-sweep-271"
    row[_SWEEP] = {
        "thin_sweep_generation": generation,
        **_generations(
            generation,
            "source_generation",
            "target_generation",
            "interval_generation",
            "layer_generation",
            "thickness_generation",
            "topology_generation",
            "orientation_generation",
            "volume_generation",
            "export_generation",
            "result_generation",
        ),
        "source_surface_id": 31,
        "result_source_surface_id": 31,
        "target_surface_id": 32,
        "result_target_surface_id": 32,
        "source_interval_count": 16,
        "result_source_interval_count": 16,
        "target_interval_count": 16,
        "result_target_interval_count": 16,
        "propagated_side_intervals": [4, 4, 4, 4],
        "result_propagated_side_intervals": [4, 4, 4, 4],
        "layer_count": 4,
        "result_layer_count": 4,
        "total_thickness_m": 2.0e-3,
        "result_total_thickness_m": 2.0e-3,
        "layer_thickness_m": [5.0e-4] * 4,
        "result_layer_thickness_m": [5.0e-4] * 4,
        "side_topology": ["structured_quad"] * 4,
        "result_side_topology": ["structured_quad"] * 4,
        "element_orientation": "source_to_target_positive",
        "result_element_orientation": "source_to_target_positive",
        "volume_owner": "volume:thin-sweep-31",
        "result_volume_owner": "volume:thin-sweep-31",
        "thin_sweep_export_sha256": "2" * 64,
        "accepted_thin_sweep_export_sha256": "2" * 64,
    }

    generation = "journal-recreate-271"
    row[_JOURNAL] = {
        "journal_generation": generation,
        **_generations(
            generation,
            "delete_generation",
            "recreate_generation",
            "group_generation",
            "block_generation",
            "sideset_generation",
            "undo_generation",
            "checkpoint_generation",
            "database_generation",
            "result_generation",
        ),
        "deleted_entity_ids": [11, 12],
        "replay_deleted_entity_ids": [11, 12],
        "recreated_entity_ids": [21, 22],
        "replay_recreated_entity_ids": [21, 22],
        "group_membership": {"group:recreated": [21, 22]},
        "replay_group_membership": {"group:recreated": [21, 22]},
        "block_membership": {"block:10": [21, 22]},
        "replay_block_membership": {"block:10": [21, 22]},
        "sideset_membership": {"sideset:20": [31, 32]},
        "replay_sideset_membership": {"sideset:20": [31, 32]},
        "undo_depth": 3,
        "replay_undo_depth": 3,
        "checkpoint_generation_id": 271,
        "replay_checkpoint_generation_id": 271,
        "database_owner": "headless:journal-recreate-271",
        "replay_database_owner": "headless:journal-recreate-271",
        "database_sha256": "3" * 64,
        "replay_database_sha256": "3" * 64,
        "journal_result_sha256": "4" * 64,
        "accepted_journal_result_sha256": "4" * 64,
    }

    generation = "adaptive-field-271"
    row[_ADAPTIVE] = {
        "adaptive_generation": generation,
        **_generations(
            generation,
            "field_generation",
            "projection_generation",
            "region_generation",
            "refinement_generation",
            "quality_generation",
            "block_generation",
            "session_generation",
            "export_generation",
            "result_generation",
        ),
        "size_field_samples_m": [1.0e-3, 7.5e-4, 5.0e-4],
        "result_size_field_samples_m": [1.0e-3, 7.5e-4, 5.0e-4],
        "projected_node_ids": [101, 102, 103],
        "result_projected_node_ids": [101, 102, 103],
        "projection_distances_m": [1.0e-6, 2.0e-6, 1.5e-6],
        "result_projection_distances_m": [1.0e-6, 2.0e-6, 1.5e-6],
        "maximum_projection_distance_m": 5.0e-6,
        "result_maximum_projection_distance_m": 5.0e-6,
        "region_boundary_ids": [31, 32, 33],
        "result_region_boundary_ids": [31, 32, 33],
        "refinement_generation_id": 271,
        "result_refinement_generation_id": 271,
        "quality_histogram": [2, 12, 36, 50],
        "result_quality_histogram": [2, 12, 36, 50],
        "block_map": {"volume:1": "block:10", "volume:2": "block:20"},
        "result_block_map": {"volume:1": "block:10", "volume:2": "block:20"},
        "session_owner": "headless:adaptive-field-271",
        "result_session_owner": "headless:adaptive-field-271",
        "adaptive_export_sha256": "5" * 64,
        "accepted_adaptive_export_sha256": "5" * 64,
    }
    return row


def test_v39_positive_public_and_source_contracts() -> None:
    row = _with_v39_coreform_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"


def test_v39_public_periodic_hex_node_pair_transform_face_orientation_jacobian_export_mismatch() -> None:
    row = _with_v39_coreform_identity(summary())
    row[_PERIODIC].update(
        {
            "pair_generation": "periodic-hex-270",
            "edge_order_generation": "periodic-hex-269",
            "result_generation": "periodic-hex-268",
            "result_periodic_node_pairs": [[1, 104], [2, 103], [3, 102], [4, 101]],
            "result_rigid_transform": [[1.0, 0.0, 0.0, -1.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
            "result_source_edge_node_order": [[2, 1]],
            "result_target_edge_node_order": [[104, 103]],
            "result_face_orientation": "same_outward_normal",
            "result_minimum_scaled_jacobian": -0.1,
            "result_block_owner": "block:old",
            "result_sideset_owners": ["sideset:old"],
            "result_database_generation_id": 270,
            "accepted_periodic_export_sha256": "9" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "periodic_hex_replays_use_current_pairs_transform_edge_order_orientation_quality_sets_database_and_export"
    ]


def test_v39_public_thin_sweep_hex_source_target_interval_propagation_layer_thickness_mismatch() -> None:
    row = _with_v39_coreform_identity(summary())
    row[_SWEEP].update(
        {
            "interval_generation": "thin-sweep-270",
            "thickness_generation": "thin-sweep-269",
            "result_generation": "thin-sweep-268",
            "result_source_surface_id": 32,
            "result_target_surface_id": 31,
            "result_source_interval_count": 8,
            "result_target_interval_count": 12,
            "result_propagated_side_intervals": [2, 3],
            "result_layer_count": 3,
            "result_total_thickness_m": 3.0e-3,
            "result_layer_thickness_m": [1.0e-3, -1.0e-3, 3.0e-3],
            "result_side_topology": ["unstructured_tri"],
            "result_element_orientation": "target_to_source_negative",
            "result_volume_owner": "volume:old",
            "accepted_thin_sweep_export_sha256": "a" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "thin_sweeps_use_current_surfaces_intervals_layers_thickness_topology_orientation_owner_and_export"
    ]


def test_v39_source_journal_delete_recreate_entity_id_group_block_undo_checkpoint_mismatch() -> None:
    row = _with_v39_coreform_identity(summary())
    row[_JOURNAL].update(
        {
            "recreate_generation": "journal-recreate-270",
            "checkpoint_generation": "journal-recreate-269",
            "result_generation": "journal-recreate-268",
            "replay_deleted_entity_ids": [12, 13],
            "replay_recreated_entity_ids": [11, 12],
            "replay_group_membership": {"group:old": [11]},
            "replay_block_membership": {"block:old": [11]},
            "replay_sideset_membership": {},
            "replay_undo_depth": 1,
            "replay_checkpoint_generation_id": 270,
            "replay_database_owner": "gui:old",
            "replay_database_sha256": "b" * 64,
            "accepted_journal_result_sha256": "c" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "journal_recreate_replays_use_current_entities_memberships_undo_checkpoint_database_and_result"
    ]


def test_v39_source_adaptive_size_field_node_projection_region_boundary_quality_export_mismatch() -> None:
    row = _with_v39_coreform_identity(summary())
    row[_ADAPTIVE].update(
        {
            "field_generation": "adaptive-field-270",
            "projection_generation": "adaptive-field-269",
            "result_generation": "adaptive-field-268",
            "result_size_field_samples_m": [1.0e-3, -7.5e-4],
            "result_projected_node_ids": [103, 102, 101],
            "result_projection_distances_m": [8.0e-6],
            "result_maximum_projection_distance_m": 1.0e-6,
            "result_region_boundary_ids": [31, 99],
            "result_refinement_generation_id": 270,
            "result_quality_histogram": [50, 36, 12, 2],
            "result_block_map": {"volume:1": "block:old"},
            "result_session_owner": "gui:old",
            "accepted_adaptive_export_sha256": "d" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "adaptive_mesh_replays_use_current_size_field_projection_region_quality_blocks_session_and_export"
    ]


def test_v39_rejects_self_consistent_periodic_edge_map_mismatch() -> None:
    row = _with_v39_coreform_identity(summary())
    row[_PERIODIC]["target_edge_node_order"] = [[102, 101], [103, 102], [104, 103], [101, 104]]
    row[_PERIODIC]["result_target_edge_node_order"] = row[_PERIODIC]["target_edge_node_order"]
    assert _public_result(row)["status"] == "needs_attention"


def test_v39_rejects_self_consistent_thin_sweep_thickness_sum_mismatch() -> None:
    row = _with_v39_coreform_identity(summary())
    row[_SWEEP]["total_thickness_m"] = 3.0e-3
    row[_SWEEP]["result_total_thickness_m"] = 3.0e-3
    assert _public_result(row)["status"] == "needs_attention"


def test_v39_rejects_self_consistent_recreated_entity_reuse() -> None:
    row = _with_v39_coreform_identity(summary())
    for key in ("recreated_entity_ids", "replay_recreated_entity_ids"):
        row[_JOURNAL][key] = [11, 22]
    for key in ("group_membership", "replay_group_membership"):
        row[_JOURNAL][key] = {"group:recreated": [11, 22]}
    for key in ("block_membership", "replay_block_membership"):
        row[_JOURNAL][key] = {"block:10": [11, 22]}
    assert _source_result(row)["status"] == "needs_attention"


def test_v39_rejects_self_consistent_nonmonotone_adaptive_size_field() -> None:
    row = _with_v39_coreform_identity(summary())
    row[_ADAPTIVE]["size_field_samples_m"] = [1.0e-3, 1.2e-3, 5.0e-4]
    row[_ADAPTIVE]["result_size_field_samples_m"] = row[_ADAPTIVE]["size_field_samples_m"]
    assert _source_result(row)["status"] == "needs_attention"
