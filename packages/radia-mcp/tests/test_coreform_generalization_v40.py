from __future__ import annotations

from test_coreform_generalization_v37 import _public_result, _source_result, summary
from test_coreform_generalization_v39 import _generations, _with_v39_coreform_identity


_MEDIAL = (
    "medial_axis_hex_decomposition_sheet_pair_thickness_topology_interval_"
    "quality_block_export_generation_identity"
)
_CHAIN = (
    "curve_chain_interval_bias_corner_parity_boundary_layer_orientation_"
    "sideset_export_generation_identity"
)
_SMOOTH = (
    "smoothing_iteration_constraint_node_motion_quality_history_checkpoint_"
    "database_result_generation_identity"
)
_NAMED = "named_entity_metadata_group_transfer_save_open_export_owner_generation_identity"
_PROMOTED_CASE_IDS = (
    "v40_public_medial_axis_hex_decomposition_sheet_pair_thickness_topology_quality_export_mismatch",
    "v40_public_curve_chain_interval_bias_corner_parity_hex_boundary_layer_export_mismatch",
    "v40_source_smoothing_iteration_constraint_node_motion_quality_history_checkpoint_mismatch",
    "v40_source_named_entity_metadata_group_transfer_save_open_export_owner_mismatch",
)


def _with_v40_coreform_identity(row: dict) -> dict:
    row = _with_v39_coreform_identity(row)

    generation = "medial-axis-hex-311"
    row[_MEDIAL] = {
        "medial_axis_generation": generation,
        **_generations(
            generation,
            "sheet_generation",
            "thickness_generation",
            "decomposition_generation",
            "topology_generation",
            "interval_generation",
            "quality_generation",
            "block_generation",
            "export_generation",
            "result_generation",
        ),
        "paired_sheet_ids": [[11, 12], [21, 22]],
        "result_paired_sheet_ids": [[11, 12], [21, 22]],
        "local_thickness_m": [0.010, 0.014],
        "result_local_thickness_m": [0.010, 0.014],
        "decomposition_cells": {
            "cell:1": [11, 12, 31, 32],
            "cell:2": [21, 22, 32, 33],
        },
        "result_decomposition_cells": {
            "cell:1": [11, 12, 31, 32],
            "cell:2": [21, 22, 32, 33],
        },
        "shared_topology_faces": [32],
        "result_shared_topology_faces": [32],
        "interval_counts": {"curve:41": 8, "curve:42": 8, "curve:43": 6},
        "result_interval_counts": {"curve:41": 8, "curve:42": 8, "curve:43": 6},
        "interval_parity": "compatible_even",
        "result_interval_parity": "compatible_even",
        "minimum_scaled_jacobian": 0.42,
        "result_minimum_scaled_jacobian": 0.42,
        "minimum_allowed_scaled_jacobian": 0.20,
        "result_minimum_allowed_scaled_jacobian": 0.20,
        "block_owner": "block:medial-axis-hex-31",
        "result_block_owner": "block:medial-axis-hex-31",
        "medial_axis_export_sha256": "6" * 64,
        "accepted_medial_axis_export_sha256": "6" * 64,
    }

    generation = "curve-chain-hex-311"
    row[_CHAIN] = {
        "curve_chain_generation": generation,
        **_generations(
            generation,
            "chain_generation",
            "interval_generation",
            "bias_generation",
            "corner_generation",
            "boundary_layer_generation",
            "orientation_generation",
            "sideset_generation",
            "export_generation",
            "result_generation",
        ),
        "curve_chain_order": [101, 102, 103, 104],
        "result_curve_chain_order": [101, 102, 103, 104],
        "chain_orientation": "counterclockwise",
        "result_chain_orientation": "counterclockwise",
        "interval_counts": [8, 10, 8, 10],
        "result_interval_counts": [8, 10, 8, 10],
        "bias_directions": ["forward", "reverse", "forward", "reverse"],
        "result_bias_directions": ["forward", "reverse", "forward", "reverse"],
        "corner_interval_sums": [18, 18, 18, 18],
        "result_corner_interval_sums": [18, 18, 18, 18],
        "boundary_layer_thickness_m": [5.0e-4, 7.5e-4, 1.0e-3],
        "result_boundary_layer_thickness_m": [5.0e-4, 7.5e-4, 1.0e-3],
        "total_boundary_layer_thickness_m": 2.25e-3,
        "result_total_boundary_layer_thickness_m": 2.25e-3,
        "element_orientation": "outward_positive",
        "result_element_orientation": "outward_positive",
        "sideset_owner": "sideset:boundary-layer-41",
        "result_sideset_owner": "sideset:boundary-layer-41",
        "curve_chain_export_sha256": "7" * 64,
        "accepted_curve_chain_export_sha256": "7" * 64,
    }

    generation = "smooth-replay-311"
    row[_SMOOTH] = {
        "smoothing_generation": generation,
        **_generations(
            generation,
            "algorithm_generation",
            "iteration_generation",
            "constraint_generation",
            "motion_generation",
            "quality_generation",
            "checkpoint_generation",
            "database_generation",
            "result_generation",
        ),
        "algorithm": "smart_laplacian",
        "replay_algorithm": "smart_laplacian",
        "iteration_count": 12,
        "replay_iteration_count": 12,
        "fixed_node_ids": [1, 2, 3, 4],
        "replay_fixed_node_ids": [1, 2, 3, 4],
        "fixed_node_displacement_m": [0.0, 0.0, 0.0, 0.0],
        "replay_fixed_node_displacement_m": [0.0, 0.0, 0.0, 0.0],
        "moved_node_ids": [101, 102, 103],
        "replay_moved_node_ids": [101, 102, 103],
        "node_displacement_m": [2.0e-4, 1.5e-4, 1.0e-4],
        "replay_node_displacement_m": [2.0e-4, 1.5e-4, 1.0e-4],
        "maximum_allowed_displacement_m": 5.0e-4,
        "replay_maximum_allowed_displacement_m": 5.0e-4,
        "quality_history": [0.18, 0.25, 0.31, 0.37],
        "replay_quality_history": [0.18, 0.25, 0.31, 0.37],
        "checkpoint_generation_id": 311,
        "replay_checkpoint_generation_id": 311,
        "database_owner": "headless:smooth-replay-311",
        "replay_database_owner": "headless:smooth-replay-311",
        "database_sha256": "8" * 64,
        "replay_database_sha256": "8" * 64,
        "smoothing_result_sha256": "9" * 64,
        "accepted_smoothing_result_sha256": "9" * 64,
    }

    generation = "named-entity-replay-311"
    row[_NAMED] = {
        "named_entity_generation": generation,
        **_generations(
            generation,
            "name_generation",
            "metadata_generation",
            "group_generation",
            "set_generation",
            "save_generation",
            "open_generation",
            "export_generation",
            "result_generation",
        ),
        "entity_names": {"volume:1": "rotor", "surface:11": "airgap_master"},
        "replay_entity_names": {"volume:1": "rotor", "surface:11": "airgap_master"},
        "metadata_attributes": {
            "volume:1": {"material": "steel", "frame": "rotating"},
            "surface:11": {"role": "periodic_master"},
        },
        "replay_metadata_attributes": {
            "volume:1": {"material": "steel", "frame": "rotating"},
            "surface:11": {"role": "periodic_master"},
        },
        "group_membership": {"group:motor": ["volume:1", "surface:11"]},
        "replay_group_membership": {"group:motor": ["volume:1", "surface:11"]},
        "block_membership": {"block:10": ["volume:1"]},
        "replay_block_membership": {"block:10": ["volume:1"]},
        "sideset_membership": {"sideset:20": ["surface:11"]},
        "replay_sideset_membership": {"sideset:20": ["surface:11"]},
        "save_generation_id": 311,
        "replay_save_generation_id": 311,
        "open_generation_id": 311,
        "replay_open_generation_id": 311,
        "export_owner": "headless:named-entity-replay-311",
        "replay_export_owner": "headless:named-entity-replay-311",
        "database_sha256": "a" * 64,
        "replay_database_sha256": "a" * 64,
        "named_entity_export_sha256": "b" * 64,
        "accepted_named_entity_export_sha256": "b" * 64,
    }
    return row


def test_v40_positive_public_and_source_contracts() -> None:
    row = _with_v40_coreform_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"


def test_v40_public_medial_axis_mismatch() -> None:
    row = _with_v40_coreform_identity(summary())
    row[_MEDIAL].update(
        {
            "sheet_generation": "medial-axis-hex-310",
            "topology_generation": "medial-axis-hex-309",
            "result_generation": "medial-axis-hex-308",
            "result_paired_sheet_ids": [[11, 21], [12, 22]],
            "result_local_thickness_m": [0.014, -0.010],
            "result_decomposition_cells": {"cell:1": [11, 21, 31]},
            "result_shared_topology_faces": [99],
            "result_interval_counts": {"curve:41": 7, "curve:42": 8},
            "result_interval_parity": "incompatible_odd",
            "result_minimum_scaled_jacobian": -0.08,
            "result_block_owner": "block:old",
            "accepted_medial_axis_export_sha256": "c" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "medial_axis_hexes_use_current_sheet_pairs_thickness_cells_topology_intervals_quality_block_and_export"
    ]


def test_v40_public_curve_chain_mismatch() -> None:
    row = _with_v40_coreform_identity(summary())
    row[_CHAIN].update(
        {
            "chain_generation": "curve-chain-hex-310",
            "corner_generation": "curve-chain-hex-309",
            "result_generation": "curve-chain-hex-308",
            "result_curve_chain_order": [101, 103, 102, 104],
            "result_chain_orientation": "clockwise",
            "result_interval_counts": [8, 9, 8, 10],
            "result_bias_directions": ["reverse"] * 4,
            "result_corner_interval_sums": [17, 18, 17, 18],
            "result_boundary_layer_thickness_m": [5.0e-4, -7.5e-4],
            "result_total_boundary_layer_thickness_m": 4.0e-3,
            "result_element_orientation": "inward_negative",
            "result_sideset_owner": "sideset:old",
            "accepted_curve_chain_export_sha256": "d" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "curve_chain_hexes_use_current_order_intervals_bias_corner_parity_boundary_layers_orientation_sideset_and_export"
    ]


def test_v40_source_smoothing_mismatch() -> None:
    row = _with_v40_coreform_identity(summary())
    row[_SMOOTH].update(
        {
            "algorithm_generation": "smooth-replay-310",
            "quality_generation": "smooth-replay-309",
            "result_generation": "smooth-replay-308",
            "replay_algorithm": "laplacian_unconstrained",
            "replay_iteration_count": 6,
            "replay_fixed_node_ids": [1, 2],
            "replay_fixed_node_displacement_m": [1.0e-3, 0.0],
            "replay_moved_node_ids": [103, 102, 101],
            "replay_node_displacement_m": [8.0e-4],
            "replay_maximum_allowed_displacement_m": 2.0e-4,
            "replay_quality_history": [0.18, 0.12, 0.10],
            "replay_checkpoint_generation_id": 310,
            "replay_database_owner": "gui:old",
            "replay_database_sha256": "e" * 64,
            "accepted_smoothing_result_sha256": "f" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "smoothing_replays_use_current_algorithm_iterations_constraints_motion_quality_checkpoint_database_and_result"
    ]


def test_v40_source_named_entity_mismatch() -> None:
    row = _with_v40_coreform_identity(summary())
    row[_NAMED].update(
        {
            "metadata_generation": "named-entity-replay-310",
            "open_generation": "named-entity-replay-309",
            "result_generation": "named-entity-replay-308",
            "replay_entity_names": {"volume:1": "stator"},
            "replay_metadata_attributes": {"volume:1": {"material": "air"}},
            "replay_group_membership": {"group:motor": ["volume:2"]},
            "replay_block_membership": {"block:99": ["volume:1"]},
            "replay_sideset_membership": {},
            "replay_save_generation_id": 310,
            "replay_open_generation_id": 309,
            "replay_export_owner": "gui:old",
            "replay_database_sha256": "0" * 64,
            "accepted_named_entity_export_sha256": "1" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "named_entity_replays_use_current_names_metadata_groups_sets_save_open_owner_database_and_export"
    ]


def test_v40_rejects_self_consistent_reused_medial_sheet() -> None:
    row = _with_v40_coreform_identity(summary())
    for key in ("paired_sheet_ids", "result_paired_sheet_ids"):
        row[_MEDIAL][key] = [[11, 12], [11, 22]]
    assert _public_result(row)["status"] == "needs_attention"


def test_v40_rejects_self_consistent_odd_corner_parity() -> None:
    row = _with_v40_coreform_identity(summary())
    for key in ("corner_interval_sums", "result_corner_interval_sums"):
        row[_CHAIN][key] = [17, 18, 17, 18]
    assert _public_result(row)["status"] == "needs_attention"


def test_v40_rejects_self_consistent_fixed_node_motion() -> None:
    row = _with_v40_coreform_identity(summary())
    for key in ("fixed_node_displacement_m", "replay_fixed_node_displacement_m"):
        row[_SMOOTH][key] = [0.0, 0.0, 1.0e-4, 0.0]
    assert _source_result(row)["status"] == "needs_attention"


def test_v40_rejects_self_consistent_named_set_owner_type() -> None:
    row = _with_v40_coreform_identity(summary())
    for key in ("sideset_membership", "replay_sideset_membership"):
        row[_NAMED][key] = {"sideset:20": ["volume:1"]}
    assert _source_result(row)["status"] == "needs_attention"
