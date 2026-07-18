from __future__ import annotations

from test_coreform_generalization_v36 import (
    _public_result,
    _source_result,
    _with_v36_sweep_transition_journal_exodus_identity,
    summary,
)


_PROMOTED_CASE_IDS = (
    "v37_public_periodic_hex_pair_transform_node_order_face_orientation_jacobian_owner_mismatch",
    "v37_public_high_order_hex_curve_midnode_projection_volume_jacobian_cad_owner_mismatch",
    "v37_source_imprint_merge_tolerance_topology_count_entity_lineage_checkpoint_mismatch",
    "v37_source_headless_batch_exit_license_fallback_journal_log_database_owner_mismatch",
)


def _with_v37_periodic_curved_imprint_batch_identity(row: dict) -> dict:
    row = _with_v36_sweep_transition_journal_exodus_identity(row)
    generation = "periodic-hex-contract-241"
    row[
        "periodic_hex_pair_transform_node_order_face_orientation_jacobian_region_export_mesh_result_generation_identity"
    ] = {
        "periodic_generation": generation,
        **{
            key: generation
            for key in (
                "pair_generation",
                "transform_generation",
                "node_order_generation",
                "orientation_generation",
                "jacobian_generation",
                "region_generation",
                "export_generation",
                "mesh_generation",
                "result_generation",
            )
        },
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
        "face_orientation": "source_outward_target_inward",
        "result_face_orientation": "source_outward_target_inward",
        "minimum_scaled_jacobian": 0.38,
        "result_minimum_scaled_jacobian": 0.38,
        "minimum_allowed_scaled_jacobian": 0.2,
        "result_minimum_allowed_scaled_jacobian": 0.2,
        "region_owner": "volume:21/periodic-pair-241",
        "result_region_owner": "volume:21/periodic-pair-241",
        "export_owner": "headless:periodic-export-241",
        "result_export_owner": "headless:periodic-export-241",
        "periodic_mesh_sha256": "1" * 64,
        "accepted_periodic_mesh_sha256": "1" * 64,
        "periodic_result_sha256": "2" * 64,
        "accepted_periodic_result_sha256": "2" * 64,
    }

    generation = "curved-hex-contract-241"
    row[
        "high_order_hex_curve_midnode_projection_volume_jacobian_order_cad_mesh_export_generation_identity"
    ] = {
        "curve_generation": generation,
        **{
            key: generation
            for key in (
                "edge_generation",
                "face_generation",
                "projection_generation",
                "volume_generation",
                "jacobian_generation",
                "order_generation",
                "cad_generation",
                "mesh_generation",
                "export_generation",
                "result_generation",
            )
        },
        "element_family": "HEX27",
        "result_element_family": "HEX27",
        "polynomial_order": 2,
        "result_polynomial_order": 2,
        "edge_midnode_roles": ["edge_mid"] * 12,
        "result_edge_midnode_roles": ["edge_mid"] * 12,
        "face_midnode_roles": ["face_mid"] * 6,
        "result_face_midnode_roles": ["face_mid"] * 6,
        "maximum_cad_projection_distance_m": 2.0e-8,
        "result_maximum_cad_projection_distance_m": 2.0e-8,
        "allowed_cad_projection_distance_m": 1.0e-6,
        "result_allowed_cad_projection_distance_m": 1.0e-6,
        "cad_volume_m3": 1.0,
        "curved_mesh_volume_m3": 1.0,
        "result_curved_mesh_volume_m3": 1.0,
        "volume_tolerance_m3": 1.0e-9,
        "result_volume_tolerance_m3": 1.0e-9,
        "minimum_high_order_jacobian": 0.31,
        "result_minimum_high_order_jacobian": 0.31,
        "minimum_allowed_high_order_jacobian": 0.1,
        "result_minimum_allowed_high_order_jacobian": 0.1,
        "cad_owner": "cad:volume21/curve-241",
        "result_cad_owner": "cad:volume21/curve-241",
        "mesh_owner": "headless:hex27-241",
        "result_mesh_owner": "headless:hex27-241",
        "curved_export_sha256": "3" * 64,
        "accepted_curved_export_sha256": "3" * 64,
    }

    generation = "imprint-merge-contract-241"
    row[
        "imprint_merge_tolerance_topology_count_entity_lineage_command_checkpoint_model_final_database_result_generation_identity"
    ] = {
        "imprint_merge_generation": generation,
        **{
            key: generation
            for key in (
                "tolerance_generation",
                "topology_generation",
                "lineage_generation",
                "command_generation",
                "checkpoint_generation",
                "model_generation",
                "final_generation",
                "database_generation",
                "result_generation",
            )
        },
        "merge_tolerance_m": 1.0e-7,
        "replay_merge_tolerance_m": 1.0e-7,
        "coincident_topology_count": 4,
        "replay_coincident_topology_count": 4,
        "merged_entity_lineage": [[11, 21], [12, 22], [13, 23], [14, 24]],
        "replay_merged_entity_lineage": [[11, 21], [12, 22], [13, 23], [14, 24]],
        "command_sequence": ["imprint all", "merge all", "compress all"],
        "replay_command_sequence": ["imprint all", "merge all", "compress all"],
        "checkpoint_owner": "headless:checkpoint-241",
        "replay_checkpoint_owner": "headless:checkpoint-241",
        "model_generation_id": 41,
        "replay_model_generation_id": 41,
        "final_topology_counts": {"volume": 2, "surface": 11, "curve": 20},
        "replay_final_topology_counts": {"volume": 2, "surface": 11, "curve": 20},
        "database_sha256": "4" * 64,
        "replay_database_sha256": "4" * 64,
        "imprint_result_sha256": "5" * 64,
        "accepted_imprint_result_sha256": "5" * 64,
    }

    generation = "headless-batch-contract-241"
    row[
        "headless_batch_exit_license_fallback_journal_log_database_command_process_result_generation_identity"
    ] = {
        "batch_generation": generation,
        **{
            key: generation
            for key in (
                "process_generation",
                "license_generation",
                "journal_generation",
                "log_generation",
                "database_generation",
                "command_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "execution_mode": "nographics_batch",
        "result_execution_mode": "nographics_batch",
        "gui_launched": False,
        "result_gui_launched": False,
        "process_exit_code": 0,
        "result_process_exit_code": 0,
        "license_warning": "License Error: No license found",
        "result_license_warning": "License Error: No license found",
        "license_fallback": "limited_mode_batch_completed",
        "result_license_fallback": "limited_mode_batch_completed",
        "journal_completion_marker": "CAEAI_BATCH_COMPLETE",
        "result_journal_completion_marker": "CAEAI_BATCH_COMPLETE",
        "log_owner": "headless:batch-log-241",
        "result_log_owner": "headless:batch-log-241",
        "database_save_generation_id": 41,
        "result_database_save_generation_id": 41,
        "command_sha256": "6" * 64,
        "result_command_sha256": "6" * 64,
        "process_owner": "coreform_cubit:-nographics:-batch",
        "result_process_owner": "coreform_cubit:-nographics:-batch",
        "batch_result_sha256": "7" * 64,
        "accepted_batch_result_sha256": "7" * 64,
    }
    return row


def test_v37_positive_public_and_source_contracts() -> None:
    row = _with_v37_periodic_curved_imprint_batch_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"


def test_v37_public_periodic_hex_pair_transform_node_order_face_orientation_jacobian_owner_mismatch() -> None:
    row = _with_v37_periodic_curved_imprint_batch_identity(summary())
    row[
        "periodic_hex_pair_transform_node_order_face_orientation_jacobian_region_export_mesh_result_generation_identity"
    ].update(
        {
            "transform_generation": "periodic-hex-contract-240",
            "result_target_surface": "surface:13",
            "result_rigid_transform": [[1.0, 0.0, 0.0, -1.0]],
            "result_periodic_node_pairs": [[1, 104], [2, 103]],
            "result_face_orientation": "same_normal",
            "result_minimum_scaled_jacobian": -0.2,
            "result_region_owner": "volume:stale",
            "result_export_owner": "gui:interactive",
            "accepted_periodic_mesh_sha256": "8" * 64,
            "accepted_periodic_result_sha256": "9" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "periodic_hex_pairs_use_current_transform_node_order_orientation_jacobian_region_export_mesh_and_result"
    ]


def test_v37_public_high_order_hex_curve_midnode_projection_volume_jacobian_cad_owner_mismatch() -> None:
    row = _with_v37_periodic_curved_imprint_batch_identity(summary())
    row[
        "high_order_hex_curve_midnode_projection_volume_jacobian_order_cad_mesh_export_generation_identity"
    ].update(
        {
            "edge_generation": "curved-hex-contract-240",
            "result_element_family": "HEX8",
            "result_polynomial_order": 1,
            "result_edge_midnode_roles": [],
            "result_face_midnode_roles": ["corner"] * 6,
            "result_maximum_cad_projection_distance_m": 1.0e-2,
            "result_curved_mesh_volume_m3": 0.8,
            "result_minimum_high_order_jacobian": -0.1,
            "result_cad_owner": "faceted:stale",
            "result_mesh_owner": "gui:stale",
            "accepted_curved_export_sha256": "a" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "high_order_curved_hexes_use_current_midnodes_projection_volume_jacobian_order_cad_mesh_and_export"
    ]


def test_v37_source_imprint_merge_tolerance_topology_count_entity_lineage_checkpoint_mismatch() -> None:
    row = _with_v37_periodic_curved_imprint_batch_identity(summary())
    row[
        "imprint_merge_tolerance_topology_count_entity_lineage_command_checkpoint_model_final_database_result_generation_identity"
    ].update(
        {
            "tolerance_generation": "imprint-merge-contract-240",
            "replay_merge_tolerance_m": 1.0e-2,
            "replay_coincident_topology_count": 0,
            "replay_merged_entity_lineage": [[11, 99]],
            "replay_command_sequence": ["merge all", "imprint all"],
            "replay_checkpoint_owner": "gui:old-checkpoint",
            "replay_model_generation_id": 40,
            "replay_final_topology_counts": {"volume": 3, "surface": 15},
            "replay_database_sha256": "b" * 64,
            "accepted_imprint_result_sha256": "c" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "imprint_merge_replays_use_current_tolerance_topology_lineage_order_checkpoint_model_database_and_result"
    ]


def test_v37_source_headless_batch_exit_license_fallback_journal_log_database_owner_mismatch() -> None:
    row = _with_v37_periodic_curved_imprint_batch_identity(summary())
    row[
        "headless_batch_exit_license_fallback_journal_log_database_command_process_result_generation_identity"
    ].update(
        {
            "process_generation": "headless-batch-contract-240",
            "result_execution_mode": "gui",
            "result_gui_launched": True,
            "result_process_exit_code": 1,
            "result_license_fallback": "fatal",
            "result_journal_completion_marker": "INCOMPLETE",
            "result_log_owner": "interactive:unknown",
            "result_database_save_generation_id": 40,
            "result_command_sha256": "d" * 64,
            "result_process_owner": "cubit_gui",
            "accepted_batch_result_sha256": "e" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "headless_batches_use_current_exit_license_fallback_journal_log_database_command_process_and_result"
    ]


def test_v37_rejects_self_consistent_reflection_as_periodic_transform() -> None:
    row = _with_v37_periodic_curved_imprint_batch_identity(summary())
    identity = row[
        "periodic_hex_pair_transform_node_order_face_orientation_jacobian_region_export_mesh_result_generation_identity"
    ]
    reflection = [
        [-1.0, 0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    identity["rigid_transform"] = reflection
    identity["result_rigid_transform"] = reflection
    assert _public_result(row)["status"] == "needs_attention"


def test_v37_accepts_known_license_warning_only_with_headless_completion() -> None:
    row = _with_v37_periodic_curved_imprint_batch_identity(summary())
    assert _source_result(row)["status"] == "ok"
    identity = row[
        "headless_batch_exit_license_fallback_journal_log_database_command_process_result_generation_identity"
    ]
    identity["journal_completion_marker"] = "INCOMPLETE"
    identity["result_journal_completion_marker"] = "INCOMPLETE"
    assert _source_result(row)["status"] == "needs_attention"
