from __future__ import annotations

from test_coreform_generalization_v37 import (
    _public_result,
    _source_result,
    _with_v37_periodic_curved_imprint_batch_identity,
    summary,
)

_PROMOTED_CASE_IDS = (
    "v38_public_midsurface_shell_thickness_normal_area_volume_reconstruction_block_export_mismatch",
    "v38_public_cohesive_crack_interface_node_pair_normal_orientation_traction_block_mismatch",
    "v38_source_virtual_geometry_small_entity_suppression_topology_map_quality_undo_checkpoint_mismatch",
    "v38_source_anisotropic_crack_metric_eigenframe_sizing_alignment_hybrid_region_export_mismatch",
)


def _with_v38_shell_cohesive_virtual_anisotropic_identity(row: dict) -> dict:
    row = _with_v37_periodic_curved_imprint_batch_identity(row)

    generation = "midsurface-shell-258"
    row[
        "midsurface_shell_facepair_thickness_normal_area_volume_block_sideset_geometry_export_result_generation_identity"
    ] = {
        "midsurface_generation": generation,
        **{
            key: generation
            for key in (
                "facepair_generation",
                "thickness_generation",
                "normal_generation",
                "area_generation",
                "volume_generation",
                "block_generation",
                "geometry_generation",
                "export_generation",
                "result_generation",
            )
        },
        "paired_face_ids": [[11, 12], [13, 14]],
        "result_paired_face_ids": [[11, 12], [13, 14]],
        "paired_face_distance_m": [2.0e-3, 3.0e-3],
        "result_paired_face_distance_m": [2.0e-3, 3.0e-3],
        "shell_thickness_m": [2.0e-3, 3.0e-3],
        "result_shell_thickness_m": [2.0e-3, 3.0e-3],
        "paired_normal_dot": [-1.0, -1.0],
        "result_paired_normal_dot": [-1.0, -1.0],
        "midsurface_area_m2": [1.0, 2.0],
        "result_midsurface_area_m2": [1.0, 2.0],
        "reconstructed_volume_m3": 8.0e-3,
        "result_reconstructed_volume_m3": 8.0e-3,
        "source_volume_m3": 8.0e-3,
        "volume_tolerance_m3": 1.0e-9,
        "shell_block": "block:shell-10",
        "result_shell_block": "block:shell-10",
        "shell_sidesets": ["sideset:top", "sideset:bottom"],
        "result_shell_sidesets": ["sideset:top", "sideset:bottom"],
        "geometry_owner": "headless:midsurface-258",
        "result_geometry_owner": "headless:midsurface-258",
        "midsurface_export_sha256": "1" * 64,
        "accepted_midsurface_export_sha256": "1" * 64,
        "midsurface_result_sha256": "2" * 64,
        "accepted_midsurface_result_sha256": "2" * 64,
    }

    generation = "cohesive-crack-258"
    row[
        "cohesive_crack_face_nodepair_front_normal_orientation_traction_block_jacobian_mesh_export_result_generation_identity"
    ] = {
        "cohesive_generation": generation,
        **{
            key: generation
            for key in (
                "face_generation",
                "nodepair_generation",
                "front_generation",
                "normal_generation",
                "orientation_generation",
                "traction_generation",
                "block_generation",
                "jacobian_generation",
                "mesh_generation",
                "export_generation",
                "result_generation",
            )
        },
        "duplicated_face_pairs": [[21, 121], [22, 122]],
        "result_duplicated_face_pairs": [[21, 121], [22, 122]],
        "duplicated_node_pairs": [[1, 101], [2, 102], [3, 103], [4, 104]],
        "result_duplicated_node_pairs": [[1, 101], [2, 102], [3, 103], [4, 104]],
        "crack_front_nodes": [1, 2],
        "result_crack_front_nodes": [1, 2],
        "interface_normal_relation": "opposed_outward_normals",
        "result_interface_normal_relation": "opposed_outward_normals",
        "cohesive_orientation": "positive_reference_orientation",
        "result_cohesive_orientation": "positive_reference_orientation",
        "traction_direction": "source_to_target",
        "result_traction_direction": "source_to_target",
        "cohesive_block": "block:cohesive-20",
        "result_cohesive_block": "block:cohesive-20",
        "minimum_scaled_jacobian": 0.25,
        "result_minimum_scaled_jacobian": 0.25,
        "minimum_allowed_scaled_jacobian": 0.1,
        "result_minimum_allowed_scaled_jacobian": 0.1,
        "mesh_owner": "headless:cohesive-258",
        "result_mesh_owner": "headless:cohesive-258",
        "cohesive_export_sha256": "3" * 64,
        "accepted_cohesive_export_sha256": "3" * 64,
        "cohesive_result_sha256": "4" * 64,
        "accepted_cohesive_result_sha256": "4" * 64,
    }

    generation = "virtual-geometry-258"
    row[
        "virtual_geometry_suppression_topology_map_inheritance_quality_undo_checkpoint_database_result_generation_identity"
    ] = {
        "virtual_generation": generation,
        **{
            key: generation
            for key in (
                "suppression_generation",
                "topology_generation",
                "inheritance_generation",
                "quality_generation",
                "undo_generation",
                "checkpoint_generation",
                "database_generation",
                "result_generation",
            )
        },
        "suppressed_curve_ids": [7, 8],
        "replay_suppressed_curve_ids": [7, 8],
        "suppressed_surface_ids": [5],
        "replay_suppressed_surface_ids": [5],
        "virtual_topology_map": {"7": 17, "8": 18, "5": 15},
        "replay_virtual_topology_map": {"7": 17, "8": 18, "5": 15},
        "block_inheritance": {"volume:1": "block:10"},
        "replay_block_inheritance": {"volume:1": "block:10"},
        "sideset_inheritance": {"surface:15": "sideset:20"},
        "replay_sideset_inheritance": {"surface:15": "sideset:20"},
        "minimum_quality_before": 0.15,
        "replay_minimum_quality_before": 0.15,
        "minimum_quality_after": 0.30,
        "replay_minimum_quality_after": 0.30,
        "undo_restored_topology": True,
        "replay_undo_restored_topology": True,
        "checkpoint_owner": "headless:virtual-checkpoint-258",
        "replay_checkpoint_owner": "headless:virtual-checkpoint-258",
        "database_sha256": "5" * 64,
        "replay_database_sha256": "5" * 64,
        "virtual_result_sha256": "6" * 64,
        "accepted_virtual_result_sha256": "6" * 64,
    }

    generation = "anisotropic-crack-258"
    row[
        "anisotropic_crack_metric_eigenframe_sizing_alignment_cohesive_hexconformity_quality_region_export_result_generation_identity"
    ] = {
        "anisotropic_generation": generation,
        **{
            key: generation
            for key in (
                "metric_generation",
                "eigenframe_generation",
                "sizing_generation",
                "alignment_generation",
                "cohesive_generation",
                "conformity_generation",
                "quality_generation",
                "region_generation",
                "export_generation",
                "result_generation",
            )
        },
        "metric_eigenvalues": [1.0, 4.0, 16.0],
        "result_metric_eigenvalues": [1.0, 4.0, 16.0],
        "metric_eigenvectors": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "result_metric_eigenvectors": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "sizing_field_m": [1.0e-3, 5.0e-4, 2.5e-4],
        "result_sizing_field_m": [1.0e-3, 5.0e-4, 2.5e-4],
        "crack_front_tangent": [1.0, 0.0, 0.0],
        "result_crack_front_tangent": [1.0, 0.0, 0.0],
        "aligned_metric_axis": 0,
        "result_aligned_metric_axis": 0,
        "cohesive_transition_face_count": 12,
        "result_cohesive_transition_face_count": 12,
        "hex_interface_conformal": True,
        "result_hex_interface_conformal": True,
        "minimum_quality": 0.24,
        "result_minimum_quality": 0.24,
        "minimum_allowed_quality": 0.1,
        "result_minimum_allowed_quality": 0.1,
        "region_owner": "volume:crack-zone-258",
        "result_region_owner": "volume:crack-zone-258",
        "anisotropic_export_sha256": "7" * 64,
        "accepted_anisotropic_export_sha256": "7" * 64,
        "anisotropic_result_sha256": "8" * 64,
        "accepted_anisotropic_result_sha256": "8" * 64,
    }
    return row


def test_v38_positive_public_and_source_contracts() -> None:
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"


def test_v38_public_midsurface_shell_thickness_normal_area_volume_reconstruction_block_export_mismatch() -> (
    None
):
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(summary())
    row[
        "midsurface_shell_facepair_thickness_normal_area_volume_block_sideset_geometry_export_result_generation_identity"
    ].update(
        {
            "thickness_generation": "midsurface-shell-257",
            "volume_generation": "midsurface-shell-256",
            "result_generation": "midsurface-shell-255",
            "result_paired_face_ids": [[11, 13]],
            "result_shell_thickness_m": [2.0e-2, 3.0e-2],
            "result_paired_normal_dot": [1.0, 1.0],
            "result_midsurface_area_m2": [2.0, 1.0],
            "result_reconstructed_volume_m3": 8.0e-2,
            "result_shell_block": "block:old",
            "result_shell_sidesets": ["sideset:old"],
            "result_geometry_owner": "gui:old",
            "accepted_midsurface_result_sha256": "9" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "midsurface_shells_use_current_facepairs_thickness_normals_area_volume_sets_geometry_and_result"
    ]


def test_v38_public_cohesive_crack_interface_node_pair_normal_orientation_traction_block_mismatch() -> (
    None
):
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(summary())
    row[
        "cohesive_crack_face_nodepair_front_normal_orientation_traction_block_jacobian_mesh_export_result_generation_identity"
    ].update(
        {
            "nodepair_generation": "cohesive-crack-257",
            "orientation_generation": "cohesive-crack-256",
            "result_generation": "cohesive-crack-255",
            "result_duplicated_face_pairs": [[21, 22]],
            "result_duplicated_node_pairs": [[1, 104], [2, 103]],
            "result_crack_front_nodes": [3, 4],
            "result_interface_normal_relation": "same_normal",
            "result_cohesive_orientation": "negative",
            "result_traction_direction": "target_to_source",
            "result_cohesive_block": "block:bulk",
            "result_minimum_scaled_jacobian": -0.2,
            "result_mesh_owner": "gui:old",
            "accepted_cohesive_result_sha256": "a" * 64,
        }
    )
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cohesive_cracks_use_current_faces_nodes_front_normals_orientation_traction_quality_mesh_and_result"
    ]


def test_v38_source_virtual_geometry_small_entity_suppression_topology_map_quality_undo_checkpoint_mismatch() -> (
    None
):
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(summary())
    row[
        "virtual_geometry_suppression_topology_map_inheritance_quality_undo_checkpoint_database_result_generation_identity"
    ].update(
        {
            "topology_generation": "virtual-geometry-257",
            "undo_generation": "virtual-geometry-256",
            "result_generation": "virtual-geometry-255",
            "replay_suppressed_curve_ids": [8, 9],
            "replay_suppressed_surface_ids": [],
            "replay_virtual_topology_map": {"8": 99},
            "replay_block_inheritance": {"volume:1": "block:old"},
            "replay_sideset_inheritance": {},
            "replay_minimum_quality_after": 0.05,
            "replay_undo_restored_topology": False,
            "replay_checkpoint_owner": "gui:old",
            "accepted_virtual_result_sha256": "b" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "virtual_geometry_replays_use_current_suppression_topology_inheritance_quality_undo_checkpoint_database_and_result"
    ]


def test_v38_source_anisotropic_crack_metric_eigenframe_sizing_alignment_hybrid_region_export_mismatch() -> (
    None
):
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(summary())
    row[
        "anisotropic_crack_metric_eigenframe_sizing_alignment_cohesive_hexconformity_quality_region_export_result_generation_identity"
    ].update(
        {
            "metric_generation": "anisotropic-crack-257",
            "conformity_generation": "anisotropic-crack-256",
            "result_generation": "anisotropic-crack-255",
            "result_metric_eigenvalues": [16.0, 4.0, -1.0],
            "result_metric_eigenvectors": [[0.0, 1.0, 0.0]],
            "result_sizing_field_m": [-1.0, 5.0e-4],
            "result_crack_front_tangent": [0.0, 1.0, 0.0],
            "result_aligned_metric_axis": 2,
            "result_cohesive_transition_face_count": 8,
            "result_hex_interface_conformal": False,
            "result_minimum_quality": -0.1,
            "result_region_owner": "volume:old",
            "accepted_anisotropic_result_sha256": "c" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "anisotropic_crack_regions_use_spd_metric_eigenframe_sizing_alignment_cohesive_hex_quality_and_result"
    ]


def test_v38_rejects_self_consistent_same_normal_shell_pair() -> None:
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(summary())
    identity = row[
        "midsurface_shell_facepair_thickness_normal_area_volume_block_sideset_geometry_export_result_generation_identity"
    ]
    identity["paired_normal_dot"] = [1.0, 1.0]
    identity["result_paired_normal_dot"] = [1.0, 1.0]
    assert _public_result(row)["status"] == "needs_attention"


def test_v38_rejects_self_consistent_same_normal_cohesive_interface() -> None:
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(summary())
    identity = row[
        "cohesive_crack_face_nodepair_front_normal_orientation_traction_block_jacobian_mesh_export_result_generation_identity"
    ]
    identity["interface_normal_relation"] = "same_normal"
    identity["result_interface_normal_relation"] = "same_normal"
    assert _public_result(row)["status"] == "needs_attention"


def test_v38_rejects_self_consistent_noninjective_virtual_topology_map() -> None:
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(summary())
    identity = row[
        "virtual_geometry_suppression_topology_map_inheritance_quality_undo_checkpoint_database_result_generation_identity"
    ]
    identity["virtual_topology_map"] = {"7": 17, "8": 17, "5": 15}
    identity["replay_virtual_topology_map"] = {"7": 17, "8": 17, "5": 15}
    assert _source_result(row)["status"] == "needs_attention"


def test_v38_rejects_self_consistent_metric_size_law_mismatch() -> None:
    row = _with_v38_shell_cohesive_virtual_anisotropic_identity(summary())
    identity = row[
        "anisotropic_crack_metric_eigenframe_sizing_alignment_cohesive_hexconformity_quality_region_export_result_generation_identity"
    ]
    identity["sizing_field_m"] = [1.0e-3, 8.0e-4, 2.5e-4]
    identity["result_sizing_field_m"] = [1.0e-3, 8.0e-4, 2.5e-4]
    assert _source_result(row)["status"] == "needs_attention"
