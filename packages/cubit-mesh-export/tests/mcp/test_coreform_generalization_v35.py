from __future__ import annotations

from test_coreform_generalization_v34 import (
    _public_result,
    _source_result,
    _with_v34_high_order_sheet_sweep_checkpoint_identity,
    summary,
)


_PROMOTED_CASE_IDS = (
    "v35_public_anisotropic_hex_metric_principal_direction_gradation_jacobian_error_mismatch",
    "v35_public_curved_highorder_boundary_normal_hausdorff_metric_volume_jacobian_mismatch",
    "v35_source_sideset_skin_remesh_adjacent_block_normal_entity_generation_mismatch",
    "v35_source_parallel_sculpt_partition_seed_stitch_ghost_qa_determinism_mismatch",
)


def _with_v35_metric_curved_sideset_sculpt_identity(row: dict) -> dict:
    row = _with_v34_high_order_sheet_sweep_checkpoint_identity(row)
    generation = "anisotropic-hex-221"
    row[
        "anisotropic_hex_metric_direction_size_gradation_alignment_jacobian_block_mesh_result_generation_identity"
    ] = {
        "anisotropic_generation": generation,
        **{key: generation for key in ("metric_generation", "direction_generation", "size_generation", "gradation_generation", "quality_generation", "block_generation", "mesh_generation", "result_generation")},
        "metric_eigenvectors": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "result_metric_eigenvectors": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "principal_sizes_m": [0.001, 0.002, 0.004],
        "result_principal_sizes_m": [0.001, 0.002, 0.004],
        "gradation_ratio": 1.3,
        "result_gradation_ratio": 1.3,
        "maximum_gradation_ratio": 1.5,
        "result_maximum_gradation_ratio": 1.5,
        "alignment_error_deg": 2.0,
        "result_alignment_error_deg": 2.0,
        "maximum_alignment_error_deg": 5.0,
        "result_maximum_alignment_error_deg": 5.0,
        "minimum_scaled_jacobian": 0.4,
        "result_minimum_scaled_jacobian": 0.4,
        "minimum_allowed_scaled_jacobian": 0.2,
        "result_minimum_allowed_scaled_jacobian": 0.2,
        "block_ids": [101, 102],
        "result_block_ids": [101, 102],
        "anisotropic_mesh_sha256": "1" * 64,
        "accepted_anisotropic_mesh_sha256": "1" * 64,
    }
    generation = "curved-boundary-221"
    row[
        "curved_highorder_boundary_normal_hausdorff_area_volume_jacobian_order_geometry_result_generation_identity"
    ] = {
        "curved_generation": generation,
        **{key: generation for key in ("normal_generation", "hausdorff_generation", "area_generation", "volume_generation", "jacobian_generation", "order_generation", "geometry_generation", "result_generation")},
        "minimum_normal_dot": 0.999,
        "result_minimum_normal_dot": 0.999,
        "hausdorff_error_m": 1.0e-5,
        "result_hausdorff_error_m": 1.0e-5,
        "maximum_hausdorff_error_m": 5.0e-5,
        "result_maximum_hausdorff_error_m": 5.0e-5,
        "cad_surface_area_m2": 1.0,
        "mesh_surface_area_m2": 0.99999,
        "surface_measure_tolerance": 1.0e-4,
        "cad_volume_m3": 0.25,
        "mesh_volume_m3": 0.24999,
        "volume_measure_tolerance": 1.0e-4,
        "minimum_curved_jacobian": 0.35,
        "result_minimum_curved_jacobian": 0.35,
        "polynomial_order": 2,
        "result_polynomial_order": 2,
        "geometry_owner": "cad:volume11/surface21",
        "result_geometry_owner": "cad:volume11/surface21",
        "curved_geometry_sha256": "2" * 64,
        "accepted_curved_geometry_sha256": "2" * 64,
    }
    generation = "sideset-remesh-221"
    row[
        "sideset_skin_remesh_adjacent_block_normal_face_multiplicity_entity_journal_result_generation_identity"
    ] = {
        "sideset_generation": generation,
        **{key: generation for key in ("skin_generation", "remesh_generation", "adjacency_generation", "normal_generation", "face_generation", "entity_generation", "journal_generation", "result_generation")},
        "sideset_id": 201,
        "restored_sideset_id": 201,
        "adjacent_block_ids": [101],
        "restored_adjacent_block_ids": [101],
        "outward_normals": [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        "restored_outward_normals": [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        "face_ids": [31, 32],
        "restored_face_ids": [31, 32],
        "face_multiplicities": [1, 1],
        "restored_face_multiplicities": [1, 1],
        "source_entity_keys": ["surface:left:0", "surface:left:1"],
        "restored_source_entity_keys": ["surface:left:0", "surface:left:1"],
        "remesh_revision": 7,
        "restored_remesh_revision": 7,
        "sideset_journal_sha256": "3" * 64,
        "replayed_sideset_journal_sha256": "3" * 64,
        "sideset_result_sha256": "4" * 64,
        "accepted_sideset_result_sha256": "4" * 64,
    }
    generation = "parallel-sculpt-221"
    row[
        "parallel_sculpt_seed_rank_owned_ghost_stitch_qa_connectivity_invocation_export_generation_identity"
    ] = {
        "sculpt_generation": generation,
        **{key: generation for key in ("seed_generation", "rank_generation", "partition_generation", "stitch_generation", "qa_generation", "connectivity_generation", "invocation_generation", "result_generation")},
        "partition_seed": 12345,
        "replay_partition_seed": 12345,
        "rank_count": 4,
        "replay_rank_count": 4,
        "owned_cell_counts": [100, 100, 100, 100],
        "replay_owned_cell_counts": [100, 100, 100, 100],
        "ghost_cell_counts": [10, 12, 12, 10],
        "replay_ghost_cell_counts": [10, 12, 12, 10],
        "stitched_interface_pair_count": 36,
        "replay_stitched_interface_pair_count": 36,
        "qa_record": ["Cubit", "2026.6", "parallel_sculpt"],
        "replay_qa_record": ["Cubit", "2026.6", "parallel_sculpt"],
        "invocation_owner": "headless:batch42",
        "replay_invocation_owner": "headless:batch42",
        "connectivity_sha256": "5" * 64,
        "replay_connectivity_sha256": "5" * 64,
        "sculpt_export_sha256": "6" * 64,
        "accepted_sculpt_export_sha256": "6" * 64,
    }
    return row


def test_v35_positive_all_four_contracts() -> None:
    row = _with_v35_metric_curved_sideset_sculpt_identity(summary())
    assert _public_result(row)["status"] == "ok"
    assert _source_result(row)["status"] == "ok"


def test_v35_public_anisotropic_hex_metric_principal_direction_gradation_jacobian_error_mismatch() -> None:
    row = _with_v35_metric_curved_sideset_sculpt_identity(summary())
    row["anisotropic_hex_metric_direction_size_gradation_alignment_jacobian_block_mesh_result_generation_identity"].update({"metric_generation": "anisotropic-hex-220", "quality_generation": "anisotropic-hex-219", "result_generation": "anisotropic-hex-218", "result_metric_eigenvectors": [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], "result_principal_sizes_m": [0.004, 0.002, 0.001], "result_gradation_ratio": 2.0, "result_maximum_gradation_ratio": 1.1, "result_alignment_error_deg": 20.0, "result_maximum_alignment_error_deg": 1.0, "result_minimum_scaled_jacobian": -0.1, "result_minimum_allowed_scaled_jacobian": 0.5, "result_block_ids": [999], "accepted_anisotropic_mesh_sha256": "7" * 64})
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["anisotropic_hexes_use_current_metric_directions_sizes_gradation_alignment_jacobian_blocks_and_mesh"]


def test_v35_public_curved_highorder_boundary_normal_hausdorff_metric_volume_jacobian_mismatch() -> None:
    row = _with_v35_metric_curved_sideset_sculpt_identity(summary())
    row["curved_highorder_boundary_normal_hausdorff_area_volume_jacobian_order_geometry_result_generation_identity"].update({"normal_generation": "curved-boundary-220", "geometry_generation": "curved-boundary-219", "result_generation": "curved-boundary-218", "result_minimum_normal_dot": -0.5, "result_hausdorff_error_m": 0.01, "result_maximum_hausdorff_error_m": 1.0e-6, "mesh_surface_area_m2": 0.8, "mesh_volume_m3": 0.2, "result_minimum_curved_jacobian": -0.2, "result_polynomial_order": 1, "result_geometry_owner": "cad:stale", "accepted_curved_geometry_sha256": "8" * 64})
    result = _public_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["curved_high_order_boundaries_use_current_normals_hausdorff_measures_jacobian_order_geometry_and_result"]


def test_v35_source_sideset_skin_remesh_adjacent_block_normal_entity_generation_mismatch() -> None:
    row = _with_v35_metric_curved_sideset_sculpt_identity(summary())
    row["sideset_skin_remesh_adjacent_block_normal_face_multiplicity_entity_journal_result_generation_identity"].update({"remesh_generation": "sideset-remesh-220", "entity_generation": "sideset-remesh-219", "result_generation": "sideset-remesh-218", "restored_sideset_id": 202, "restored_adjacent_block_ids": [102, 103], "restored_outward_normals": [[-1.0, 0.0, 0.0]], "restored_face_ids": [32, 31, 31], "restored_face_multiplicities": [2, 1], "restored_source_entity_keys": ["surface:right:0"], "restored_remesh_revision": 6, "replayed_sideset_journal_sha256": "9" * 64, "accepted_sideset_result_sha256": "a" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["skinned_sidesets_use_current_remesh_adjacency_normals_faces_entities_journal_and_result"]


def test_v35_source_parallel_sculpt_partition_seed_stitch_ghost_qa_determinism_mismatch() -> None:
    row = _with_v35_metric_curved_sideset_sculpt_identity(summary())
    row["parallel_sculpt_seed_rank_owned_ghost_stitch_qa_connectivity_invocation_export_generation_identity"].update({"seed_generation": "parallel-sculpt-220", "connectivity_generation": "parallel-sculpt-219", "result_generation": "parallel-sculpt-218", "replay_partition_seed": 54321, "replay_rank_count": 3, "replay_owned_cell_counts": [120, 120, 120], "replay_ghost_cell_counts": [0, 0, 0], "replay_stitched_interface_pair_count": 0, "replay_qa_record": ["Cubit", "2025.8", "serial"], "replay_invocation_owner": "gui:interactive", "replay_connectivity_sha256": "b" * 64, "accepted_sculpt_export_sha256": "c" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["parallel_sculpt_uses_current_seed_ranks_owned_ghost_stitch_qa_connectivity_invocation_and_export"]


def test_v35_rejects_self_consistent_nonorthogonal_metric() -> None:
    row = _with_v35_metric_curved_sideset_sculpt_identity(summary())
    identity = row["anisotropic_hex_metric_direction_size_gradation_alignment_jacobian_block_mesh_result_generation_identity"]
    identity["metric_eigenvectors"] = [[1.0, 0.0, 0.0]] * 3
    identity["result_metric_eigenvectors"] = [[1.0, 0.0, 0.0]] * 3
    assert _public_result(row)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_curved_measure_drift() -> None:
    row = _with_v35_metric_curved_sideset_sculpt_identity(summary())
    row["curved_highorder_boundary_normal_hausdorff_area_volume_jacobian_order_geometry_result_generation_identity"]["mesh_volume_m3"] = 0.2
    assert _public_result(row)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_internal_skin_face() -> None:
    row = _with_v35_metric_curved_sideset_sculpt_identity(summary())
    identity = row["sideset_skin_remesh_adjacent_block_normal_face_multiplicity_entity_journal_result_generation_identity"]
    identity["adjacent_block_ids"] = [101, 102]
    identity["restored_adjacent_block_ids"] = [101, 102]
    assert _source_result(row)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_gui_sculpt_invocation() -> None:
    row = _with_v35_metric_curved_sideset_sculpt_identity(summary())
    identity = row["parallel_sculpt_seed_rank_owned_ghost_stitch_qa_connectivity_invocation_export_generation_identity"]
    identity["invocation_owner"] = "gui:interactive"
    identity["replay_invocation_owner"] = "gui:interactive"
    assert _source_result(row)["status"] == "needs_attention"
