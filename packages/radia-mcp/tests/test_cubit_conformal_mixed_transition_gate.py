import copy
import json

import pytest

from radia_mcp.cubit.server import (
    cubit_conformal_hex_pyramid_tet_interface_gate,
    cubit_mixed_transition_source_gate,
)


def summary() -> dict:
    return {
        "version": "2025.12",
        "source_kind": "source_native_local_mixed_element_journal_path_adapted",
        "source_journal": "01_Tet_Hex_Pyramid.jou",
        "source_sha256": "a" * 64,
        "source_commands": [
            "brick x 2 y 1 z 1",
            "webcut volume 1 with plane xplane imprint merge",
            "volume 1 scheme map",
            "mesh volume 1",
            "volume 2 scheme tetmesh",
            "mesh volume 2",
            "block 1 add hex all",
            "block 2 add pyramid all",
            "block 3 add tet all",
            "block 4 add tri all",
            "block 5 add face all",
            "volume all scale 0.001",
        ],
        "execution_mode": "headless_combined_journal_then_python_inventory",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "element_counts": {"hex": 1, "pyramid": 1, "tet": 10, "wedge": 0},
        "per_volume_element_counts": {
            "1": {"hex": 1, "pyramid": 0, "tet": 0, "wedge": 0},
            "2": {"hex": 0, "pyramid": 1, "tet": 10, "wedge": 0},
        },
        "quality": {
            "hex": {"scaled_jacobian": {"count": 1, "min": 1.0}},
            "pyramid": {"scaled_jacobian": {"count": 1, "min": 1.0}},
            "tet": {"scaled_jacobian": {"count": 10, "min": 0.3826}},
        },
        "mesh_identity": {
            "generation": "mesh-generation-42",
            "sha256": "d" * 64,
        },
        "quality_report_identity": {
            "mesh_generation": "mesh-generation-42",
            "mesh_sha256": "d" * 64,
            "report_sha256": "e" * 64,
        },
        "quality_scope_identity": {
            "mesh_volume_ids": ["1", "2"],
            "minimum_quality_volume_ids": ["1", "2"],
            "histogram_volume_ids": ["1", "2"],
            "histogram_owned_element_counts": {"hex": 1, "pyramid": 1, "tet": 10},
        },
        "partition_aggregation": {
            "aggregation_policy": "owned_elements_only",
            "reported_global_owned_counts": {"hex": 1, "pyramid": 1, "tet": 10},
            "partitions": [
                {
                    "partition_id": 0,
                    "owned_counts": {"hex": 1, "pyramid": 0, "tet": 4},
                    "ghost_counts": {"hex": 0, "pyramid": 1, "tet": 2},
                },
                {
                    "partition_id": 1,
                    "owned_counts": {"hex": 0, "pyramid": 1, "tet": 6},
                    "ghost_counts": {"hex": 1, "pyramid": 0, "tet": 0},
                },
            ],
        },
        "boundary_sets": [
            {
                "name": "interface",
                "mesh_generation": "mesh-generation-42",
                "mesh_sha256": "d" * 64,
                "entity_ids": [7],
                "connectivity_sha256": "f" * 64,
            }
        ],
        "interface_surfaces": [
            {
                "surface_id": 7,
                "adjacent_volumes": [1, 2],
                "face_ids": [1],
                "face_connectivity": [[2, 1, 3, 4]],
            }
        ],
        "interface_face_ownership": [
            {"face_id": 1, "node_count": 4, "hex_owners": [1], "pyramid_owners": [1]}
        ],
        "matched_pyramid_count": 1,
        "geometry": {
            "cad_total_volume_m3": 2.0e-9,
            "analytic_total_volume_m3": 2.0e-9,
            "element_volume_source": "independent_gmsh_v41_coordinate_reconstruction",
        },
        "gmsh_export": {"bytes": 1233, "sha256": "b" * 64},
        "gmsh_inventory": {
            "status": "ok",
            "mesh_format": "4.1",
            "binary": False,
            "connectivity_mismatches": [],
            "volume_family_counts": {"hex": 1, "pyramid": 1, "tet": 10},
        },
        "gmsh_volume_inventory": {
            "family_counts": {"hex": 1, "pyramid": 1, "tet": 10},
            "family_volumes_m3": {
                "hex": 1.0e-9,
                "pyramid": 2.3570233333333336e-10,
                "tet": 7.642976666666667e-10,
            },
            "total_volume_m3": 2.0e-9,
        },
        "quality_probe": {
            "command_supported": False,
            "diagnostic": "Unknown metric name volume",
            "failure_interpretation": "unsupported_api_not_zero_quality",
            "fallback": "per_element_scaled_jacobian_by_family",
            "families": ["hex", "pyramid", "tet"],
        },
        "process": {
            "exit_code": 3,
            "unexpected_error_lines": [],
            "known_headless_diagnostics_only": True,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
        "export_artifacts": {
            "required": ["mixed.msh", "mixed.vol"],
            "artifacts": [
                {
                    "name": "mixed.msh",
                    "fresh": True,
                    "bytes": 1233,
                    "sha256": "b" * 64,
                },
                {
                    "name": "mixed.vol",
                    "fresh": True,
                    "bytes": 900,
                    "sha256": "c" * 64,
                },
            ],
        },
        "replay_identity": {
            "pinned_journal_sha256": "a" * 64,
            "pinned_source_model_sha256": "c" * 64,
            "replayed_journal_sha256": "a" * 64,
            "replayed_source_model_sha256": "c" * 64,
        },
        "export_manifest": {
            "invocation_id": "batch-invocation-42",
            "model_generation": "model-generation-42",
            "artifacts": [
                {
                    "name": "mixed.msh",
                    "sha256": "b" * 64,
                    "model_generation": "model-generation-42",
                    "invocation_id": "batch-invocation-42",
                },
                {
                    "name": "mixed.vol",
                    "sha256": "c" * 64,
                    "model_generation": "model-generation-42",
                    "invocation_id": "batch-invocation-42",
                },
            ],
        },
        "batch_invocation": {
            "invocation_id": "batch-invocation-42",
            "process_start_utc": "2026-07-16T02:00:00Z",
            "log": {
                "invocation_id": "batch-invocation-42",
                "process_start_utc": "2026-07-16T02:00:00Z",
                "sha256": "1" * 64,
            },
            "exports_invocation_id": "batch-invocation-42",
        },
        "operation_dag_identity": {
            "final_model_generation": "model-generation-42",
            "final_operation_sequence": 6,
            "export_model_generation": "model-generation-42",
            "export_after_operation_sequence": 6,
        },
        "length_scale_identity": {
            "source_geometry_unit": "mm",
            "export_geometry_unit": "m",
            "declared_source_to_export_scale": 0.001,
            "scale_application_stages": ["source-command"],
            "effective_scale": 0.001,
        },
        "signed_jacobian_identity": {
            "mesh_generation": "mesh-generation-42",
            "minimum_signed_jacobian": 0.31,
            "maximum_signed_jacobian": 1.24,
            "interior_sign_change_count": 0,
            "absolute_volume_matches_cad": True,
        },
        "coordinate_scale_identity": {
            "source_geometry_unit": "mm",
            "export_coordinate_unit": "m",
            "coordinate_scale_to_si": 0.001,
            "volume_scale_to_si": 1.0e-9,
            "coordinate_scale_generation": "scale-generation-42",
            "volume_scale_generation": "scale-generation-42",
        },
        "exodus_connectivity_identity": {
            "connectivity_permutation_generation": "exodus-ordering-42",
            "sideset_face_ordinal_generation": "exodus-ordering-42",
            "permuted_connectivity_sha256": "a" * 64,
            "sideset_connectivity_sha256": "a" * 64,
            "target_ordering": "solver-target-ordering-v1",
        },
        "quality_report_generation_identity": {
            "final_mesh_generation": "mesh-generation-42",
            "quality_report_mesh_generation": "mesh-generation-42",
            "final_smoothing_sequence": 8,
            "quality_report_after_operation_sequence": 8,
        },
        "timing_breakdown_s": {
            "source_replay": 0.2,
            "mesh_inventory": 0.3,
            "gmsh_export": 0.2,
            "independent_validation": 0.1,
        },
    }


def test_interface_gate_accepts_minimal_conformal_transition_without_hex_dominance():
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(summary()))
    assert result["status"] == "ok"
    assert result["element_counts"]["tet"] > result["element_counts"]["hex"]
    assert result["checks"]["each_interface_quad_has_one_hex_and_one_pyramid_owner"] is True
    assert result["gmsh_reconstructed_volume_relative_error"] == 0.0


def test_interface_gate_rejects_count_only_false_positive_without_dual_ownership():
    row = copy.deepcopy(summary())
    row["interface_face_ownership"][0]["pyramid_owners"] = []
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["gmsh_volume_families_match_live_inventory"] is True
    assert result["checks"]["each_interface_quad_has_one_hex_and_one_pyramid_owner"] is False


def test_interface_gate_rejects_reconstructed_volume_drift_at_small_si_scale():
    row = copy.deepcopy(summary())
    row["gmsh_volume_inventory"]["total_volume_m3"] = 1.8e-9
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["gmsh_reconstructed_volume_relative_error"] == pytest.approx(0.1)


def test_source_gate_accepts_classified_headless_exit_and_quality_fallback():
    result = json.loads(cubit_mixed_transition_source_gate(summary()))
    assert result["status"] == "ok"
    assert result["warnings"] == []
    assert result["process_exit_code"] == 3
    assert result["checks"]["unsupported_aggregate_quality_probe_is_diagnosed"] is True
    assert result["checks"]["independent_interface_gate_passed"] is True


def test_source_gate_rejects_exit_code_allowlist_without_semantic_evidence():
    row = copy.deepcopy(summary())
    row["process"]["known_headless_diagnostics_only"] = False
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["nonzero_exit_is_semantically_classified"] is False


def test_source_gate_rejects_unsupported_quality_query_misread_as_zero_quality():
    row = copy.deepcopy(summary())
    row["quality_probe"]["failure_interpretation"] = "zero_quality"
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["unsupported_aggregate_quality_probe_is_diagnosed"] is False


def test_server_rejects_missing_independent_gmsh_volume_inventory():
    row = copy.deepcopy(summary())
    del row["gmsh_volume_inventory"]
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "invalid_input"
    assert "gmsh_volume_inventory" in result["error"]


def test_legacy_source_identity_is_accepted_with_warnings():
    row = summary()
    row.pop("export_artifacts")
    row.pop("replay_identity")
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "ok"
    assert set(result["warnings"]) == {
        "per_artifact_export_freshness_not_recorded",
        "journal_model_replay_identity_not_recorded",
    }


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_hex_orientation_negative_jacobian",
        "v7_public_pyramid_transition_nonmanifold",
    ],
)
def test_generalization_v7_public(case_id: str):
    row = summary()
    if case_id == "v7_public_hex_orientation_negative_jacobian":
        row["quality"]["hex"]["scaled_jacobian"]["min"] = -0.01
        expected = "all_volume_families_above_quality_threshold"
    else:
        row["interface_face_ownership"][0]["tet_owners"] = [7]
        row["interface_surfaces"][0]["face_incidence_count"] = 3
        expected = "interface_quads_are_two_sided_manifold"
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_source_batch_partial_export_success",
        "v7_source_journal_model_digest_mismatch",
    ],
)
def test_generalization_v7_source(case_id: str):
    row = summary()
    if case_id == "v7_source_batch_partial_export_success":
        row["process"]["exit_code"] = 0
        row["export_artifacts"]["artifacts"][1]["fresh"] = False
        expected = "all_required_export_artifacts_are_fresh"
    else:
        row["replay_identity"]["replayed_source_model_sha256"] = "d" * 64
        expected = "journal_and_source_model_identity_match_replay"
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize("duplicate_location", ["required", "artifacts"])
def test_source_gate_rejects_duplicate_export_artifact_names(duplicate_location: str):
    row = summary()
    if duplicate_location == "required":
        row["export_artifacts"]["required"][1] = "mixed.msh"
    else:
        row["export_artifacts"]["artifacts"][1]["name"] = "mixed.msh"
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_required_export_artifacts_are_fresh"] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v8_public_quality_report_older_than_mesh",
        "v8_public_sideset_generation_mixed_after_remesh",
    ],
)
def test_generalization_v8_public(case_id: str):
    row = summary()
    if case_id == "v8_public_quality_report_older_than_mesh":
        row["quality_report_identity"]["mesh_generation"] = "mesh-generation-41"
        row["quality_report_identity"]["mesh_sha256"] = "2" * 64
        expected = "quality_report_matches_current_mesh_generation"
    else:
        row["boundary_sets"][0]["mesh_generation"] = "mesh-generation-41"
        row["boundary_sets"][0]["mesh_sha256"] = "2" * 64
        expected = "boundary_sets_match_current_mesh_generation"
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v8_source_export_manifest_mixed_generations",
        "v8_source_batch_log_prior_invocation_identity",
    ],
)
def test_generalization_v8_source(case_id: str):
    row = summary()
    if case_id == "v8_source_export_manifest_mixed_generations":
        row["export_manifest"]["artifacts"][1][
            "model_generation"
        ] = "model-generation-41"
        expected = "export_manifest_uses_one_model_and_invocation_generation"
    else:
        row["batch_invocation"]["log"].update(
            {
                "invocation_id": "batch-invocation-41",
                "process_start_utc": "2026-07-16T01:00:00Z",
            }
        )
        expected = "batch_log_and_exports_share_invocation_identity"
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v9_public_quality_histogram_entity_scope_mismatch",
        "v9_public_partition_ghost_elements_double_counted",
    ],
)
def test_generalization_v9_public(case_id: str):
    row = summary()
    if case_id == "v9_public_quality_histogram_entity_scope_mismatch":
        row["quality_scope_identity"]["histogram_volume_ids"] = ["1"]
        expected = "quality_histogram_covers_the_complete_mesh_scope"
    else:
        row["partition_aggregation"].update(
            {
                "aggregation_policy": "owned_plus_ghost_elements",
                "reported_global_owned_counts": {
                    "hex": 2,
                    "pyramid": 2,
                    "tet": 12,
                },
            }
        )
        expected = "partition_aggregation_excludes_ghost_elements"
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v9_source_export_precedes_final_geometry_operation",
        "v9_source_export_length_scale_applied_twice",
    ],
)
def test_generalization_v9_source(case_id: str):
    row = summary()
    if case_id == "v9_source_export_precedes_final_geometry_operation":
        row["operation_dag_identity"].update(
            {
                "export_model_generation": "model-generation-41",
                "export_after_operation_sequence": 5,
            }
        )
        expected = "exports_follow_the_final_geometry_operation"
    else:
        row["length_scale_identity"].update(
            {
                "scale_application_stages": ["source-command", "export-manifest"],
                "effective_scale": 1.0e-6,
            }
        )
        expected = "length_scale_is_applied_exactly_once"
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v10_public_signed_jacobian_folded_hex_abs_volume_passes",
        "v10_public_export_coordinate_scale_mesh_cad_mismatch",
    ],
)
def test_generalization_v10_public(case_id: str):
    row = summary()
    if case_id == "v10_public_signed_jacobian_folded_hex_abs_volume_passes":
        row["signed_jacobian_identity"].update(
            {
                "minimum_signed_jacobian": -0.08,
                "interior_sign_change_count": 3,
            }
        )
        expected = "signed_jacobians_remain_positive_inside_high_order_hexes"
    else:
        row["coordinate_scale_identity"].update(
            {
                "coordinate_scale_to_si": 1.0,
                "coordinate_scale_generation": "scale-generation-41",
            }
        )
        expected = "mesh_coordinates_and_volume_use_one_length_scale"
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v10_source_exodus_connectivity_permuted_sideset_stale",
        "v10_source_quality_report_before_final_smoothing",
    ],
)
def test_generalization_v10_source(case_id: str):
    row = summary()
    if case_id == "v10_source_exodus_connectivity_permuted_sideset_stale":
        row["exodus_connectivity_identity"].update(
            {
                "sideset_face_ordinal_generation": "exodus-ordering-41",
                "sideset_connectivity_sha256": "b" * 64,
            }
        )
        expected = "exodus_sidesets_follow_connectivity_permutation"
    else:
        row["quality_report_generation_identity"].update(
            {
                "quality_report_mesh_generation": "mesh-generation-41",
                "quality_report_after_operation_sequence": 7,
            }
        )
        expected = "quality_report_follows_final_smoothing_generation"
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


def _with_v11_ownership(row: dict) -> dict:
    row["high_order_shared_face_orientation_identity"] = {
        "mesh_generation": "mesh-generation-43",
        "left_face_node_ids": [1, 2, 3, 4, 5, 6, 7, 8],
        "right_face_node_ids": [1, 4, 3, 2, 8, 7, 6, 5],
        "left_element_generation": "mesh-generation-43",
        "right_element_generation": "mesh-generation-43",
    }
    row["live_cad_mesh_identity"] = {
        "live_cad_sha256": "1" * 64,
        "mesh_source_cad_sha256": "1" * 64,
        "live_cad_generation": "cad-generation-43",
        "mesh_source_cad_generation": "cad-generation-43",
        "live_cad_volume": 12.5,
        "mesh_reference_cad_volume": 12.5,
    }
    row["block_material_map_identity"] = {
        "final_mesh_generation": "mesh-generation-43",
        "material_map_mesh_generation": "mesh-generation-43",
        "block_table_sha256": "2" * 64,
        "material_map_block_table_sha256": "2" * 64,
        "unmapped_block_ids": [],
    }
    row["parallel_sculpt_completion_identity"] = {
        "expected_rank_count": 4,
        "finalized_rank_ids": [0, 1, 2, 3],
        "rank_artifact_sha256": ["3" * 64, "4" * 64, "5" * 64, "6" * 64],
        "rank_manifest_generation": "sculpt-generation-43",
        "global_aggregation_generation": "sculpt-generation-43",
    }
    return row


@pytest.mark.parametrize(
    "case_id",
    [
        "v11_public_high_order_hex_shared_face_orientation_mismatch",
        "v11_public_mesh_volume_live_cad_digest_mismatch",
    ],
)
def test_generalization_v11_public(case_id: str):
    row = _with_v11_ownership(summary())
    if case_id == "v11_public_high_order_hex_shared_face_orientation_mismatch":
        row["high_order_shared_face_orientation_identity"][
            "right_face_node_ids"
        ] = [1, 2, 3, 4, 5, 6, 7, 8]
        expected = "high_order_hex_shared_faces_have_reciprocal_orientation"
    else:
        row["live_cad_mesh_identity"].update(
            {
                "mesh_source_cad_sha256": "7" * 64,
                "mesh_source_cad_generation": "cad-generation-42",
            }
        )
        expected = "mesh_manifest_matches_live_cad_identity"
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v11_source_block_material_map_previous_export_generation",
        "v11_source_parallel_sculpt_partial_rank_artifact_aggregated",
    ],
)
def test_generalization_v11_source(case_id: str):
    row = _with_v11_ownership(summary())
    if case_id == "v11_source_block_material_map_previous_export_generation":
        row["block_material_map_identity"].update(
            {
                "material_map_mesh_generation": "mesh-generation-42",
                "material_map_block_table_sha256": "8" * 64,
            }
        )
        expected = "block_material_map_matches_final_mesh_generation"
    else:
        row["parallel_sculpt_completion_identity"].update(
            {
                "finalized_rank_ids": [0, 1, 3],
                "rank_artifact_sha256": ["3" * 64, "4" * 64, "6" * 64],
            }
        )
        expected = "parallel_sculpt_waits_for_every_rank_artifact"
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False
