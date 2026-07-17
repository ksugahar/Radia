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


def _with_v12_topology_generation_identity(row: dict) -> dict:
    row = _with_v11_ownership(row)
    row["webcut_sideset_topology_identity"] = {
        "final_geometry_operation_sequence": 9,
        "sideset_capture_after_operation_sequence": 9,
        "final_geometry_generation": "geometry-generation-44",
        "sideset_geometry_generation": "geometry-generation-44",
        "final_topology_generation": "topology-generation-44",
        "sideset_topology_generation": "topology-generation-44",
        "sideset_surface_ids": [21, 22],
        "resolved_surface_ids": [21, 22],
        "sideset_connectivity_sha256": "7" * 64,
        "resolved_connectivity_sha256": "7" * 64,
    }
    row["high_order_curved_node_identity"] = {
        "export_order": 2,
        "final_mesh_generation": "mesh-generation-44",
        "edge_node_mesh_generation": "mesh-generation-44",
        "face_node_mesh_generation": "mesh-generation-44",
        "curving_generation": "curving-generation-44",
        "export_curving_generation": "curving-generation-44",
        "high_order_edge_node_count": 48,
        "high_order_face_node_count": 36,
    }
    row["headless_sideset_manifest_identity"] = {
        "batch_invocation_id": "batch-invocation-44",
        "manifest_invocation_id": "batch-invocation-44",
        "final_webcut_operation_sequence": 9,
        "manifest_capture_after_operation_sequence": 9,
        "final_topology_generation": "topology-generation-44",
        "manifest_topology_generation": "topology-generation-44",
        "manifest_surface_ids": [21, 22],
        "live_surface_ids": [21, 22],
        "manifest_connectivity_sha256": "8" * 64,
        "live_connectivity_sha256": "8" * 64,
    }
    row["netgen_high_order_export_identity"] = {
        "export_order": 2,
        "final_mesh_generation": "mesh-generation-44",
        "higher_order_node_mesh_generation": "mesh-generation-44",
        "export_model_generation": "model-generation-44",
        "active_model_generation": "model-generation-44",
        "higher_order_node_count": 84,
        "netgen_export_sha256": "9" * 64,
    }
    return row


def _with_v13_block_orientation_identity(row: dict) -> dict:
    row = _with_v12_topology_generation_identity(row)
    row["hex_block_material_topology_identity"] = {
        "final_imprint_generation": "imprint-generation-45", "final_topology_generation": "topology-generation-45",
        "block_topology_generation": "topology-generation-45", "material_assignment_topology_generation": "topology-generation-45",
        "block_volume_ids": [1, 2, 3], "material_assignment_volume_ids": [1, 2, 3],
        "block_material_map_sha256": "c" * 64, "resolved_material_map_sha256": "c" * 64,
    }
    row["pyramid_transition_face_orientation_identity"] = {
        "transition_generation": "transition-generation-45", "pyramid_face_generation": "transition-generation-45", "hex_face_generation": "transition-generation-45",
        "shared_face_node_ids": [101, 102, 103, 104], "pyramid_face_node_ids": [101, 102, 103, 104], "hex_face_node_ids": [104, 103, 102, 101], "opposed_outward_normal_dot": -1.0,
    }
    row["headless_block_material_manifest_identity"] = {
        "batch_invocation_id": "batch-invocation-45", "manifest_invocation_id": "batch-invocation-45", "final_imprint_generation": "imprint-generation-45",
        "active_topology_generation": "topology-generation-45", "manifest_topology_generation": "topology-generation-45",
        "manifest_volume_ids": [1, 2, 3], "live_volume_ids": [1, 2, 3], "manifest_material_map_sha256": "d" * 64, "live_material_map_sha256": "d" * 64,
    }
    row["mesh_export_transition_orientation_identity"] = {
        "export_generation": "mesh-export-45", "pyramid_face_export_generation": "mesh-export-45", "hex_face_export_generation": "mesh-export-45",
        "shared_face_node_ids": [101, 102, 103, 104], "pyramid_face_node_ids": [101, 102, 103, 104], "hex_face_node_ids": [104, 103, 102, 101], "opposed_outward_normal_dot": -1.0,
    }
    return row


def _with_v14_high_order_id_width_identity(row: dict) -> dict:
    row = _with_v13_block_orientation_identity(row)
    row["high_order_hex_curved_node_ordering_identity"] = {
        "high_order_mesh_generation": "high-order-mesh-46",
        "curved_geometry_generation": "curved-geometry-46",
        "element_geometry_generation": "curved-geometry-46",
        "element_type": "hex20_serendipity",
        "export_element_type": "hex20_serendipity",
        "node_ordering_convention": "cubit_hex20",
        "export_node_ordering_convention": "cubit_hex20",
        "canonical_node_ids": list(range(101, 121)),
        "export_node_ids": list(range(101, 121)),
        "canonical_node_order_sha256": "a" * 64,
        "export_node_order_sha256": "a" * 64,
    }
    row["sideset_outward_normal_merge_identity"] = {
        "final_merge_generation": "merge-generation-46",
        "sideset_topology_generation": "merge-generation-46",
        "normal_owner_topology_generation": "merge-generation-46",
        "sideset_face_ids": [201, 202, 203],
        "normal_owner_face_ids": [201, 202, 203],
        "owner_volume_ids": [1, 1, 2],
        "resolved_owner_volume_ids": [1, 1, 2],
        "outward_normal_signs": [1, 1, 1],
    }
    row["journal_entity_id_map_reset_identity"] = {
        "reset_generation": "reset-generation-46",
        "journal_replay_reset_generation": "reset-generation-46",
        "entity_id_map_reset_generation": "reset-generation-46",
        "entity_kinds": ["volume", "surface", "curve"],
        "requested_entity_ids": [1, 7, 19],
        "resolved_entity_ids": [1, 7, 19],
        "entity_id_map_sha256": "b" * 64,
        "resolved_entity_id_map_sha256": "b" * 64,
    }
    row["exodus_entity_id_width_identity"] = {
        "export_generation": "exodus-export-46",
        "decoder_export_generation": "exodus-export-46",
        "declared_entity_id_width_bits": 64,
        "decoder_entity_id_width_bits": 64,
        "integer_storage_type": "int64",
        "maximum_entity_id": 4294967311,
        "decoded_maximum_entity_id": 4294967311,
        "entity_id_stream_sha256": "c" * 64,
        "decoded_entity_id_stream_sha256": "c" * 64,
    }
    return row


def _with_v15_interface_quadrature_identity(row: dict) -> dict:
    row = _with_v14_high_order_id_width_identity(row)
    row["mixed_interface_smoothing_orientation_identity"] = {
        "smoothing_generation": "smoothing-47",
        "interface_face_orientation_generation": "smoothing-47",
        "interface_topology_generation": "interface-topology-47",
        "orientation_topology_generation": "interface-topology-47",
        "hex_interface_face_ids": [301, 302],
        "pyramid_interface_face_ids": [401, 402],
        "paired_orientation_products": [-1, -1],
        "interface_pair_sha256": "1" * 64,
        "oriented_interface_pair_sha256": "1" * 64,
    }
    row["high_order_jacobian_quadrature_identity"] = {
        "high_order_mesh_generation": "high-order-mesh-47",
        "quality_evaluation_mesh_generation": "high-order-mesh-47",
        "element_order": 2,
        "required_jacobian_exactness_degree": 4,
        "jacobian_quadrature_exactness_degree": 4,
        "quadrature_rule_generation": "jacobian-quadrature-47",
        "quality_evaluation_quadrature_generation": "jacobian-quadrature-47",
        "element_geometry_sha256": "2" * 64,
        "quality_evaluation_geometry_sha256": "2" * 64,
    }
    row["imprint_merge_tolerance_unit_identity"] = {
        "geometry_generation": "geometry-47",
        "imprint_geometry_generation": "geometry-47",
        "merge_geometry_generation": "geometry-47",
        "tolerance_generation": "tolerance-47",
        "imprint_tolerance_generation": "tolerance-47",
        "merge_tolerance_generation": "tolerance-47",
        "model_length_unit": "mm",
        "imprint_tolerance_unit": "mm",
        "merge_tolerance_unit": "mm",
        "imprint_tolerance_value": 1.0e-6,
        "merge_tolerance_value": 1.0e-6,
        "tolerance_si_m": 1.0e-9,
    }
    row["exodus_block_sideset_renumber_identity"] = {
        "renumber_generation": "renumber-47",
        "block_map_generation": "renumber-47",
        "sideset_map_generation": "renumber-47",
        "block_ids": [11, 12],
        "exported_block_ids": [11, 12],
        "sideset_ids": [21, 22, 23],
        "exported_sideset_ids": [21, 22, 23],
        "entity_map_sha256": "3" * 64,
        "exported_entity_map_sha256": "3" * 64,
    }
    return row


def _with_v16_sweep_jacobian_and_source_identity(row: dict) -> dict:
    row = _with_v15_interface_quadrature_identity(row)
    row["hex_sweep_vertex_correspondence_heal_identity"] = {
        "geometry_heal_generation": "heal-50",
        "sweep_geometry_heal_generation": "heal-50",
        "source_vertex_map_heal_generation": "heal-50",
        "target_vertex_map_heal_generation": "heal-50",
        "source_vertex_ids": [101, 102, 103, 104],
        "target_vertex_ids": [201, 202, 203, 204],
        "sweep_source_vertex_ids": [101, 102, 103, 104],
        "sweep_target_vertex_ids": [201, 202, 203, 204],
        "vertex_correspondence_sha256": "5" * 64,
        "sweep_vertex_correspondence_sha256": "5" * 64,
    }
    row["transition_jacobian_parent_orientation_identity"] = {
        "transition_mesh_generation": "transition-mesh-50",
        "jacobian_mesh_generation": "transition-mesh-50",
        "parent_orientation_generation": "transition-mesh-50",
        "parent_orientation_convention": "right_handed_positive",
        "jacobian_orientation_convention": "right_handed_positive",
        "minimum_signed_jacobian": 0.125,
        "minimum_absolute_jacobian": 0.125,
        "parent_orientation_sha256": "6" * 64,
        "jacobian_parent_orientation_sha256": "6" * 64,
    }
    row["journal_entity_id_imprint_identity"] = {
        "imprint_generation": "imprint-50",
        "journal_entity_generation": "imprint-50",
        "resolved_entity_generation": "imprint-50",
        "journal_volume_ids": [11, 12],
        "resolved_volume_ids": [11, 12],
        "journal_surface_ids": [21, 22, 23],
        "resolved_surface_ids": [21, 22, 23],
        "entity_table_sha256": "7" * 64,
        "resolved_entity_table_sha256": "7" * 64,
    }
    row["exodus_sideset_outward_normal_topology_identity"] = {
        "topology_generation": "topology-50",
        "sideset_map_topology_generation": "topology-50",
        "normal_ownership_topology_generation": "topology-50",
        "sideset_ids": [31, 32],
        "normal_ownership_sideset_ids": [31, 32],
        "normal_orientation": "outward",
        "exported_normal_orientation": "outward",
        "normal_ownership_sha256": "8" * 64,
        "exported_normal_ownership_sha256": "8" * 64,
    }
    return row


def _with_v17_periodic_transition_and_export_identity(row: dict) -> dict:
    row = _with_v16_sweep_jacobian_and_source_identity(row)
    row["periodic_hex_node_pair_transform_frame_identity"] = {
        "mesh_generation": "periodic-hex-mesh-51",
        "node_pair_mesh_generation": "periodic-hex-mesh-51",
        "periodic_transform_frame_mesh_generation": "periodic-hex-mesh-51",
        "periodic_transform_generation": "periodic-transform-51",
        "node_pair_periodic_transform_generation": "periodic-transform-51",
        "source_node_ids": [101, 102, 103, 104],
        "target_node_ids": [201, 202, 203, 204],
        "paired_source_node_ids": [101, 102, 103, 104],
        "paired_target_node_ids": [201, 202, 203, 204],
        "coordinate_frame": "global_cartesian",
        "node_pair_coordinate_frame": "global_cartesian",
        "transform_matrix_sha256": "1" * 64,
        "applied_transform_matrix_sha256": "1" * 64,
    }
    row["pyramid_transition_face_diagonal_convention_identity"] = {
        "transition_mesh_generation": "transition-mesh-51",
        "tet_neighbor_mesh_generation": "transition-mesh-51",
        "hex_neighbor_mesh_generation": "transition-mesh-51",
        "pyramid_face_ids": [301, 302],
        "tet_neighbor_face_ids": [301, 302],
        "hex_neighbor_face_ids": [301, 302],
        "diagonal_convention": "canonical_node_0_to_2",
        "tet_neighbor_diagonal_convention": "canonical_node_0_to_2",
        "hex_neighbor_diagonal_convention": "canonical_node_0_to_2",
        "transition_face_connectivity_sha256": "2" * 64,
        "neighbor_face_connectivity_sha256": "2" * 64,
    }
    row["block_attribute_material_id_merge_identity"] = {
        "final_merge_generation": "merge-51",
        "block_attribute_merge_generation": "merge-51",
        "material_id_map_merge_generation": "merge-51",
        "block_ids": [11, 12],
        "block_attribute_block_ids": [11, 12],
        "material_ids": [101, 102],
        "exported_material_ids": [101, 102],
        "block_material_map_sha256": "3" * 64,
        "exported_block_material_map_sha256": "3" * 64,
    }
    row["high_order_exodus_node_permutation_export_order_identity"] = {
        "export_order_generation": "exodus-order-51",
        "permutation_table_export_order_generation": "exodus-order-51",
        "writer_export_order_generation": "exodus-order-51",
        "element_order": 2,
        "source_node_order": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "permutation_table": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "written_node_order": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "node_permutation_sha256": "4" * 64,
        "written_node_permutation_sha256": "4" * 64,
    }
    return row


def _with_v18_ordinal_bias_and_quality_identity(row: dict) -> dict:
    row = _with_v17_periodic_transition_and_export_identity(row)
    row["hex_sideset_outward_face_ordinal_volume_reorder_identity"] = {
        "mesh_generation": "hex-mesh-52",
        "volume_connectivity_reorder_generation": "hex-reorder-52",
        "face_ordinal_mesh_generation": "hex-mesh-52",
        "normal_ownership_mesh_generation": "hex-mesh-52",
        "face_ordinal_connectivity_reorder_generation": "hex-reorder-52",
        "normal_connectivity_reorder_generation": "hex-reorder-52",
        "element_ids": [101, 102],
        "face_ordinals": [5, 6],
        "exported_element_ids": [101, 102],
        "exported_face_ordinals": [5, 6],
        "outward_normal_signs": [1, 1],
        "exported_outward_normal_signs": [1, 1],
        "element_face_map_sha256": "1" * 64,
        "exported_element_face_map_sha256": "1" * 64,
    }
    row["sweep_layer_bias_source_curve_orientation_generation_identity"] = {
        "sweep_generation": "sweep-52",
        "source_curve_orientation_generation": "curve-orientation-52",
        "layer_bias_sweep_generation": "sweep-52",
        "layer_bias_curve_orientation_generation": "curve-orientation-52",
        "source_curve_ids": [31, 32],
        "source_curve_orientations": [1, -1],
        "biased_curve_ids": [31, 32],
        "biased_curve_orientations": [1, -1],
        "interval_counts": [4, 8],
        "biased_interval_counts": [4, 8],
        "bias_factors": [1.2, 1.5],
        "applied_bias_factors": [1.2, 1.5],
        "curve_bias_map_sha256": "2" * 64,
        "applied_curve_bias_map_sha256": "2" * 64,
    }
    row["exodus_sideset_element_face_topology_generation_identity"] = {
        "mesh_generation": "mesh-52",
        "exodus_export_generation": "exodus-52",
        "topology_ordinal_mesh_generation": "mesh-52",
        "writer_mesh_generation": "mesh-52",
        "topology_ordinal_export_generation": "exodus-52",
        "writer_export_generation": "exodus-52",
        "element_ids": [101, 102],
        "element_face_topology_ordinals": [5, 6],
        "written_element_ids": [101, 102],
        "written_element_face_topology_ordinals": [5, 6],
        "element_face_topology_sha256": "3" * 64,
        "written_element_face_topology_sha256": "3" * 64,
    }
    row["high_order_quality_reference_coordinate_generation_identity"] = {
        "mesh_generation": "high-order-mesh-52",
        "element_order_generation": "element-order-52",
        "reference_node_mesh_generation": "high-order-mesh-52",
        "quality_mesh_generation": "high-order-mesh-52",
        "reference_node_element_order_generation": "element-order-52",
        "quality_element_order_generation": "element-order-52",
        "element_order": 2,
        "reference_node_count": 10,
        "quality_reference_node_count": 10,
        "jacobian_sampling_rule": "tet10_reference_nodes_and_interior",
        "quality_jacobian_sampling_rule": "tet10_reference_nodes_and_interior",
        "reference_coordinates_sha256": "4" * 64,
        "quality_reference_coordinates_sha256": "4" * 64,
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


def test_v12_public_hex_sideset_after_webcut_topology_generation_mismatch():
    row = _with_v12_topology_generation_identity(summary())
    row["webcut_sideset_topology_identity"].update(
        {
            "sideset_topology_generation": "topology-generation-43",
            "sideset_connectivity_sha256": "a" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["sidesets_follow_final_webcut_topology"] is False


def test_v12_public_high_order_export_curved_node_generation_mismatch():
    row = _with_v12_topology_generation_identity(summary())
    row["high_order_curved_node_identity"].update(
        {
            "edge_node_mesh_generation": "mesh-generation-43",
            "face_node_mesh_generation": "mesh-generation-43",
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "high_order_export_nodes_match_current_mesh_and_curving_generation"
        ]
        is False
    )


def test_v12_source_hex_sideset_after_webcut_topology_generation_mismatch():
    row = _with_v12_topology_generation_identity(summary())
    row["headless_sideset_manifest_identity"].update(
        {
            "manifest_capture_after_operation_sequence": 8,
            "manifest_topology_generation": "topology-generation-43",
            "manifest_connectivity_sha256": "b" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "headless_sideset_manifest_follows_final_webcut_topology"
        ]
        is False
    )


def test_v12_source_high_order_export_curved_node_generation_mismatch():
    row = _with_v12_topology_generation_identity(summary())
    row["netgen_high_order_export_identity"].update(
        {
            "higher_order_node_mesh_generation": "mesh-generation-43",
            "export_model_generation": "model-generation-43",
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "netgen_high_order_nodes_match_current_mesh_and_model_generation"
        ]
        is False
    )


def test_v13_public_hex_block_material_assignment_after_imprint_generation_mismatch():
    row = _with_v13_block_orientation_identity(summary())
    row["hex_block_material_topology_identity"].update({"material_assignment_topology_generation": "topology-generation-44", "resolved_material_map_sha256": "e" * 64})
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["hex_block_materials_follow_final_imprint_topology"] is False


def test_v13_public_pyramid_transition_face_orientation_mismatch():
    row = _with_v13_block_orientation_identity(summary())
    row["pyramid_transition_face_orientation_identity"].update({"hex_face_node_ids": [101, 102, 103, 104], "opposed_outward_normal_dot": 1.0})
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["pyramid_hex_transition_faces_have_opposed_orientation"] is False


def test_v13_source_hex_block_material_assignment_after_imprint_generation_mismatch():
    row = _with_v13_block_orientation_identity(summary())
    row["headless_block_material_manifest_identity"].update({"manifest_topology_generation": "topology-generation-44", "manifest_material_map_sha256": "f" * 64})
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["headless_block_material_manifest_follows_final_imprint"] is False


def test_v13_source_pyramid_transition_face_orientation_mismatch():
    row = _with_v13_block_orientation_identity(summary())
    row["mesh_export_transition_orientation_identity"].update({"hex_face_node_ids": [101, 102, 103, 104], "opposed_outward_normal_dot": 1.0})
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["mesh_export_transition_faces_have_opposed_orientation"] is False


def test_v14_public_high_order_hex_curved_node_ordering_mismatch():
    row = _with_v14_high_order_id_width_identity(summary())
    export_ids = list(range(101, 121))
    export_ids[8], export_ids[9] = export_ids[9], export_ids[8]
    row["high_order_hex_curved_node_ordering_identity"].update(
        {
            "export_node_ordering_convention": "vtk_quadratic_hex",
            "export_node_ids": export_ids,
            "export_node_order_sha256": "d" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["curved_high_order_hex_uses_canonical_node_ordering"] is False


def test_v14_public_sideset_outward_normal_after_merge_generation_mismatch():
    row = _with_v14_high_order_id_width_identity(summary())
    row["sideset_outward_normal_merge_identity"].update(
        {
            "normal_owner_topology_generation": "merge-generation-45",
            "resolved_owner_volume_ids": [2, 1, 2],
            "outward_normal_signs": [-1, 1, 1],
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["merged_sideset_normals_follow_final_topology_owners"] is False


def test_v14_source_journal_entity_id_map_previous_reset_generation():
    row = _with_v14_high_order_id_width_identity(summary())
    row["journal_entity_id_map_reset_identity"].update(
        {
            "entity_id_map_reset_generation": "reset-generation-45",
            "resolved_entity_ids": [2, 8, 20],
            "resolved_entity_id_map_sha256": "d" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["journal_entity_ids_follow_current_reset_generation"] is False


def test_v14_source_exodus_64bit_entity_id_width_truncation():
    row = _with_v14_high_order_id_width_identity(summary())
    row["exodus_entity_id_width_identity"].update(
        {
            "decoder_entity_id_width_bits": 32,
            "integer_storage_type": "int32",
            "decoded_maximum_entity_id": 15,
            "decoded_entity_id_stream_sha256": "d" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["exodus_decoder_preserves_declared_64bit_entity_ids"] is False


def test_v15_positive_interface_quadrature_and_source_lineage():
    row = _with_v15_interface_quadrature_identity(summary())
    public = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    source = json.loads(cubit_mixed_transition_source_gate(row))
    assert public["status"] == "ok"
    assert source["status"] == "ok"


def test_v15_public_hex_pyramid_interface_face_orientation_after_smoothing_mismatch():
    row = _with_v15_interface_quadrature_identity(summary())
    row["mixed_interface_smoothing_orientation_identity"].update(
        {
            "interface_face_orientation_generation": "smoothing-46",
            "paired_orientation_products": [1, -1],
            "oriented_interface_pair_sha256": "4" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "smoothed_hex_pyramid_interface_faces_keep_opposed_orientation"
        ]
        is False
    )


def test_v15_public_high_order_jacobian_quadrature_rule_generation_mismatch():
    row = _with_v15_interface_quadrature_identity(summary())
    row["high_order_jacobian_quadrature_identity"].update(
        {
            "jacobian_quadrature_exactness_degree": 2,
            "quality_evaluation_quadrature_generation": "jacobian-quadrature-46",
            "quality_evaluation_geometry_sha256": "4" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["high_order_jacobian_uses_current_sufficient_quadrature"]
        is False
    )


def test_v15_source_imprint_merge_tolerance_length_unit_basis_mismatch():
    row = _with_v15_interface_quadrature_identity(summary())
    row["imprint_merge_tolerance_unit_identity"].update(
        {
            "merge_tolerance_generation": "tolerance-46",
            "merge_tolerance_unit": "m",
            "merge_tolerance_value": 1.0e-6,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["imprint_merge_tolerances_share_one_physical_length_basis"]
        is False
    )


def test_v15_source_exodus_block_sideset_map_previous_renumber_generation():
    row = _with_v15_interface_quadrature_identity(summary())
    row["exodus_block_sideset_renumber_identity"].update(
        {
            "block_map_generation": "renumber-46",
            "exported_block_ids": [1, 2],
            "exported_entity_map_sha256": "4" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "exodus_block_sideset_maps_follow_current_renumber_generation"
        ]
        is False
    )


def test_v16_positive_sweep_jacobian_and_source_lineage():
    row = _with_v16_sweep_jacobian_and_source_identity(summary())
    public = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    source = json.loads(cubit_mixed_transition_source_gate(row))
    assert public["status"] == "ok"
    assert source["status"] == "ok"


def test_v16_public_hex_sweep_source_target_vertex_map_after_heal_mismatch():
    row = _with_v16_sweep_jacobian_and_source_identity(summary())
    row["hex_sweep_vertex_correspondence_heal_identity"].update(
        {
            "source_vertex_map_heal_generation": "heal-49",
            "sweep_source_vertex_ids": [101, 103, 102, 104],
            "sweep_vertex_correspondence_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["hex_sweep_uses_post_heal_vertex_correspondence"] is False


def test_v16_public_transition_element_jacobian_parent_orientation_convention_mismatch():
    row = _with_v16_sweep_jacobian_and_source_identity(summary())
    row["transition_jacobian_parent_orientation_identity"].update(
        {
            "parent_orientation_generation": "transition-mesh-49",
            "jacobian_orientation_convention": "left_handed_negative",
            "minimum_signed_jacobian": -0.125,
            "jacobian_parent_orientation_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["transition_jacobian_uses_current_parent_orientation"]
        is False
    )


def test_v16_source_journal_entity_ids_previous_imprint_generation():
    row = _with_v16_sweep_jacobian_and_source_identity(summary())
    row["journal_entity_id_imprint_identity"].update(
        {
            "journal_entity_generation": "imprint-49",
            "resolved_surface_ids": [22, 23, 24],
            "resolved_entity_table_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["journal_entity_ids_follow_final_imprint_generation"] is False


def test_v16_source_exodus_sideset_outward_normal_topology_generation_mismatch():
    row = _with_v16_sweep_jacobian_and_source_identity(summary())
    row["exodus_sideset_outward_normal_topology_identity"].update(
        {
            "normal_ownership_topology_generation": "topology-49",
            "exported_normal_orientation": "inward",
            "exported_normal_ownership_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["exodus_sideset_normals_follow_current_topology"] is False


def test_v17_positive_periodic_transition_and_export_identity():
    row = _with_v17_periodic_transition_and_export_identity(summary())
    public = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    source = json.loads(cubit_mixed_transition_source_gate(row))
    assert public["status"] == "ok"
    assert source["status"] == "ok"


def test_v17_public_periodic_hex_node_pair_transform_frame_generation_mismatch():
    row = _with_v17_periodic_transition_and_export_identity(summary())
    row["periodic_hex_node_pair_transform_frame_identity"].update(
        {
            "periodic_transform_frame_mesh_generation": "periodic-hex-mesh-50",
            "node_pair_periodic_transform_generation": "periodic-transform-50",
            "node_pair_coordinate_frame": "periodic_local_previous",
            "applied_transform_matrix_sha256": "5" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["periodic_hex_node_pairs_use_current_transform_frame"]
        is False
    )


def test_v17_public_pyramid_transition_face_diagonal_convention_neighbor_mismatch():
    row = _with_v17_periodic_transition_and_export_identity(summary())
    row["pyramid_transition_face_diagonal_convention_identity"].update(
        {
            "hex_neighbor_diagonal_convention": "canonical_node_1_to_3",
            "neighbor_face_connectivity_sha256": "5" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "pyramid_transition_neighbors_share_one_face_diagonal_convention"
        ]
        is False
    )


def test_v17_source_block_attribute_material_id_previous_merge_generation():
    row = _with_v17_periodic_transition_and_export_identity(summary())
    row["block_attribute_material_id_merge_identity"].update(
        {
            "block_attribute_merge_generation": "merge-50",
            "exported_material_ids": [102, 101],
            "exported_block_material_map_sha256": "5" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["block_material_attributes_follow_final_merge_generation"]
        is False
    )


def test_v17_source_high_order_exodus_node_permutation_previous_export_order_generation():
    row = _with_v17_periodic_transition_and_export_identity(summary())
    row["high_order_exodus_node_permutation_export_order_identity"].update(
        {
            "permutation_table_export_order_generation": "exodus-order-50",
            "permutation_table": [1, 3, 2, 4, 5, 6, 7, 8, 9, 10],
            "written_node_permutation_sha256": "5" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "high_order_exodus_nodes_use_current_export_order_permutation"
        ]
        is False
    )


def test_v18_positive_ordinal_bias_and_quality_identity():
    row = _with_v18_ordinal_bias_and_quality_identity(summary())
    public = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    source = json.loads(cubit_mixed_transition_source_gate(row))
    assert public["status"] == "ok"
    assert source["status"] == "ok"
    assert public["checks"][
        "hex_sideset_face_ordinals_and_normals_follow_connectivity_reorder"
    ]
    assert public["checks"][
        "biased_sweep_layers_follow_current_source_curve_orientation"
    ]
    assert source["checks"][
        "exodus_sideset_ordinals_follow_current_mesh_and_export_topology"
    ]
    assert source["checks"][
        "high_order_quality_uses_current_reference_coordinates_and_order"
    ]


def test_v18_public_hex_sideset_outward_face_ordinal_volume_reorder_mismatch():
    row = _with_v18_ordinal_bias_and_quality_identity(summary())
    row["hex_sideset_outward_face_ordinal_volume_reorder_identity"].update(
        {
            "face_ordinal_connectivity_reorder_generation": "hex-reorder-51",
            "normal_connectivity_reorder_generation": "hex-reorder-51",
            "exported_face_ordinals": [6, 5],
            "exported_outward_normal_signs": [-1, 1],
            "exported_element_face_map_sha256": "5" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "hex_sideset_face_ordinals_and_normals_follow_connectivity_reorder"
    ] is False


def test_v18_public_sweep_layer_bias_source_curve_orientation_generation_mismatch():
    row = _with_v18_ordinal_bias_and_quality_identity(summary())
    row["sweep_layer_bias_source_curve_orientation_generation_identity"].update(
        {
            "layer_bias_curve_orientation_generation": "curve-orientation-51",
            "biased_curve_orientations": [-1, 1],
            "biased_interval_counts": [8, 4],
            "applied_bias_factors": [1.5, 1.2],
            "applied_curve_bias_map_sha256": "5" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "biased_sweep_layers_follow_current_source_curve_orientation"
    ] is False


def test_v18_source_exodus_sideset_element_face_topology_generation_mismatch():
    row = _with_v18_ordinal_bias_and_quality_identity(summary())
    row["exodus_sideset_element_face_topology_generation_identity"].update(
        {
            "topology_ordinal_mesh_generation": "mesh-51",
            "topology_ordinal_export_generation": "exodus-51",
            "written_element_face_topology_ordinals": [6, 5],
            "written_element_face_topology_sha256": "5" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "exodus_sideset_ordinals_follow_current_mesh_and_export_topology"
    ] is False


def test_v18_source_high_order_quality_reference_coordinate_generation_mismatch():
    row = _with_v18_ordinal_bias_and_quality_identity(summary())
    row["high_order_quality_reference_coordinate_generation_identity"].update(
        {
            "reference_node_element_order_generation": "element-order-51",
            "quality_element_order_generation": "element-order-51",
            "quality_reference_node_count": 4,
            "quality_jacobian_sampling_rule": "tet4_corner_nodes",
            "quality_reference_coordinates_sha256": "5" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "high_order_quality_uses_current_reference_coordinates_and_order"
    ] is False


def _with_v19_instance_layer_partition_namespace_identity(row):
    row = _with_v18_ordinal_bias_and_quality_identity(row)
    row["periodic_hex_node_pair_transform_instance_generation_identity"] = {
        "mesh_generation": "periodic-hex-mesh-53",
        "node_pair_mesh_generation": "periodic-hex-mesh-53",
        "volume_instance_generation": "volume-instance-53",
        "periodic_transform_volume_instance_generation": "volume-instance-53",
        "node_pair_volume_instance_generation": "volume-instance-53",
        "source_node_ids": [11, 12, 13],
        "target_node_ids": [21, 22, 23],
        "paired_source_node_ids": [11, 12, 13],
        "paired_target_node_ids": [21, 22, 23],
        "transform_translation_m": [0.1, 0.0, 0.0],
        "paired_transform_translation_m": [0.1, 0.0, 0.0],
        "node_pair_transform_sha256": "5" * 64,
        "applied_node_pair_transform_sha256": "5" * 64,
    }
    row["hex_boundary_layer_thickness_surface_normal_generation_identity"] = {
        "geometry_generation": "healed-geometry-53",
        "surface_normal_geometry_generation": "healed-geometry-53",
        "boundary_layer_geometry_generation": "healed-geometry-53",
        "boundary_layer_generation": "boundary-layer-53",
        "thickness_boundary_layer_generation": "boundary-layer-53",
        "collapse_direction_boundary_layer_generation": "boundary-layer-53",
        "surface_ids": [31, 32],
        "surface_normal_signs": [1, -1],
        "applied_surface_ids": [31, 32],
        "applied_collapse_direction_signs": [1, -1],
        "layer_thickness_m": [0.001, 0.0015],
        "applied_layer_thickness_m": [0.001, 0.0015],
        "surface_layer_map_sha256": "6" * 64,
        "applied_surface_layer_map_sha256": "6" * 64,
    }
    row["partition_ghost_element_owner_shared_node_map_identity"] = {
        "partition_generation": "partition-53",
        "ghost_owner_partition_generation": "partition-53",
        "shared_node_partition_generation": "partition-53",
        "partition_ids": [0, 1],
        "element_ids": [101, 102],
        "element_owner_partition_ids": [0, 1],
        "ghost_element_ids": [102, 101],
        "ghost_owner_partition_ids": [1, 0],
        "shared_node_ids": [41, 42],
        "shared_node_partition_pairs": [[0, 1], [0, 1]],
        "partition_ownership_sha256": "7" * 64,
        "exported_partition_ownership_sha256": "7" * 64,
    }
    row["exodus_block_id_namespace_qa_record_mesh_generation_identity"] = {
        "mesh_generation": "mesh-53",
        "block_namespace_mesh_generation": "mesh-53",
        "qa_record_mesh_generation": "mesh-53",
        "exodus_export_generation": "exodus-53",
        "writer_export_generation": "exodus-53",
        "block_ids": [10, 20],
        "block_names": ["rotor", "stator"],
        "written_block_ids": [10, 20],
        "written_block_names": ["rotor", "stator"],
        "qa_record": ["radia-mcp", "v19", "2026-07-17", "headless"],
        "written_qa_record": ["radia-mcp", "v19", "2026-07-17", "headless"],
        "block_namespace_sha256": "8" * 64,
        "written_block_namespace_sha256": "8" * 64,
    }
    return row


def test_v19_positive_instance_layer_partition_namespace_identity():
    row = _with_v19_instance_layer_partition_namespace_identity(summary())
    assert json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))["status"] == "ok"
    assert json.loads(cubit_mixed_transition_source_gate(row))["status"] == "ok"


def test_v19_public_periodic_hex_node_pair_transform_instance_generation_mismatch():
    row = _with_v19_instance_layer_partition_namespace_identity(summary())
    row["periodic_hex_node_pair_transform_instance_generation_identity"].update(
        {
            "periodic_transform_volume_instance_generation": "volume-instance-52",
            "node_pair_volume_instance_generation": "volume-instance-52",
            "paired_target_node_ids": [22, 21, 23],
            "paired_transform_translation_m": [0.0, 0.1, 0.0],
            "applied_node_pair_transform_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "periodic_hex_pairs_follow_current_volume_instance_transform"
    ] is False


def test_v19_public_hex_boundary_layer_thickness_surface_normal_generation_mismatch():
    row = _with_v19_instance_layer_partition_namespace_identity(summary())
    row["hex_boundary_layer_thickness_surface_normal_generation_identity"].update(
        {
            "surface_normal_geometry_generation": "healed-geometry-52",
            "thickness_boundary_layer_generation": "boundary-layer-52",
            "applied_collapse_direction_signs": [-1, 1],
            "applied_layer_thickness_m": [0.0015, 0.001],
            "applied_surface_layer_map_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "hex_boundary_layers_follow_current_healed_surface_normals"
    ] is False


def test_v19_source_partition_ghost_element_owner_shared_node_map_mismatch():
    row = _with_v19_instance_layer_partition_namespace_identity(summary())
    row["partition_ghost_element_owner_shared_node_map_identity"].update(
        {
            "ghost_owner_partition_generation": "partition-52",
            "shared_node_partition_generation": "partition-52",
            "ghost_owner_partition_ids": [0, 1],
            "shared_node_partition_pairs": [[1, 0], [1, 0]],
            "exported_partition_ownership_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "partition_ghosts_and_shared_nodes_use_current_owner_map"
    ] is False


def test_v19_source_exodus_block_id_namespace_qa_record_mesh_generation_mismatch():
    row = _with_v19_instance_layer_partition_namespace_identity(summary())
    row["exodus_block_id_namespace_qa_record_mesh_generation_identity"].update(
        {
            "block_namespace_mesh_generation": "mesh-52",
            "qa_record_mesh_generation": "mesh-52",
            "written_block_ids": [20, 10],
            "written_qa_record": ["radia-mcp", "v18", "2026-07-16", "headless"],
            "written_block_namespace_sha256": "9" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "exodus_blocks_and_qa_use_current_mesh_namespace"
    ] is False


def _with_v20_jacobian_interface_journal_vol_identity(row):
    row = _with_v19_instance_layer_partition_namespace_identity(row)
    row["high_order_hex_jacobian_node_order_coordinate_scale_identity"] = {
        "mesh_generation": "hex-mesh-62",
        "curving_mesh_generation": "hex-mesh-62",
        "jacobian_mesh_generation": "hex-mesh-62",
        "node_order_generation": "hex-order-62",
        "jacobian_node_order_generation": "hex-order-62",
        "coordinate_scale_generation": "coordinate-scale-62",
        "jacobian_coordinate_scale_generation": "coordinate-scale-62",
        "element_order": 2,
        "jacobian_element_order": 2,
        "corner_node_ids": [1, 2, 3, 4, 5, 6, 7, 8],
        "jacobian_corner_node_ids": [1, 2, 3, 4, 5, 6, 7, 8],
        "reference_corner_coordinates": [
            [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
        ],
        "jacobian_reference_corner_coordinates": [
            [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
        ],
        "coordinate_scale_m": 0.001,
        "jacobian_coordinate_scale_m": 0.001,
        "minimum_scaled_jacobian": 0.72,
        "jacobian_table_sha256": "a" * 64,
        "evaluated_jacobian_table_sha256": "a" * 64,
    }
    row["tet_hex_pyramid_interface_face_orientation_conformity_identity"] = {
        "interface_generation": "transition-interface-62",
        "tet_mesh_generation": "tet-mesh-62",
        "interface_tet_mesh_generation": "tet-mesh-62",
        "hex_mesh_generation": "hex-mesh-62",
        "interface_hex_mesh_generation": "hex-mesh-62",
        "pyramid_transition_generation": "pyramid-transition-62",
        "interface_pyramid_transition_generation": "pyramid-transition-62",
        "quad_face_node_ids": [[11, 12, 13, 14]],
        "interface_quad_face_node_ids": [[11, 12, 13, 14]],
        "pyramid_base_node_ids": [[11, 12, 13, 14]],
        "interface_pyramid_base_node_ids": [[11, 12, 13, 14]],
        "pyramid_apex_node_ids": [21],
        "interface_pyramid_apex_node_ids": [21],
        "face_orientation_signs": [1],
        "interface_face_orientation_signs": [1],
        "interface_conformity_sha256": "b" * 64,
        "exported_interface_conformity_sha256": "b" * 64,
    }
    row["journal_transaction_undo_entity_id_reuse_generation_identity"] = {
        "journal_generation": "journal-62",
        "transaction_journal_generation": "journal-62",
        "entity_table_journal_generation": "journal-62",
        "group_table_journal_generation": "journal-62",
        "transaction_id": "transaction-62",
        "replay_transaction_id": "transaction-62",
        "reset_epoch": 7,
        "replay_reset_epoch": 7,
        "undo_depth": 1,
        "replay_undo_depth": 1,
        "created_entity_ids": [101, 102],
        "replay_created_entity_ids": [101, 102],
        "group_entity_ids": [101, 102],
        "replay_group_entity_ids": [101, 102],
        "transaction_entity_table_sha256": "c" * 64,
        "replay_transaction_entity_table_sha256": "c" * 64,
    }
    row["netgen_vol_element_block_order_curving_generation_identity"] = {
        "mesh_generation": "hybrid-mesh-62",
        "writer_mesh_generation": "hybrid-mesh-62",
        "export_generation": "netgen-export-62",
        "writer_export_generation": "netgen-export-62",
        "curving_generation": "curving-62",
        "writer_curving_generation": "curving-62",
        "element_block_ids": [10, 20, 30],
        "element_block_types": ["tet", "hex", "pyramid"],
        "writer_element_block_types": ["tet", "hex", "pyramid"],
        "element_orders": [2, 2, 1],
        "writer_element_orders": [2, 2, 1],
        "curving_node_counts": [10, 20, 5],
        "writer_curving_node_counts": [10, 20, 5],
        "element_block_table_sha256": "d" * 64,
        "writer_element_block_table_sha256": "d" * 64,
    }
    return row


def test_v20_positive_jacobian_interface_journal_vol_identity():
    row = _with_v20_jacobian_interface_journal_vol_identity(summary())
    assert json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))["status"] == "ok"
    assert json.loads(cubit_mixed_transition_source_gate(row))["status"] == "ok"


def test_v20_public_high_order_hex_jacobian_node_order_coordinate_scale_generation_mismatch():
    row = _with_v20_jacobian_interface_journal_vol_identity(summary())
    row["high_order_hex_jacobian_node_order_coordinate_scale_identity"].update(
        {
            "jacobian_mesh_generation": "hex-mesh-61",
            "jacobian_node_order_generation": "hex-order-61",
            "jacobian_coordinate_scale_generation": "coordinate-scale-61",
            "jacobian_corner_node_ids": [1, 4, 3, 2, 5, 8, 7, 6],
            "jacobian_coordinate_scale_m": 1.0,
            "evaluated_jacobian_table_sha256": "f" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "high_order_hex_jacobians_use_current_node_order_and_coordinate_scale"
    ] is False


def test_v20_public_tet_hex_pyramid_interface_face_orientation_conformity_generation_mismatch():
    row = _with_v20_jacobian_interface_journal_vol_identity(summary())
    row["tet_hex_pyramid_interface_face_orientation_conformity_identity"].update(
        {
            "interface_tet_mesh_generation": "tet-mesh-61",
            "interface_pyramid_transition_generation": "pyramid-transition-61",
            "interface_quad_face_node_ids": [[11, 14, 13, 12]],
            "interface_pyramid_base_node_ids": [[11, 14, 13, 12]],
            "interface_face_orientation_signs": [-1],
            "exported_interface_conformity_sha256": "f" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "tet_hex_pyramid_interface_uses_current_face_orientation_and_conformity"
    ] is False


def test_v20_source_journal_transaction_undo_entity_id_reuse_generation_mismatch():
    row = _with_v20_jacobian_interface_journal_vol_identity(summary())
    row["journal_transaction_undo_entity_id_reuse_generation_identity"].update(
        {
            "entity_table_journal_generation": "journal-61",
            "group_table_journal_generation": "journal-61",
            "replay_transaction_id": "transaction-61",
            "replay_reset_epoch": 6,
            "replay_created_entity_ids": [101, 103],
            "replay_group_entity_ids": [101, 103],
            "replay_transaction_entity_table_sha256": "f" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "journal_replay_uses_current_transaction_undo_and_entity_ids"
    ] is False


def test_v20_source_netgen_vol_element_block_order_curving_generation_mismatch():
    row = _with_v20_jacobian_interface_journal_vol_identity(summary())
    row["netgen_vol_element_block_order_curving_generation_identity"].update(
        {
            "writer_mesh_generation": "hybrid-mesh-61",
            "writer_curving_generation": "curving-61",
            "writer_element_block_types": ["tet", "pyramid", "hex"],
            "writer_element_orders": [1, 1, 2],
            "writer_curving_node_counts": [4, 5, 20],
            "writer_element_block_table_sha256": "f" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "netgen_export_uses_current_block_order_and_curving_generation"
    ] is False


def _with_v21_sweep_quality_group_aprepro_identity(row):
    row = _with_v20_jacobian_interface_journal_vol_identity(row)
    row["hex_sweep_source_target_face_vertex_twist_generation_identity"] = {
        "sweep_generation": "hex-sweep-71",
        "source_face_sweep_generation": "hex-sweep-71",
        "target_face_sweep_generation": "hex-sweep-71",
        "vertex_map_sweep_generation": "hex-sweep-71",
        "twist_path_sweep_generation": "hex-sweep-71",
        "source_face_id": 11,
        "mapped_source_face_id": 11,
        "target_face_id": 21,
        "mapped_target_face_id": 21,
        "source_vertex_ids": [101, 102, 103, 104],
        "mapped_source_vertex_ids": [101, 102, 103, 104],
        "target_vertex_ids": [201, 202, 203, 204],
        "mapped_target_vertex_ids": [201, 202, 203, 204],
        "twist_path_vertex_ids": [101, 201, 102, 202, 103, 203, 104, 204],
        "mapped_twist_path_vertex_ids": [101, 201, 102, 202, 103, 203, 104, 204],
        "face_vertex_map_sha256": "1" * 64,
        "applied_face_vertex_map_sha256": "1" * 64,
        "twist_path_sha256": "2" * 64,
        "applied_twist_path_sha256": "2" * 64,
    }
    row["quality_histogram_metric_element_set_unit_generation_identity"] = {
        "mesh_generation": "quality-mesh-71",
        "metric_mesh_generation": "quality-mesh-71",
        "element_set_mesh_generation": "quality-mesh-71",
        "coordinate_unit_mesh_generation": "quality-mesh-71",
        "metric_name": "scaled_jacobian",
        "evaluated_metric_name": "scaled_jacobian",
        "coordinate_unit": "m",
        "evaluated_coordinate_unit": "m",
        "element_ids": [301, 302, 303, 304],
        "evaluated_element_ids": [301, 302, 303, 304],
        "metric_values": [0.62, 0.74, 0.83, 0.91],
        "evaluated_metric_values": [0.62, 0.74, 0.83, 0.91],
        "histogram_bin_edges": [0.0, 0.5, 0.75, 1.0],
        "histogram_counts": [0, 2, 2],
        "evaluated_histogram_counts": [0, 2, 2],
        "quality_table_sha256": "3" * 64,
        "evaluated_quality_table_sha256": "3" * 64,
    }
    row["block_sideset_group_entity_merge_renumber_generation_identity"] = {
        "topology_generation": "topology-71",
        "block_topology_generation": "topology-71",
        "sideset_topology_generation": "topology-71",
        "group_topology_generation": "topology-71",
        "renumber_topology_generation": "topology-71",
        "merge_transaction_generation": "merge-71",
        "group_merge_transaction_generation": "merge-71",
        "block_ids": [10, 20],
        "exported_block_ids": [10, 20],
        "block_entity_ids": [[401, 402], [403]],
        "exported_block_entity_ids": [[401, 402], [403]],
        "sideset_ids": [30, 40],
        "exported_sideset_ids": [30, 40],
        "sideset_entity_ids": [[501, 502], [503, 504]],
        "exported_sideset_entity_ids": [[501, 502], [503, 504]],
        "group_entity_ids": [401, 402, 403, 501, 502, 503, 504],
        "exported_group_entity_ids": [401, 402, 403, 501, 502, 503, 504],
        "ownership_table_sha256": "4" * 64,
        "exported_ownership_table_sha256": "4" * 64,
    }
    row[
        "aprepro_include_variable_expansion_working_directory_generation_identity"
    ] = {
        "journal_transaction_generation": "journal-71",
        "variable_table_transaction_generation": "journal-71",
        "include_expansion_transaction_generation": "journal-71",
        "working_directory_transaction_generation": "journal-71",
        "working_directory": "model/input",
        "replay_working_directory": "model/input",
        "variable_names": ["radius", "height", "intervals"],
        "expanded_variable_names": ["radius", "height", "intervals"],
        "variable_values": [0.025, 0.08, 12.0],
        "expanded_variable_values": [0.025, 0.08, 12.0],
        "include_paths": ["common/units.inc", "mesh/sweep.inc"],
        "expanded_include_paths": ["common/units.inc", "mesh/sweep.inc"],
        "variable_table_sha256": "5" * 64,
        "expanded_variable_table_sha256": "5" * 64,
        "include_tree_sha256": "6" * 64,
        "expanded_include_tree_sha256": "6" * 64,
    }
    return row


def test_v21_positive_sweep_quality_group_aprepro_identity():
    row = _with_v21_sweep_quality_group_aprepro_identity(summary())
    assert json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))["status"] == "ok"
    assert json.loads(cubit_mixed_transition_source_gate(row))["status"] == "ok"


def test_v21_public_hex_sweep_source_target_face_vertex_twist_generation_mismatch():
    row = _with_v21_sweep_quality_group_aprepro_identity(summary())
    row["hex_sweep_source_target_face_vertex_twist_generation_identity"].update(
        {
            "target_face_sweep_generation": "hex-sweep-70",
            "vertex_map_sweep_generation": "hex-sweep-69",
            "twist_path_sweep_generation": "hex-sweep-68",
            "mapped_target_face_id": 22,
            "mapped_source_vertex_ids": [104, 103, 102, 101],
            "mapped_twist_path_vertex_ids": [101, 202, 102, 203, 103, 204, 104, 201],
            "applied_face_vertex_map_sha256": "a" * 64,
            "applied_twist_path_sha256": "b" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "hex_sweep_uses_current_source_target_vertex_map_and_twist_path"
    ] is False


def test_v21_public_quality_histogram_metric_element_set_unit_generation_mismatch():
    row = _with_v21_sweep_quality_group_aprepro_identity(summary())
    row["quality_histogram_metric_element_set_unit_generation_identity"].update(
        {
            "metric_mesh_generation": "quality-mesh-70",
            "element_set_mesh_generation": "quality-mesh-69",
            "coordinate_unit_mesh_generation": "quality-mesh-68",
            "evaluated_metric_name": "aspect_ratio",
            "evaluated_coordinate_unit": "mm",
            "evaluated_element_ids": [301, 303, 305],
            "evaluated_metric_values": [1.2, 2.1, 3.4],
            "evaluated_histogram_counts": [1, 1, 1],
            "evaluated_quality_table_sha256": "c" * 64,
        }
    )
    result = json.loads(cubit_conformal_hex_pyramid_tet_interface_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "quality_histogram_uses_current_metric_element_set_and_units"
    ] is False


def test_v21_source_block_sideset_group_entity_merge_renumber_generation_mismatch():
    row = _with_v21_sweep_quality_group_aprepro_identity(summary())
    row["block_sideset_group_entity_merge_renumber_generation_identity"].update(
        {
            "group_topology_generation": "topology-70",
            "renumber_topology_generation": "topology-69",
            "group_merge_transaction_generation": "merge-70",
            "exported_block_entity_ids": [[401, 405], [403]],
            "exported_sideset_entity_ids": [[501], [503, 505]],
            "exported_group_entity_ids": [401, 403, 405, 501, 503, 505],
            "exported_ownership_table_sha256": "d" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "block_sideset_groups_use_current_merge_and_renumber_generation"
    ] is False


def test_v21_source_aprepro_include_variable_expansion_working_directory_generation_mismatch():
    row = _with_v21_sweep_quality_group_aprepro_identity(summary())
    row[
        "aprepro_include_variable_expansion_working_directory_generation_identity"
    ].update(
        {
            "variable_table_transaction_generation": "journal-70",
            "include_expansion_transaction_generation": "journal-69",
            "working_directory_transaction_generation": "journal-68",
            "replay_working_directory": "archive/input",
            "expanded_variable_values": [25.0, 80.0, 10.0],
            "expanded_include_paths": ["old/units.inc", "mesh/tet.inc"],
            "expanded_variable_table_sha256": "e" * 64,
            "expanded_include_tree_sha256": "f" * 64,
        }
    )
    result = json.loads(cubit_mixed_transition_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "aprepro_replay_uses_current_variables_includes_and_working_directory"
    ] is False
