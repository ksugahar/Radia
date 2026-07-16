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
