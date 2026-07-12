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
