import copy
import json

from radia_mcp.cubit.server import (
    cubit_helical_conductor_source_gate,
    cubit_region_owned_mixed_mesh_gate,
)


def summary() -> dict:
    counts = {"hex": 120, "tet": 900, "pyramid": 40, "wedge": 0}
    interfaces = [
        {
            "surface_id": volume,
            "adjacent_volumes": [volume, 6],
            "face_count": 24,
            "quad_count": 24,
            "tri_count": 0,
            "area": 1.0,
        }
        for volume in range(1, 6)
    ]
    return {
        "source_kind": "source_native_local_helical_conductor_journal_adapted_to_hex_tet_pyramid",
        "source_journal": "helical_conductor.jou",
        "source_sha256": "a" * 64,
        "execution_mode": "headless_combined_journal_then_python_inventory",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "geometry_parameters": {"a_mm": 100.0, "b_mm": 15.0, "h_mm": 50.0, "turns": 4},
        "conductor_volumes": [1, 2, 3, 4, 5],
        "air_volume": 6,
        "element_counts": counts,
        "per_volume_element_counts": {
            str(volume): {"hex": 24, "tet": 0, "pyramid": 0, "wedge": 0}
            for volume in range(1, 6)
        }
        | {"6": {"hex": 0, "tet": 900, "pyramid": 40, "wedge": 0}},
        "quality": {
            "hex": {"scaled_jacobian": {"count": 120, "min": 0.62}},
            "tet": {"scaled_jacobian": {"count": 900, "min": 0.12}},
            "pyramid": {"scaled_jacobian": {"count": 40, "min": 0.14}},
        },
        "conductor_air_interfaces": interfaces,
        "geometry": {"volume_relative_error": 1.5e-8},
        "gmsh_export": {
            "bytes": 1024,
            "sha256": "b" * 64,
            "header": {
                "version": "4.1",
                "file_type": 0,
                "has_entities_section": True,
                "has_nodes_section": True,
                "has_elements_section": True,
            },
        },
        "gmsh_inventory": {
            "status": "ok",
            "volume_family_counts": {"hex": 120, "tet": 900, "pyramid": 40},
            "connectivity_mismatches": [],
        },
        "timing": {
            "geometry_boolean_s": 1.0,
            "conductor_sweep_hex_s": 4.0,
            "air_tet_pyramid_s": 2.0,
            "gmsh_export_inventory_s": 1.0,
        },
        "process": {
            "exit_code": 4,
            "error_categories": [
                "headless_startup_diagnostics",
                "hpc_tet_attempt_failed",
                "standard_tet_pyramid_fallback_completed",
                "session_error_summary",
            ],
            "unexpected_error_lines": [],
            "hpc_fallback_completed": True,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
            "process_exit_policy": "classified_hpc_fallback_plus_fresh_artifact",
        },
    }


def test_region_gate_accepts_air_dominated_tet_count_with_hex_owned_conductors():
    result = json.loads(cubit_region_owned_mixed_mesh_gate(summary()))
    assert result["status"] == "ok"
    assert result["element_counts"]["tet"] > result["element_counts"]["hex"]
    assert result["checks"]["conductor_regions_are_hex_only"] is True


def test_region_gate_rejects_missing_pyramids_and_interface_owner_drift():
    row = summary()
    row["element_counts"]["pyramid"] = 0
    row["per_volume_element_counts"]["6"]["pyramid"] = 0
    row["quality"]["pyramid"]["scaled_jacobian"]["count"] = 0
    row["conductor_air_interfaces"][0]["adjacent_volumes"] = [1, 7]
    result = json.loads(cubit_region_owned_mixed_mesh_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["air_region_is_tet_with_pyramid_transition"] is False
    assert result["checks"]["one_conformal_quad_interface_per_conductor_region"] is False


def test_region_gate_rejects_quality_count_drift_and_old_gmsh():
    row = summary()
    row["quality"]["tet"]["scaled_jacobian"]["count"] -= 1
    row["gmsh_export"]["header"]["version"] = "2.2"
    row["gmsh_inventory"]["volume_family_counts"]["tet"] -= 1
    result = json.loads(cubit_region_owned_mixed_mesh_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["quality_count_matches_topology"] is False
    assert result["checks"]["gmsh_ascii_v41_handoff_complete"] is False
    assert result["checks"]["parsed_gmsh_topology_matches_cubit"] is False


def test_source_gate_accepts_classified_hpc_fallback_with_independent_gate():
    result = json.loads(cubit_helical_conductor_source_gate(summary()))
    assert result["status"] == "ok"
    assert result["process_exit_code"] == 4
    assert result["checks"]["nonzero_exit_has_classified_completed_fallback"] is True


def test_source_gate_rejects_exit_code_only_allowlist_and_unexpected_error():
    row = summary()
    row["process"]["error_categories"] = ["session_error_summary"]
    row["process"]["unexpected_error_lines"] = ["ERROR: unclassified failure"]
    result = json.loads(cubit_helical_conductor_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["nonzero_exit_has_classified_completed_fallback"] is False


def test_source_gate_rejects_stale_artifact_and_public_mesh_failure():
    row = summary()
    row["process"]["result_artifact_fresh"] = False
    row["quality"]["pyramid"]["scaled_jacobian"]["min"] = 0.01
    result = json.loads(cubit_helical_conductor_source_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["fresh_artifact_and_no_owned_process_leak"] is False
    assert result["checks"]["independent_region_owned_mesh_gate_passed"] is False


def test_server_returns_invalid_input_for_missing_quality_family():
    row = copy.deepcopy(summary())
    row["quality"]["pyramid"] = {}
    result = json.loads(cubit_region_owned_mixed_mesh_gate(row))
    assert result["status"] == "invalid_input"
    assert "quality.pyramid.scaled_jacobian" in result["error"]
