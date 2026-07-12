import copy
import json

from radia_mcp.cubit.server import (
    cubit_embedded_pipe_source_recovery_gate,
    cubit_embedded_region_mixed_transition_gate,
)


def live_summary() -> dict:
    checks = {
        "source_native_forum_workflow_recorded": True,
        "source_hex_stage_completed": True,
        "version_drift_unmeshed_set_recorded": True,
        "all_source_volumes_recovered": True,
        "pipe_block_survived_decomposition": True,
        "soil_block_survived_decomposition": True,
        "hex_dominant_mesh_present": True,
        "tet_recovery_mesh_present": True,
        "pyramid_transition_present": True,
        "unsupported_wedge_absent": True,
        "shared_pipe_soil_interfaces_meshed": True,
        "recovery_transition_surfaces_meshed": True,
        "positive_scaled_jacobian": True,
        "brick_union_volume_conserved": True,
        "swept_pipe_volume_matches_path_area": True,
        "gmsh_v41_export_created": True,
    }
    quad_interface = {
        "surface_id": 10,
        "adjacent_volumes": [1, 2],
        "face_count": 64,
        "quad_count": 64,
        "tri_count": 0,
    }
    return {
        "schema": "coreform-cubit.embedded-pipe-mixed-transition-live.v1",
        "pass": True,
        "version": "2025.12",
        "source_kind": "source_native_coreform_forum_journal_with_version_recovery",
        "source_journal": "embedded_pipe.jou",
        "source_sha256": "a" * 64,
        "source_url": "https://forum.coreform.com/t/meshing-an-embedded-pipe/2759/1",
        "execution_mode": "headless_combined_journal_then_python_inventory",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "recovery": {
            "unmeshed_soil_before": [8, 12, 17, 21, 26, 28, 31, 33],
            "fallback_scheme": "tetmesh",
            "fallback_size": 0.05,
            "pyramid_transition_required": True,
        },
        "volume_count": 33,
        "unmeshed_volumes_after": [],
        "preexisting_hex_count": 409980,
        "element_counts": {
            "hex": 409980,
            "tet": 132113,
            "pyramid": 6814,
            "wedge": 0,
        },
        "quality": {
            "hex": {"scaled_jacobian": {"count": 409980, "min": 0.7368}},
            "tet": {"scaled_jacobian": {"count": 132113, "min": 0.1545}},
            "pyramid": {"scaled_jacobian": {"count": 6814, "min": 0.4883}},
        },
        "pipe_soil_interfaces": [quad_interface],
        "recovery_transition_surfaces": [
            {**quad_interface, "surface_id": 20, "adjacent_volumes": [8, 9]}
        ],
        "geometry": {
            "cad_volume": 7.0000054158,
            "analytic_brick_volume": 7.0,
            "pipe_volume": 0.210358990434947,
            "analytic_pipe_volume": 0.210358990434947,
        },
        "gmsh_export": {
            "bytes": 43209355,
            "sha256": "b" * 64,
            "header": {
                "version": "4.1",
                "file_type": 0,
                "has_entities_section": True,
                "has_nodes_section": True,
                "has_elements_section": True,
            },
        },
        "timing": {
            "source_hex_journal_s": 5.5,
            "tet_pyramid_recovery_s": 1.5,
            "gmsh_export_and_parse_s": 3.5,
            "total_s": 11.3,
        },
        "process": {
            "exit_code": 24,
            "error_categories": [
                "headless_startup_diagnostics",
                "acis_webcut_version_drift",
                "source_unmeshed_volumes_recovered",
                "session_error_summary",
            ],
            "unexpected_error_count": 0,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
        "checks": checks,
    }


def test_public_gate_accepts_verified_hex_tet_pyramid_transition():
    result = json.loads(cubit_embedded_region_mixed_transition_gate(live_summary()))
    assert result["status"] == "ok"
    assert result["quality_minima"]["tet"] == 0.1545
    assert result["checks"]["gmsh_ascii_v41_handoff_complete"] is True
    assert result["checks"]["hex_tet_transition_is_conformal_quad_mesh"] is True


def test_public_gate_rejects_missing_transition_and_quality_count_drift():
    row = live_summary()
    row["element_counts"]["pyramid"] = 0
    row["quality"]["tet"]["scaled_jacobian"]["count"] -= 1
    result = json.loads(cubit_embedded_region_mixed_transition_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["hex_led_mixed_topology"] is False
    assert result["checks"]["quality_count_matches_topology"] is False


def test_public_gate_rejects_unmeshed_volume_and_binary_or_old_gmsh():
    row = live_summary()
    row["unmeshed_volumes_after"] = [31]
    row["gmsh_export"]["header"]["version"] = "2.2"
    row["gmsh_export"]["header"]["file_type"] = 1
    result = json.loads(cubit_embedded_region_mixed_transition_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_cad_volumes_meshed"] is False
    assert result["checks"]["gmsh_ascii_v41_handoff_complete"] is False


def test_source_gate_accepts_semantically_explained_nonzero_exit():
    result = json.loads(cubit_embedded_pipe_source_recovery_gate(live_summary()))
    assert result["status"] == "ok"
    assert result["process_exit_code"] == 24
    assert result["checks"]["nonzero_exit_semantically_explained"] is True
    assert result["checks"]["independent_mixed_transition_gate_passed"] is True


def test_source_gate_rejects_exit_code_only_allowlisting_and_unexpected_error():
    row = live_summary()
    row["process"]["error_categories"] = ["session_error_summary"]
    row["process"]["unexpected_error_count"] = 1
    result = json.loads(cubit_embedded_pipe_source_recovery_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["nonzero_exit_semantically_explained"] is False


def test_source_gate_rejects_changed_version_drift_set_and_stale_artifact():
    row = live_summary()
    row["recovery"]["unmeshed_soil_before"] = [8, 12]
    row["process"]["result_artifact_fresh"] = False
    result = json.loads(cubit_embedded_pipe_source_recovery_gate(row))
    assert result["status"] == "needs_attention"
    assert result["checks"]["version_drift_volume_set_matches"] is False
    assert result["checks"]["fresh_result_and_no_owned_process_leak"] is False


def test_server_returns_invalid_input_for_empty_quality_family():
    row = live_summary()
    row["quality"]["pyramid"] = {}
    result = json.loads(cubit_embedded_region_mixed_transition_gate(row))
    assert result["status"] == "invalid_input"
    assert "quality.pyramid.scaled_jacobian" in result["error"]
