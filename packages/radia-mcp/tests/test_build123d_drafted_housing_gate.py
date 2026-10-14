import copy
import json

from radia_mcp.build123d.server import (
    build123d_drafted_housing_cross_kernel_gate,
    build123d_drafted_housing_source_replay_gate,
)


def cross_kernel_summary():
    native = {
        "volume_mm3": 46882.6848405,
        "area_mm2": 14647.5747291,
        "bbox_extent_mm": [98.35, 56.85, 26.0],
        "solid_count": 1,
        "shell_count": 1,
        "face_count": 37,
        "edge_count": 85,
        "vertex_count": 51,
        "is_valid": True,
    }
    run = {
        "native": native,
        "tessellated": {"volume_mm3": 46870.53, "area_mm2": 14647.16},
        "self_roundtrip": {"volume_mm3": 46881.337, "area_mm2": 14647.491, "solid_count": 1, "is_valid": True},
    }
    rows = []
    for run_index in (1, 2):
        for mode in ("noheal", "heal"):
            rows.append({"run_index": run_index, "import_mode": mode, "volume_count": 1, "positive_volume_count": 1, "surface_count": 37, "curve_count": 85, "vertex_count": 51, "total_volume_mm3": 46883.0776, "total_area_mm2": 14646.7964})
    return {
        "runs": [dict(copy.deepcopy(run), run_index=1), dict(copy.deepcopy(run), run_index=2)],
        "external_rows": rows,
        "mesh_evidence": {"tet_count": 137697, "node_count": 29976, "connectivity_orders": [4], "min_scaled_jacobian": 0.1075, "cad_volume_relative_error": 0.00128},
        "gmsh_inventory": {"status": "ok", "mesh_format": "4.1", "binary": False, "node_count": 29976, "element_family_counts": {"tet": 137697}, "connectivity_mismatches": []},
    }


def source_summary():
    return {
        "source": {"kind": "upstream_native_example_tag_v0.10.0", "version": "0.10.0", "sha256": "a" * 64, "commit": "b" * 40, "preserved": True, "display_stubbed_only": True},
        "source_contract": {"operations": ["draft", "fillet", "counterbore", "through_holes"], "parameters": {"draft_angle_deg": 4.0, "mounting_hole_count": 2, "counterbore_depth_mm": 14.0}},
        "native_replay_count": 2,
        "step_sha256": ["c" * 64, "d" * 64],
        "external_execution": {"mode": "python_api_headless", "headless_flags": ["-nographics", "-batch"], "gui_daemon_enabled": False, "result_artifact_fresh": True, "owned_processes_remaining": 0},
        "gmsh_companions": {"msh": True, "geo": True, "geo_opt": True, "msh_opt": True},
        "cross_kernel_gate_status": "ok",
        "solver_ready": True,
        "timing_breakdown_s": {"native_brep_and_step": 8.0, "external_import_and_mesh": 22.0, "mcp_protocol": 5.0, "verification_and_commit": 10.0},
    }


def test_cross_kernel_gate_accepts_bounded_multi_kernel_spread_and_tet_mesh():
    result = json.loads(build123d_drafted_housing_cross_kernel_gate(json.dumps(cross_kernel_summary())))
    assert result["status"] == "ok"
    assert result["solver_ready"] is True
    assert result["checks"]["gmsh_v41_contains_the_complete_tet_block"] is True


def test_cross_kernel_gate_rejects_missing_tets_and_false_positive_volume():
    row = cross_kernel_summary()
    row["external_rows"][0]["positive_volume_count"] = 0
    row["mesh_evidence"]["tet_count"] = 0
    row["gmsh_inventory"]["element_family_counts"]["tet"] = 0
    result = json.loads(build123d_drafted_housing_cross_kernel_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["solver_ready"] is False
    assert result["checks"]["external_imports_are_positive_single_solids"] is False
    assert result["checks"]["tet_mesh_is_positive_first_order"] is False


def test_source_replay_gate_accepts_tagged_headless_companion_contract():
    result = json.loads(build123d_drafted_housing_source_replay_gate(json.dumps(source_summary())))
    assert result["status"] == "ok"
    assert result["checks"]["drafted_perforated_source_contract_recorded"] is True


def test_source_replay_gate_rejects_geometry_edit_gui_and_missing_geo_opt():
    row = source_summary()
    row["source"]["preserved"] = False
    row["external_execution"]["gui_daemon_enabled"] = True
    row["gmsh_companions"]["geo_opt"] = False
    row["cross_kernel_gate_status"] = "needs_attention"
    result = json.loads(build123d_drafted_housing_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_preserved_except_display_stub"] is False
    assert result["checks"]["headless_external_cad_execution"] is False
    assert result["checks"]["gmsh_mesh_and_launch_companions_recorded"] is False
    assert result["checks"]["cross_kernel_gate_accepts_solver_handoff"] is False
