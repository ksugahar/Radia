import copy
import json

from radia_mcp.build123d.server import (
    build123d_dual_api_prismatic_pattern_gate,
    build123d_dual_api_source_replay_gate,
)


def good():
    native_row = {
        "volume": 42462.8633886947,
        "area": 88463.3001628255,
        "bbox_extent": [35.0, 1000.0, 7.5],
        "solid_count": 1,
        "face_count": 178,
        "edge_count": 528,
        "vertex_count": 352,
        "is_valid": True,
        "self_roundtrip": {"is_valid": True},
    }
    native = {
        "expected_volume": 42462.863388694714,
        "records": {
            "builder": copy.deepcopy(native_row),
            "algebra": copy.deepcopy(native_row),
            "algebra_centered": copy.deepcopy(native_row),
        },
    }
    external_row = {
        "volume": 42462.91,
        "bbox_extent": [35.0, 1000.0, 7.5],
        "volume_count": 1,
        "surface_count": 178,
        "curve_count": 528,
        "vertex_count": 352,
    }
    external = {
        "records": {
            "builder": copy.deepcopy(external_row),
            "algebra": copy.deepcopy(external_row),
            "algebra_centered": copy.deepcopy(external_row),
        }
    }
    return {"native": native, "external": external}


def source_good():
    return {
        "upstream_commit": "a" * 40,
        "sources": {
            "builder": {"sha256": "b" * 64, "preserved": True},
            "algebra": {"sha256": "c" * 64, "preserved": True},
        },
        "official_assertion_reproduced": True,
        "viewer_suppressed": True,
        "source_replay_mode": "runpy_with_viewer_stub",
        "artifacts": {
            "builder": {"sha256": "d" * 64},
            "algebra": {"sha256": "e" * 64},
            "algebra_centered": {"sha256": "f" * 64, "derived_control": True},
        },
        "external_execution_mode": "python_api_headless",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "external_process_exit_code": 3,
        "startup_diagnostics": [
            "ERROR: Could not open file: C:/x/plugins",
            "ERROR: Could not open file: -commandplugindir",
            "ERROR: Could not open file: -nojournal",
        ],
        "script_error_lines": [],
        "result_artifact_fresh": True,
        "owned_processes_remaining": 0,
        "public_gate_status": "ok",
    }


def test_accepts_native_exact_and_external_ppm_contract():
    result = json.loads(build123d_dual_api_prismatic_pattern_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["checks"]["native_dual_api_topology_match"] is True


def test_rejects_external_volume_and_topology_drift():
    row = good()
    row["external"]["records"]["algebra"]["volume"] *= 1.00001
    row["external"]["records"]["algebra"]["surface_count"] -= 1
    result = json.loads(build123d_dual_api_prismatic_pattern_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["external_topology_matches_across_exports"] is False


def test_source_replay_gate_accepts_live_contract_and_rejects_rewritten_source():
    assert json.loads(build123d_dual_api_source_replay_gate(json.dumps(source_good())))["status"] == "ok"
    bad = source_good()
    bad["sources"]["algebra"]["preserved"] = False
    result = json.loads(build123d_dual_api_source_replay_gate(json.dumps(bad)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_files_preserved"] is False
