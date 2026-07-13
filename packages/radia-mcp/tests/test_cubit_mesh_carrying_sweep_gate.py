import copy
import json

from radia_mcp.cubit.mesh_carrying_sweep_gate import (
    cubit_mesh_carrying_straight_sweep_gate,
    cubit_mesh_carrying_straight_sweep_source_replay_gate,
)
from radia_mcp.cubit.server import (
    cubit_mesh_carrying_straight_sweep_gate as mcp_public_gate,
    cubit_mesh_carrying_straight_sweep_source_replay_gate as mcp_source_gate,
)


def _public_summary():
    return {
        "sweep_mode": "vector",
        "command": "sweep surface 1 vector 0 0 1 distance 6 include_mesh",
        "source_quad_count": 24,
        "source_node_count": 35,
        "source_xy_count": 35,
        "element_counts": {"hex": 144, "tet": 0, "pyramid": 0, "wedge": 0},
        "node_count": 245,
        "z_levels": [0, 1, 2, 3, 4, 5, 6],
        "sweep_interval_count": 6,
        "xy_column_count": 35,
        "xy_column_depths": [7],
        "complete_xy_column_count": 35,
        "cad_total_volume": 576.0,
        "expected_volume": 576.0,
        "quality": {"scaled_jacobian": {"min": 1.0}},
        "gmsh_export": {
            "mesh_format": "4.1",
            "binary": False,
            "node_count": 245,
            "hex_count": 144,
            "other_volume_count": 0,
        },
    }


def _source_summary():
    public = cubit_mesh_carrying_straight_sweep_gate(_public_summary())
    return {
        "source_kind": "installed-official-help-command-with-synthetic-replay",
        "source_sha256": "a" * 64,
        "source_contract": (
            "include_mesh carries an already meshed planar source into a volume"
        ),
        "execution_mode": "combined_journal_headless",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "process": {
            "acceptable": True,
            "result_artifact_fresh": True,
            "unexpected_error_lines": [],
        },
        "api_entity_contract": {
            "surface_mesh_entity": "face",
            "quad_connectivity_size": 4,
            "quad_alias_count": 0,
        },
        "public_gate": public,
        "negative_control": {
            "command": "sweep surface 1 vector 0 0 1 distance 6",
            "volume_ids": [1],
            "cad_total_volume": 576.0,
            "element_counts": {"hex": 0, "tet": 0, "pyramid": 0, "wedge": 0},
        },
        "timing_breakdown_s": {
            "headless_replay": 1.0,
            "process_classification": 0.1,
            "independent_export_parse": 0.2,
            "artifact_finalization": 0.1,
        },
    }


def test_accepts_live_straight_mesh_carrying_sweep():
    result = cubit_mesh_carrying_straight_sweep_gate(_public_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["expected_hex_count"] == 144
    assert json.loads(mcp_public_gate(_public_summary()))["status"] == "ok"


def test_rejects_missing_include_mesh_even_when_counts_look_valid():
    bad = _public_summary()
    bad["command"] = bad["command"].replace(" include_mesh", "")
    result = cubit_mesh_carrying_straight_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["command_carries_existing_mesh"] is False


def test_rejects_topology_quality_volume_and_export_drift():
    bad = copy.deepcopy(_public_summary())
    bad["element_counts"]["hex"] = 143
    bad["node_count"] = 244
    bad["xy_column_depths"] = [6, 7]
    bad["quality"]["scaled_jacobian"]["min"] = 0.1
    bad["cad_total_volume"] = 500.0
    bad["gmsh_export"]["mesh_format"] = "2.2"
    result = cubit_mesh_carrying_straight_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) >= 6


def test_accepts_source_replay_and_no_mesh_control():
    summary = _source_summary()
    result = cubit_mesh_carrying_straight_sweep_source_replay_gate(summary)
    assert result["status"] == "ok"
    assert json.loads(mcp_source_gate(summary))["status"] == "ok"


def test_rejects_unclassified_process_wrong_api_and_meshed_negative():
    bad = copy.deepcopy(_source_summary())
    bad["process"]["acceptable"] = False
    bad["process"]["unexpected_error_lines"] = ["ERROR: sweep failed"]
    bad["api_entity_contract"]["surface_mesh_entity"] = "quad"
    bad["negative_control"]["element_counts"]["hex"] = 144
    result = cubit_mesh_carrying_straight_sweep_source_replay_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["fresh_classified_process"] is False
    assert result["checks"]["surface_mesh_uses_face_entities"] is False
    assert result["checks"]["negative_has_no_volume_elements"] is False
