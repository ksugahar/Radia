import copy
import json

import pytest

from radia_mcp.cubit.ato_sculpt_gate import (
    cubit_ato_levelset_sculpt_source_replay_gate,
    cubit_levelset_sculpt_hex_validation_gate,
)
from radia_mcp.cubit.server import (
    cubit_ato_levelset_sculpt_source_replay_gate as mcp_source_gate,
    cubit_levelset_sculpt_hex_validation_gate as mcp_public_gate,
)


def _public_summary():
    return {
        "source_tet_count": 52472,
        "iso_triangle_count": 10852,
        "mbg_volume": 29.998377346,
        "mesh_series": [
            {
                "label": "coarse",
                "cell_size": 0.5,
                "hex_count": 804,
                "other_volume_count": 0,
                "source_minimum_quality": 0.18425,
            },
            {
                "label": "fine",
                "cell_size": 0.125,
                "hex_count": 22473,
                "other_volume_count": 0,
                "source_minimum_quality": 0.233243,
            },
        ],
        "mesh_ready": True,
        "solver_ready": False,
        "disposition": "accept_mesh_require_physics_label_handoff",
    }


def _source_summary():
    public = cubit_levelset_sculpt_hex_validation_gate(_public_summary())
    commands = [
        "set dev on",
        'import mesh "small_bracket.exo" nodal_var "LSD" no_geom',
        'create tri iso tet all nodal_var "LSD"',
        'export mesh "small_bracket_iso.e" block 3 4 overwrite',
        "reset",
        'import mesh "small_bracket_iso.e" no_geom',
        "node in tri in block 3 position fixed",
        "smooth tri in block 4 target free mesh iteration 5",
        "create mesh geom tri all",
        "create vol surf all",
        'save cub5 "small_bracket_mbg.cub5" overwrite journal',
        'open "small_bracket_mbg.cub5"',
        "sculpt volume all processors 1 size 0.5",
        'open "small_bracket_mbg.cub5"',
        "sculpt volume all processors 1 size 0.125",
    ]
    return {
        "source_kind": "installed-official-help-ato-levelset-exodus",
        "source_name": "small_bracket.exo",
        "source_sha256": "a" * 64,
        "source_doc_name": "ato_to_mesh.htm",
        "source_doc_sha256": "b" * 64,
        "binary_name": "coreform_cubit.com",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "source_mesh": {"node_count": 10021, "tet_count": 52472},
        "iso_surface": {"block_3_tri_count": 2704, "block_4_tri_count": 8148},
        "command_trace": commands,
        "obsolete_stl_negative_control": {
            "command_returned": True,
            "artifact_exists": False,
            "console_diagnostics": [
                "Could not export the data to a file.",
                "stl is not a valid type of file to be exported.",
            ],
            "migration": "export iso blocks to Exodus, reconstruct MBG, then sculpt volume",
        },
        "process": {
            "acceptable": True,
            "result_artifact_fresh": True,
            "unexpected_error_lines": [],
            "owned_processes_remaining": 0,
        },
        "public_gate": public,
        "deterministic_replay": {"repeat_count": 2, "stable_fields_match": True},
        "timing_breakdown_s": {
            "two_headless_native_replays": 1.0,
            "process_and_replay_classification": 0.1,
            "independent_quality_checks": 0.1,
            "result_finalization": 0.1,
        },
    }


def test_accepts_refined_all_hex_mesh_without_solver_overclaim():
    result = cubit_levelset_sculpt_hex_validation_gate(_public_summary())
    assert result["status"] == "ok"
    assert result["mesh_ready"] is True
    assert result["solver_ready"] is False
    assert json.loads(mcp_public_gate(_public_summary()))["status"] == "ok"


def test_rejects_quality_and_solver_overclaims():
    bad = copy.deepcopy(_public_summary())
    bad["mesh_series"][1]["source_minimum_quality"] = 0.19
    bad["solver_ready"] = True
    result = cubit_levelset_sculpt_hex_validation_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_quality_crosses_configured_floor"] is False
    assert result["checks"]["mesh_ready_without_solver_ready_overclaim"] is False


def test_accepts_official_ato_mbg_sculpt_replay():
    summary = _source_summary()
    result = cubit_ato_levelset_sculpt_source_replay_gate(summary)
    assert result["status"] == "ok"
    assert result["mesh_ready"] is True
    assert result["solver_ready"] is False
    assert json.loads(mcp_source_gate(summary))["status"] == "ok"


def test_rejects_missing_levelset_migration_and_hidden_stl_failure():
    bad = copy.deepcopy(_source_summary())
    bad["command_trace"] = [
        command for command in bad["command_trace"] if "nodal_var" not in command
    ]
    bad["obsolete_stl_negative_control"]["console_diagnostics"] = []
    bad["obsolete_stl_negative_control"]["artifact_exists"] = True
    bad["process"]["unexpected_error_lines"] = ["unclassified failure"]
    result = cubit_ato_levelset_sculpt_source_replay_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["mbg_sculpt_migration_order_preserved"] is False
    assert result["checks"]["obsolete_stl_failure_is_observed"] is False
    assert result["checks"]["headless_process_errors_are_classified"] is False


def test_rejects_stale_source_counts_and_unmatched_replay():
    bad = copy.deepcopy(_source_summary())
    bad["source_mesh"]["tet_count"] -= 1
    bad["iso_surface"]["block_4_tri_count"] -= 1
    bad["deterministic_replay"]["stable_fields_match"] = False
    result = cubit_ato_levelset_sculpt_source_replay_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_levelset_counts_recorded"] is False
    assert result["checks"]["documented_iso_blocks_recorded"] is False
    assert result["checks"]["two_replays_match"] is False


def test_rejects_gui_daemon_and_owned_process_leak():
    bad = copy.deepcopy(_source_summary())
    bad["gui_daemon_enabled"] = True
    bad["process"]["owned_processes_remaining"] = 1
    result = cubit_ato_levelset_sculpt_source_replay_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["headless_execution_recorded"] is False
    assert result["checks"]["headless_process_errors_are_classified"] is False


@pytest.mark.parametrize(
    "case_id",
    ["v6_source_projection_replay_drift", "v6_source_public_gate_disagreement"],
)
def test_generalization_v6_source(case_id):
    bad = copy.deepcopy(_source_summary())
    if case_id == "v6_source_projection_replay_drift":
        bad["deterministic_replay"]["stable_fields_match"] = False
    else:
        bad["public_gate"]["status"] = "needs_attention"
    assert cubit_ato_levelset_sculpt_source_replay_gate(bad)["status"] == "needs_attention"
