import copy
import json

from radia_mcp.cubit.power_tools_replay_gate import (
    cubit_partial_volume_hex_diagnosis_gate,
    cubit_power_tools_cleanup_source_replay_gate,
)
from radia_mcp.cubit.server import (
    cubit_partial_volume_hex_diagnosis_gate as mcp_public_gate,
    cubit_power_tools_cleanup_source_replay_gate as mcp_source_gate,
)


def _public_summary():
    return {
        "source_volume_count": 1,
        "prepared_volume_ids": [1, 2, 3],
        "expected_unmeshed_volume_ids": [3],
        "volume_element_counts": {
            "1": {"hex": 505, "tet": 0, "pyramid": 0, "wedge": 0},
            "2": {"hex": 509, "tet": 0, "pyramid": 0, "wedge": 0},
            "3": {"hex": 0, "tet": 0, "pyramid": 0, "wedge": 0},
        },
        "element_counts": {"hex": 1014, "tet": 0, "pyramid": 0, "wedge": 0},
        "quality": {"scaled_jacobian": {"minimum": 0.179}},
        "exported_block_volume_ids": [1, 2, 3],
        "gmsh_export": {
            "mesh_format": "4.1",
            "binary": False,
            "hex_count": 1014,
            "other_volume_count": 0,
        },
        "solver_ready": False,
        "disposition": "reject_partial_or_low_quality",
    }


def _source_summary():
    public = cubit_partial_volume_hex_diagnosis_gate(_public_summary())
    commands = [
        'import acis "knuckle.sat"',
        "healer autoheal body all",
        "split surface 22",
        "split surface 10",
        "webcut volume 1 with plane normal to curve 35 close_to vertex 51",
        "webcut volume 1 with plane normal to curve 35 close_to vertex 49",
        "remove surface 17 extend",
        "remove surface 15 extend",
        "tweak surface 16 offset -0.9",
        "imprint volume all",
        "merge volume all",
        "composite create surface 9 27",
        "composite create surface 6 24",
        "composite create surface 11 25 35",
        "volume all scheme auto",
        "mesh volume all",
    ]
    return {
        "source_kind": "installed-official-help-power-tools-cad",
        "source_name": "knuckle.sat",
        "source_sha256": "a" * 64,
        "source_doc_sha256": {str(index): "b" * 64 for index in range(1, 12)},
        "binary_name": "coreform_cubit.com",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "command_trace": commands,
        "raw_control": {
            "volume_count": 1,
            "volume_scheme": "",
            "element_counts": {"hex": 0, "tet": 0, "pyramid": 0, "wedge": 0},
        },
        "cleanup_replay": {
            "volume_count": 3,
            "volume_schemes": {"1": "sweep", "2": "sweep", "3": ""},
            "element_counts": {"hex": 1014, "tet": 0, "pyramid": 0, "wedge": 0},
        },
        "explicit_sweep_attempt": {
            "command_returned": True,
            "volume_id": 3,
            "volume_hex_count": 0,
            "console_diagnostics": [
                "Internal loop(s) do not have a corner or end node.",
                "Volumes 3 meshing unsuccessful using scheme: sweep",
            ],
        },
        "process": {
            "acceptable": True,
            "result_artifact_fresh": True,
            "unexpected_error_lines": [],
        },
        "public_gate": public,
        "deterministic_replay": {"repeat_count": 2, "stable_fields_match": True},
        "timing_breakdown_s": {
            "headless_replays": 1.0,
            "process_classification": 0.1,
            "independent_export_parse": 0.1,
            "artifact_finalization": 0.1,
        },
    }


def test_accepts_truthful_partial_mesh_diagnosis_without_promotion():
    result = cubit_partial_volume_hex_diagnosis_gate(_public_summary())
    assert result["status"] == "ok"
    assert result["solver_ready"] is False
    assert result["metrics"]["unmeshed_volume_ids"] == [3]
    assert json.loads(mcp_public_gate(_public_summary()))["status"] == "ok"


def test_rejects_solver_ready_overclaim_and_hidden_unmeshed_volume():
    bad = copy.deepcopy(_public_summary())
    bad["solver_ready"] = True
    bad["expected_unmeshed_volume_ids"] = []
    result = cubit_partial_volume_hex_diagnosis_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["partial_mesh_is_explicit"] is False
    assert result["checks"]["solver_handoff_is_suppressed"] is False


def test_rejects_count_quality_and_export_inconsistency():
    bad = copy.deepcopy(_public_summary())
    bad["volume_element_counts"]["2"]["hex"] = 508
    bad["quality"]["scaled_jacobian"]["minimum"] = 0.3
    bad["gmsh_export"]["hex_count"] = 999
    result = cubit_partial_volume_hex_diagnosis_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["owned_counts_match_global_counts"] is False
    assert result["checks"]["quality_rejection_is_measured"] is False
    assert result["checks"]["gmsh_ascii_v41_matches_produced_mesh"] is False


def test_accepts_official_cleanup_trace_and_observed_failure():
    summary = _source_summary()
    result = cubit_power_tools_cleanup_source_replay_gate(summary)
    assert result["status"] == "ok"
    assert result["solver_ready"] is False
    assert json.loads(mcp_source_gate(summary))["status"] == "ok"


def test_rejects_missing_cleanup_step_and_inferred_failure():
    bad = copy.deepcopy(_source_summary())
    bad["command_trace"] = [
        command for command in bad["command_trace"] if "tweak surface" not in command
    ]
    bad["explicit_sweep_attempt"]["command_returned"] = False
    bad["explicit_sweep_attempt"]["console_diagnostics"] = []
    bad["process"]["unexpected_error_lines"] = ["unclassified failure"]
    result = cubit_power_tools_cleanup_source_replay_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["official_cleanup_order_preserved"] is False
    assert result["checks"]["explicit_sweep_failure_is_observed_not_inferred"] is False
    assert result["checks"]["process_errors_are_classified"] is False
