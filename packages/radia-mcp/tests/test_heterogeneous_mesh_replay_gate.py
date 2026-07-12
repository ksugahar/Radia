import copy
import json

from radia_mcp.radia_ngsolve.heterogeneous_mesh_replay_gate import (
    heterogeneous_part_mesh_replay_gate,
)
from radia_mcp.radia_ngsolve.server import heterogeneous_part_mesh_replay_gate as mcp_gate


def good_summary() -> dict:
    replays = []
    for name in ("fresh_a", "fresh_b"):
        replays.append(
            {
                "replay_id": name,
                "nodes": 5919,
                "elements": 31857,
                "source_preserved": True,
                "temporary_work_copy": True,
                "has_mesh_before": False,
                "has_mesh_after": True,
                "has_mesh_any_part_after": True,
                "has_result_after": False,
                "changed_artifacts": ["mesh_input", "mesh_log", "mesh_output"],
                "pass_marker": True,
                "owned_processes_after": 0,
            }
        )
    return {
        "reference_evidence": [
            {
                "evidence_id": "archive_a",
                "nodes": 71442,
                "elements": 278634,
                "independent": True,
            },
            {
                "evidence_id": "archive_b",
                "nodes": 71442,
                "elements": 278634,
                "independent": True,
            },
        ],
        "mesh_routes": [
            {"route": "extruded_part_mesh", "observed": True},
            {"route": "external_part_mesh", "observed": True},
            {"route": "automatic_volume_mesh", "observed": True},
        ],
        "warning": {
            "code": "duplicate_part_healing_condition",
            "count": 1,
            "observed_in_report": True,
            "disposition": "source_configuration_requires_review",
        },
        "live_replays": replays,
        "classification": "reproducible_remesh_drift_not_solver_ready",
        "solver_ready": False,
        "maximum_reference_live_relative_error": 0.05,
    }


def test_accepts_complete_reproducible_drift_diagnosis_without_solver_ready_claim():
    result = heterogeneous_part_mesh_replay_gate(good_summary())
    assert result["status"] == "ok"
    assert result["solver_ready"] is False
    assert result["metrics"]["element_relative_drift"] > 0.88
    assert json.loads(mcp_gate(json.dumps(good_summary())))["status"] == "ok"


def test_rejects_false_solver_ready_claim_and_ignored_warning():
    bad = copy.deepcopy(good_summary())
    bad["solver_ready"] = True
    bad["warning"]["disposition"] = "ignore"
    result = heterogeneous_part_mesh_replay_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_warning_is_classified_not_ignored"] is False
    assert result["checks"]["material_drift_is_not_promoted_to_solver_ready"] is False


def test_rejects_nonrepeatable_live_counts_and_missing_external_route():
    bad = copy.deepcopy(good_summary())
    bad["live_replays"][1]["elements"] = 34000
    bad["mesh_routes"][1]["observed"] = False
    result = heterogeneous_part_mesh_replay_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["fresh_replay_mesh_counts_are_deterministic"] is False
    assert result["checks"]["heterogeneous_mesh_routes_are_explicit"] is False
