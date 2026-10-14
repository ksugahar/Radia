from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.server import two_body_force_magnitude_replay_gate
from radia_mcp.radia_ngsolve.two_body_force_magnitude_gate import (
    two_body_force_magnitude_replay_gate as gate,
)


def _summary() -> dict:
    row = {
        "body_a_force_magnitude_n": 4.900356631387878e-4,
        "body_b_force_magnitude_n": 4.898997777088313e-4,
        "current_a": 1.0,
        "flux_wb": 3.63798496161428e-6,
        "element_count": 987900,
        "vertex_count": 173363,
        "solver_runtime_s": 491.8,
        "fresh_result": True,
        "force_unit": "N",
        "flux_unit": "Wb",
        "current_unit": "A",
    }
    return {
        "force_quantity": "unsigned_magnitude",
        "action_reaction_sign_inferred": False,
        "commanded_current_a": 1.0,
        "runs": [dict(row, replay=1), dict(row, replay=2, solver_runtime_s=495.3)],
    }


def test_gate_accepts_fresh_unsigned_force_balance_replays() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_body_balance_relative_error"] < 3.0e-4


def test_gate_rejects_sign_overclaim_stale_result_and_imbalance() -> None:
    summary = copy.deepcopy(_summary())
    summary["action_reaction_sign_inferred"] = True
    summary["runs"][1]["fresh_result"] = False
    summary["runs"][1]["body_b_force_magnitude_n"] *= 0.5
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "action_reaction_sign_not_inferred" in result["issues"]
    assert "two_body_force_magnitudes_balance" in result["issues"]
    assert "both_results_are_fresh_solver_outputs" in result["issues"]


def test_mcp_wrapper_accepts_and_reports_invalid_input() -> None:
    assert json.loads(two_body_force_magnitude_replay_gate(json.dumps(_summary())))["status"] == "ok"
    assert json.loads(two_body_force_magnitude_replay_gate("{}"))["status"] == "invalid_input"
