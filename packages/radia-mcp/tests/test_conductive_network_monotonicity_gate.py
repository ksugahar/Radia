import copy
import json

from radia_mcp.radia_ngsolve.conductive_network_monotonicity_gate import (
    conductive_network_resistance_monotonicity_gate,
)
from radia_mcp.radia_ngsolve.server import (
    conductive_network_resistance_monotonicity_gate as mcp_gate,
)


def _summary():
    cases = []
    for case_id, count, resistance, dof, solve_s in (
        ("one", 1, 10.0, 1000, 1.0),
        ("four", 4, 8.0, 4000, 4.0),
        ("five", 5, 7.8, 4500, 4.5),
    ):
        current = 1.0
        cases.append(
            {
                "case_id": case_id,
                "contacting_conductor_count": count,
                "effective_resistance_ohm": resistance,
                "terminal_power_W": resistance * current * current,
                "joule_loss_W": resistance * current * current,
                "current_balance_relative_error": 1e-9,
                "adaptive_relative_error": 1e-3,
                "final_log10_residual": -7.0,
                "final_dof": dof,
                "solve_s": solve_s,
            }
        )
    return {
        "cases": cases,
        "replay_max_relative_error": 1e-12,
        "timing_breakdown_s": {"inventory": 1.0, "solve": 10.0, "read": 1.0, "verify": 1.0},
    }


def test_accepts_rayleigh_monotone_contact_family_and_dispatches():
    result = conductive_network_resistance_monotonicity_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["minimum_drop_to_adaptive_error_ratio"] > 20.0
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_resistance_increase_and_power_imbalance():
    bad = copy.deepcopy(_summary())
    bad["cases"][2]["effective_resistance_ohm"] = 8.5
    bad["cases"][2]["terminal_power_W"] = 20.0
    result = conductive_network_resistance_monotonicity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["effective_resistance_strictly_decreases"] is False
    assert result["checks"]["all_topologies_are_internally_valid"] is False


def test_rejects_drop_inside_adaptive_error_and_replay_drift():
    bad = _summary()
    bad["cases"][2]["effective_resistance_ohm"] = 7.99
    bad["replay_max_relative_error"] = 1e-4
    result = conductive_network_resistance_monotonicity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_resistance_drops_exceed_adaptive_error"] is False
    assert result["checks"]["independent_replay_is_deterministic"] is False
