from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.energy_budgeted_trace_kkt_gate import (
    energy_budgeted_trace_kkt_gate,
)
from radia_mcp.radia_ngsolve.server import energy_budgeted_trace_kkt_gate as mcp_gate


def _solver(algorithm: str, stationarity: float = 1.0e-8) -> dict:
    return {
        "algorithm": algorithm,
        "exitflag": 1,
        "objective": 2.0,
        "energy": 0.4,
        "constraint_value": -1.0e-11,
        "dual_lambda": 3.0,
        "stationarity_inf": stationarity,
        "complementarity_abs": 3.0e-11,
        "solution_relative_error": 2.0e-8,
        "objective_relative_error": 3.0e-9,
    }


def _payload() -> dict:
    return {
        "mesh": {
            "volume_element": "tet4",
            "boundary_element": "tri3",
            "volume_basis": "H1_P1",
            "trace_basis": "scalar_P1",
            "points": 100,
            "trace_dofs": 60,
            "tets": 300,
            "triangles": 120,
        },
        "problem": {"energy_budget": 0.4, "budget_fraction": 0.35},
        "analytic": {
            "energy": 0.4,
            "constraint_value": 0.0,
            "stationarity_inf": 1.0e-12,
            "complementarity_abs": 0.0,
        },
        "metrics": {
            "objective_gradient_relative_error": 1.0e-9,
            "constraint_gradient_relative_error": 2.0e-9,
            "solver_pair_solution_relative_difference": 4.0e-8,
        },
        "solvers": [_solver("sqp"), _solver("interior-point")],
    }


def test_energy_budgeted_trace_gate_accepts_full_kkt_closure() -> None:
    result = energy_budgeted_trace_kkt_gate(_payload())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_energy_budgeted_trace_gate_rejects_nonstationary_solver() -> None:
    payload = copy.deepcopy(_payload())
    payload["solvers"][0]["stationarity_inf"] = 1.0e-3
    result = energy_budgeted_trace_kkt_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["both_solvers_stationary"] is False


def test_energy_budgeted_trace_gate_rejects_negative_dual() -> None:
    payload = copy.deepcopy(_payload())
    payload["solvers"][1]["dual_lambda"] = -1.0
    result = energy_budgeted_trace_kkt_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["both_solvers_positive_dual"] is False


def test_energy_budgeted_trace_mcp_tool_dispatches() -> None:
    result = json.loads(mcp_gate(_payload()))
    assert result["status"] == "ok"
