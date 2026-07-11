import copy
import json

from radia_mcp.topology_optimization.server import topology_opt_simplex_stationarity_audit_gate
from radia_mcp.topology_optimization.simplex_stationarity_gate import evaluate_simplex_stationarity_audit


def _summary():
    return {
        "parameter_unit": "normalized",
        "objective_unit": "normalized",
        "reference": {"x": [0.0, -0.5], "objective": -0.25},
        "tolerances": {"gradient_norm": 1.0e-5, "objective_gap": 1.0e-10, "parameter_distance": 1.0e-5},
        "methods": [
            {
                "method_id": "pathological_simplex",
                "role": "candidate",
                "reported_converged": True,
                "x": [0.0, 0.0],
                "objective": 0.0,
                "gradient": [0.0, 1.0],
                "function_evaluations": 95,
            },
            {
                "method_id": "independent_control",
                "role": "control",
                "reported_converged": True,
                "x": [0.0, -0.5],
                "objective": -0.25,
                "gradient": [0.0, 0.0],
                "function_evaluations": 200,
            },
        ],
    }


def test_simplex_stationarity_audit_detects_false_convergence():
    result = evaluate_simplex_stationarity_audit(_summary())
    assert result["status"] == "ok"
    assert result["false_convergence_method_ids"] == ["pathological_simplex"]
    assert result["false_convergence_candidate_ids"] == ["pathological_simplex"]
    assert result["accepted_method_ids"] == ["independent_control"]
    assert result["accepted_control_method_ids"] == ["independent_control"]
    assert json.loads(topology_opt_simplex_stationarity_audit_gate(json.dumps(_summary())))["status"] == "ok"


def test_simplex_stationarity_audit_rejects_missing_accepted_control():
    bad = copy.deepcopy(_summary())
    bad["methods"][1]["gradient"] = [0.0, 0.1]
    result = evaluate_simplex_stationarity_audit(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["accepted_independent_control_present"] is False


def test_simplex_stationarity_audit_rejects_duplicate_method_ids():
    bad = copy.deepcopy(_summary())
    bad["methods"][1]["method_id"] = bad["methods"][0]["method_id"]
    result = evaluate_simplex_stationarity_audit(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["method_ids_unique"] is False
