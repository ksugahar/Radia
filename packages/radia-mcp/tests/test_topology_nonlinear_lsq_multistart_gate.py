import copy
import json

from radia_mcp.topology_optimization.nonlinear_lsq_multistart_gate import (
    evaluate_nonlinear_lsq_multistart,
)
from radia_mcp.topology_optimization.server import (
    topology_opt_nonlinear_lsq_multistart_gate,
)


def _run(start):
    return {
        "start": start,
        "solution": [0.2, -0.3],
        "residual": [1.0e-12, -2.0e-12],
        "residual_norm": 5.0**0.5 * 1.0e-12,
        "resnorm": 5.0e-24,
        "exitflag": 1,
        "iterations": 6,
        "func_count": 21,
        "projected_gradient_inf_norm": 3.0e-12,
        "solver_jacobian_relative_error": 2.0e-9,
        "finite_difference_jacobian_relative_error": 5.0e-11,
    }


def _summary():
    return {
        "parameter_unit": "normalized",
        "residual_unit": "normalized",
        "contract": {
            "objective": "0.5*||r||^2",
            "gradient": "J^T*r",
            "resnorm_identity": "resnorm=||r||^2",
            "solver_jacobian_semantics": "residual_jacobian",
        },
        "legacy_starts": [[0.0, 0.0]] * 8,
        "runs": [_run(start) for start in ([0.0, 1.0], [1.0, 0.0], [-1.0, 0.5], [0.5, -1.0])],
        "tolerances": {
            "residual_norm": 1.0e-8,
            "projected_gradient_inf_norm": 1.0e-8,
            "solver_jacobian_relative_error": 1.0e-5,
            "finite_difference_jacobian_relative_error": 1.0e-7,
            "resnorm_identity_relative_error": 1.0e-10,
        },
    }


def test_nonlinear_lsq_multistart_accepts_corrected_diverse_runs():
    result = evaluate_nonlinear_lsq_multistart(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["legacy_zero_multiplicative_start_collapse_detected"] is True
    assert result["checks"]["corrected_starts_are_distinct"] is True
    assert json.loads(topology_opt_nonlinear_lsq_multistart_gate(json.dumps(_summary())))["status"] == "ok"


def test_nonlinear_lsq_multistart_rejects_duplicate_starts_and_bad_jacobian():
    bad = copy.deepcopy(_summary())
    bad["runs"][1]["start"] = bad["runs"][0]["start"]
    bad["runs"][2]["solver_jacobian_relative_error"] = 0.2
    result = evaluate_nonlinear_lsq_multistart(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["corrected_starts_are_distinct"] is False
    assert result["checks"]["solver_jacobians_match_analytic"] is False


def test_nonlinear_lsq_multistart_rejects_wrong_objective_gradient_contract():
    bad = _summary()
    bad["contract"]["gradient"] = "2*J^T*r"
    result = evaluate_nonlinear_lsq_multistart(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["least_squares_contract_recorded"] is False
