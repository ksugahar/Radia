import copy
import json

from radia_mcp.radia_ngsolve.autodiff_harmonic_balance_gate import (
    autodiff_harmonic_balance_convergence_gate,
)
from radia_mcp.radia_ngsolve.server import (
    autodiff_harmonic_balance_convergence_gate as mcp_gate,
)


def _summary():
    return {
        "problem": {
            "periodic_point_count": 1024,
            "odd_harmonic_maximum": 61,
            "variable_count": 63,
            "periodic_endpoint_duplicated": False,
        },
        "false_convergence_control": {
            "mean_residual_abs": 1.7e-17,
            "rms_residual": 2.0**-0.5,
            "mean_criterion_threshold": 1.0e-9,
        },
        "jacobian": {
            "ad_analytic_relative_error": 3.3e-17,
            "ad_complex_step_relative_error": 1.8e-16,
            "rows": 1024,
            "columns": 63,
        },
        "exact_recovery": {
            "ad_coefficient_relative_error": 6.8e-13,
            "reference_coefficient_relative_error": 1.1e-16,
            "ad_reference_coefficient_relative_error": 6.8e-13,
            "ad_replay_coefficient_relative_error": 0.0,
            "ad_final_rms_residual": 1.7e-13,
            "reference_final_rms_residual": 2.8e-17,
            "reference_exitflag": 1,
        },
        "target_fit": {
            "initial_rms_residual": 0.7071,
            "final_rms_residual": 6.6e-5,
            "final_gradient_inf_norm": 5.2e-10,
        },
        "stopping": {
            "residual_rms_used": True,
            "gradient_inf_used": True,
            "residual_mean_only_used": False,
        },
        "timing_breakdown_s": {
            "setup": 0.01,
            "jacobian": 0.08,
            "solvers": 0.18,
            "packaging": 0.01,
        },
    }


def test_accepts_triple_checked_ad_harmonic_balance_and_dispatches():
    result = autodiff_harmonic_balance_convergence_gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["mean_only_rule_is_a_demonstrated_false_convergence"]
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_mean_only_stopping_and_duplicate_periodic_endpoint():
    bad = copy.deepcopy(_summary())
    bad["problem"]["periodic_endpoint_duplicated"] = True
    bad["stopping"]["residual_rms_used"] = False
    bad["stopping"]["residual_mean_only_used"] = True
    result = autodiff_harmonic_balance_convergence_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["periodic_grid_has_no_duplicate_endpoint"] is False
    assert (
        result["checks"]["rms_and_gradient_stopping_replace_mean_only_stopping"]
        is False
    )


def test_rejects_jacobian_drift_solver_disagreement_and_replay_drift():
    bad = copy.deepcopy(_summary())
    bad["jacobian"]["ad_complex_step_relative_error"] = 1.0e-4
    bad["exact_recovery"]["ad_reference_coefficient_relative_error"] = 1.0e-3
    bad["exact_recovery"]["ad_replay_coefficient_relative_error"] = 1.0e-5
    result = autodiff_harmonic_balance_convergence_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["ad_matches_complex_step_jacobian"] is False
    assert result["checks"]["independent_nonlinear_solvers_agree"] is False
    assert result["checks"]["fresh_ad_replay_is_deterministic"] is False
