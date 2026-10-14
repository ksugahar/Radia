import copy
import json

from radia_mcp.differential_forms.gauge_invariance_gate import evaluate_gauge_invariance
from radia_mcp.differential_forms.server import differential_forms_gauge_invariance_gate


def summary():
    b_real = [[2.993e-3, 2.994e-3, 2.995e-3], [-8e-8, -8e-8, -8e-8], [-4e-8, -4e-8, -4e-8]]
    b_imag = [[-9.60e-5, -9.60e-5, -9.60e-5], [2.7e-8, 2.4e-8, 2.1e-8], [-1.04e-6, -1.03e-6, -1.02e-6]]
    return {
        "analysis": "frequency_domain",
        "units": {"magnetic_flux_density": "T", "length": "m", "power": "W", "frequency": "Hz"},
        "frequency_hz": 60.0,
        "radius_m": 1.25e-4,
        "skin_depth_m": 3.0698e-4,
        "applied_flux_density_T": 1.0e-3,
        "air_conductivity_S_per_m": 0.0,
        "b_without_gauge": {"real": b_real, "imag": b_imag},
        "b_with_gauge": {
            "real": [[value * (1.0 - 5.0e-7) for value in row] for row in b_real],
            "imag": [[value * (1.0 - 5.0e-7) for value in row] for row in b_imag],
        },
        "sphere_loss_without_gauge_W": 9.11608e-14,
        "sphere_loss_with_gauge_W": 9.11607e-14,
        "air_loss_without_gauge_W": 3.77e-16,
        "air_loss_with_gauge_W": 2.56e-44,
        "potential_comparison_policy": "not_gated_gauge_dependent",
        "tolerances": {
            "magnetic_field_relative_error": 1.0e-5,
            "conductor_loss_relative_error": 1.0e-5,
            "air_loss_reduction_ratio": 1.0e-20,
            "weak_skin_phase_lag_rad": 0.1,
        },
    }


def test_accepts_physical_field_and_loss_invariance():
    result = evaluate_gauge_invariance(summary())
    assert result["status"] == "ok"
    assert result["checks"]["vector_potential_is_not_used_as_invariant"] is True
    assert json.loads(differential_forms_gauge_invariance_gate(json.dumps(summary())))["status"] == "ok"


def test_rejects_changed_physical_field_and_loss():
    bad = copy.deepcopy(summary())
    bad["b_with_gauge"]["real"][0] = [value * 0.8 for value in bad["b_with_gauge"]["real"][0]]
    bad["sphere_loss_with_gauge_W"] *= 0.8
    result = evaluate_gauge_invariance(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["magnetic_field_is_gauge_invariant"] is False
    assert result["checks"]["conductor_loss_is_gauge_invariant"] is False


def test_rejects_treating_vector_potential_as_invariant():
    bad = summary()
    bad["potential_comparison_policy"] = "compare_A_directly"
    result = evaluate_gauge_invariance(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["vector_potential_is_not_used_as_invariant"] is False
