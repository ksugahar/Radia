import copy
import json
import math

import pytest

from radia_mcp.radia_ngsolve.hartmann_profile_gate import hartmann_profile_gate
from radia_mcp.radia_ngsolve.server import hartmann_profile_gate as mcp_gate


def good_summary():
    hartmann = [1, 10, 100]
    analytic_center = [
        (math.cosh(value) - 1.0)
        / (math.cosh(value) - math.sinh(value) / value)
        for value in hartmann
    ]
    return {
        "units": {
            "hartmann_number": "1",
            "normalized_profile": "1",
            "magnetic_flux_density": "T",
        },
        "hartmann_numbers": hartmann,
        "profile_sample_counts": [201, 401, 801],
        "profile_max_abs_errors": [0.001, 0.002, 0.020],
        "profile_rms_errors": [0.0005, 0.0010, 0.0020],
        "profile_symmetry_errors": [1.0e-13, 2.0e-13, 4.0e-13],
        "normalized_analytic_averages": [1.0, 1.0, 1.0],
        "normalized_center_velocity_fem": [
            analytic_center[0] * 1.001,
            analytic_center[1] * 0.999,
            analytic_center[2] * 0.999,
        ],
        "normalized_center_velocity_analytic": analytic_center,
        "boundary_layer_fractions": [0.9, 0.45, 0.05],
        "magnetic_flux_density_T": [0.001, 0.01, 0.1],
        "gate_tolerances": {
            "profile_linf": 0.025,
            "profile_rms": 0.003,
            "symmetry": 1.0e-8,
            "average": 1.0e-4,
            "center": 0.005,
            "field_scaling": 1.0e-10,
            "minimum_profile_sample_count": 101,
        },
    }


def test_accepts_hartmann_profile_and_boundary_layer_thinning():
    result = hartmann_profile_gate(good_summary())
    assert result["status"] == "ok"
    assert result["checks"]["boundary_layer_thins_monotonically"] is True
    assert json.loads(mcp_gate(json.dumps(good_summary())))["status"] == "ok"


def test_rejects_solver_complete_profile_drift_and_non_thinning_layer():
    bad = copy.deepcopy(good_summary())
    bad["profile_max_abs_errors"][-1] = 0.08
    bad["boundary_layer_fractions"][-1] = 0.7
    result = hartmann_profile_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["independent_profile_linf_errors_pass"] is False
    assert result["checks"]["boundary_layer_thins_monotonically"] is False


def test_rejects_relaxed_policy_and_nonfinite_values():
    relaxed = good_summary()
    relaxed["gate_tolerances"]["profile_linf"] = 0.1
    with pytest.raises(ValueError, match="policy maximum"):
        hartmann_profile_gate(relaxed)

    nonfinite = good_summary()
    nonfinite["profile_rms_errors"][0] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        hartmann_profile_gate(nonfinite)
