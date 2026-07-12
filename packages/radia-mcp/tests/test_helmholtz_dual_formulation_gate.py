import copy
import json
import math

from radia_mcp.radia_ngsolve.helmholtz_dual_formulation_gate import (
    helmholtz_dual_formulation_axis_gate,
)
from radia_mcp.radia_ngsolve.server import helmholtz_dual_formulation_axis_gate as mcp_gate


def good_summary():
    axis = [(-10 + index) * 0.05 for index in range(21)]
    amplitude = -5.6e-9
    primary = [amplitude / (1.0 + (value / 0.35) ** 4) for value in axis]
    gradient = [
        amplitude
        * (-4.0 * value**3 / 0.35**4)
        / (1.0 + (value / 0.35) ** 4) ** 2
        for value in axis
    ]
    return {
        "units": {"axis": "m", "field": "T", "gradient": "T/m"},
        "axis_m": axis,
        "primary_field_T": primary,
        "secondary_field_T": [value * 0.9998 for value in primary],
        "primary_gradient_T_per_m": gradient,
        "secondary_gradient_T_per_m": [value * 0.999 for value in gradient],
        "gate_tolerances": {
            "maximum_field_symmetry_relative": 0.002,
            "maximum_formulation_field_relative_error": 0.002,
            "maximum_formulation_gradient_relative_error": 0.01,
            "maximum_center_gradient_normalized": 1.0e-4,
            "maximum_center_curvature_normalized": 0.002,
            "maximum_gradient_odd_symmetry_relative": 0.05,
            "maximum_edge_to_center_abs_ratio": 0.6,
        },
    }


def test_accepts_symmetric_flat_dual_formulation_profile():
    result = helmholtz_dual_formulation_axis_gate(good_summary())
    assert result["status"] == "ok"
    assert result["checks"]["center_curvature_is_small"] is True
    assert json.loads(mcp_gate(json.dumps(good_summary())))["status"] == "ok"


def test_rejects_formulation_drift_and_center_curvature():
    bad = copy.deepcopy(good_summary())
    bad["secondary_field_T"][10] *= 0.95
    bad["primary_field_T"][9] *= 0.98
    bad["primary_field_T"][11] *= 0.98
    result = helmholtz_dual_formulation_axis_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["field_formulations_agree"] is False
    assert result["checks"]["center_curvature_is_small"] is False


def test_rejects_nonfinite_profile():
    bad = good_summary()
    bad["primary_gradient_T_per_m"][4] = math.nan
    try:
        helmholtz_dual_formulation_axis_gate(bad)
    except ValueError as exc:
        assert "must be finite" in str(exc)
    else:
        raise AssertionError("nonfinite gradient was accepted")
