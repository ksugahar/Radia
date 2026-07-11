from __future__ import annotations

from radia_mcp.radia_ngsolve.field_profile_gate import (
    dual_formulation_symmetric_field_profile_gate,
)


def _summary():
    return {
        "component_id": "B_axis",
        "field_unit": "T",
        "axis_unit": "m",
        "sample_count": 201,
        "axis_min": -1.0,
        "axis_max": 1.0,
        "center_axis": 0.0,
        "center_value_a": -5.61756e-9,
        "center_value_b": -5.61627e-9,
        "profile_relative_l2_difference": 9.42e-4,
        "center_relative_difference": 2.31e-4,
        "symmetry_relative_a": 2.51e-4,
        "symmetry_relative_b": 2.63e-4,
    }


def test_dual_formulation_profile_gate_accepts_agreement_and_symmetry():
    result = dual_formulation_symmetric_field_profile_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_dual_formulation_profile_gate_rejects_single_point_only_or_asymmetry():
    bad = _summary()
    bad["sample_count"] = 1
    bad["profile_relative_l2_difference"] = 0.2
    bad["symmetry_relative_b"] = 0.3
    result = dual_formulation_symmetric_field_profile_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "sample_count_sufficient",
        "profile_formulations_agree",
        "formulation_b_is_symmetric",
    }


def test_dual_formulation_profile_gate_rejects_trivial_zero_field():
    bad = _summary()
    bad["center_value_a"] = 0.0
    bad["center_value_b"] = 0.0
    assert dual_formulation_symmetric_field_profile_gate(bad)["checks"]["center_field_nonzero"] is False
