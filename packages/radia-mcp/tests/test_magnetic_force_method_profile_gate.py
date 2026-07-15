from __future__ import annotations

import copy
import json

import pytest

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from radia_mcp.radia_ngsolve.server import (
    magnetic_force_method_profile_gate as mcp_magnetic_force_method_profile_gate,
)


def _summary() -> dict:
    return {
        "quantity_dimension": "3d_total",
        "force_unit": "N",
        "position_unit": "m",
        "comparison_axis": "z",
        "positions": [0.0, 0.001, 0.002, 0.003, 0.004, 0.005],
        "moving_body_element_force": [10.0, 11.0, 12.5, 14.0, 16.0, 18.0],
        "closed_surface_maxwell_stress_force": [10.2, 11.2, 12.7, 14.2, 16.2, 18.2],
        "independent_closed_surface_force": [10.1, 11.1, 12.6, 14.1, 16.1, 18.1],
        "all_body_element_force": [3.0, 3.1, 3.2, 3.4, 3.6, 3.8],
        "replay": {
            "parsed_max_abs": 0.0,
            "binary_nonlog_outputs_exact": True,
        },
    }


def test_accepts_pinned_target_body_and_closed_surface_profiles() -> None:
    result = magnetic_force_method_profile_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_method_relative_difference"] < 0.05
    assert result["metrics"]["minimum_selection_scope_relative_difference"] > 0.25


def test_mcp_tool_dispatches_json() -> None:
    result = json.loads(mcp_magnetic_force_method_profile_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
    assert result["policy"] == "magnetic_force_method_profile_gate_v1"


def test_rejects_unpinned_all_body_selection() -> None:
    summary = copy.deepcopy(_summary())
    summary["all_body_element_force"] = list(summary["moving_body_element_force"])
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["selection_scope_is_materially_distinct"] is False
    assert result["checks"]["all_body_control_is_not_target_body_force"] is False


def test_rejects_method_disagreement_and_nonexact_replay() -> None:
    summary = copy.deepcopy(_summary())
    summary["closed_surface_maxwell_stress_force"][2] = 20.0
    summary["replay"]["parsed_max_abs"] = 1.0e-6
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["independent_stress_replay_within_tolerance"] is False
    assert result["checks"]["parsed_replay_is_exact_enough"] is False


def test_rejects_profile_length_mismatch() -> None:
    summary = _summary()
    summary["all_body_element_force"].pop()
    with pytest.raises(ValueError, match="same length"):
        magnetic_force_method_profile_gate(summary)


def test_rejects_independent_surface_outlier_and_replay_drift() -> None:
    summary = copy.deepcopy(_summary())
    summary["independent_closed_surface_force"][2] *= 1.5
    summary["replay"]["parsed_max_abs"] = 1.0e-4
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["independent_stress_replay_within_tolerance"] is False
    assert result["checks"]["parsed_replay_is_exact_enough"] is False


@pytest.mark.parametrize(
    "case_id",
    ["quantity_dimension", "force_unit", "selection_control", "stress_method", "binary_replay"],
)
def test_counterfactual_curriculum90_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "quantity_dimension":
        summary["quantity_dimension"] = "2d_per_length"
    elif case_id == "force_unit":
        summary["force_unit"] = "mN"
    elif case_id == "selection_control":
        summary["all_body_element_force"] = list(summary["moving_body_element_force"])
    elif case_id == "stress_method":
        summary["closed_surface_maxwell_stress_force"][2] *= 2.0
    else:
        summary["replay"]["binary_nonlog_outputs_exact"] = False
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_generalization_v3s_rejects_unsupported_position_unit() -> None:
    summary = copy.deepcopy(_summary())
    summary["position_unit"] = "inch"
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v4_invalid_comparison_axis", "v4_position_order", "v4_element_force_sign", "v4_element_force_nonfinite", "v4_missing_target_signal"],
)
def test_counterfactual_curriculum90_v4_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v4_invalid_comparison_axis":
        summary["comparison_axis"] = "unsupported"
    elif case_id == "v4_position_order":
        summary["positions"][2] = summary["positions"][1]
    elif case_id == "v4_element_force_sign":
        summary["moving_body_element_force"][1] *= -1.0
    elif case_id == "v4_element_force_nonfinite":
        summary["moving_body_element_force"][4] = float("nan")
    else:
        summary["moving_body_element_force"] = [0.0] * len(summary["moving_body_element_force"])
    result = json.loads(mcp_magnetic_force_method_profile_gate(json.dumps(summary)))
    assert result["status"] in {"needs_attention", "invalid_input"}


def test_generalization_v5_rejects_short_target_force_profile() -> None:
    summary = copy.deepcopy(_summary())
    summary["moving_body_element_force"].pop()
    with pytest.raises(ValueError, match="same length"):
        magnetic_force_method_profile_gate(summary)


@pytest.mark.parametrize(
    "case_id",
    ["v6_public_closed_surface_force_drift", "v6_public_force_dimension_unit_mismatch"],
)
def test_generalization_v6_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v6_public_closed_surface_force_drift":
        summary["closed_surface_maxwell_stress_force"][3] *= 1.30
    else:
        summary["force_unit"] = "N/m"
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"
