from __future__ import annotations

import copy
import json

import numpy as np

from radia_mcp.radia_ngsolve.passive_axial_bearing_gate import (
    passive_axial_bearing_stiffness_gate as gate,
)
from radia_mcp.radia_ngsolve.server import passive_axial_bearing_stiffness_gate


def _summary() -> dict:
    position = np.linspace(-0.04, 0.04, 41)
    stiffness = 1.5e5
    saturation = 1.0 + (position / 0.012) ** 4
    force0 = stiffness * position / saturation
    force1 = -force0 + 1.0e-4 * np.max(np.abs(force0))
    return {
        "position_m": position.tolist(),
        "force_object_0_n": force0.tolist(),
        "force_object_1_n": force1.tolist(),
        "expected_axial_stability": "unstable",
        "fresh_replay_relative_error": 5.0e-14,
    }


def test_gate_accepts_symmetric_axially_unstable_force_curve() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["observed_axial_stability"] == "unstable"
    assert result["metrics"]["center_signed_stiffness_df_dx_n_per_m"] > 0.0


def test_gate_rejects_stability_overclaim_and_broken_force_balance() -> None:
    summary = copy.deepcopy(_summary())
    summary["expected_axial_stability"] = "stable"
    summary["force_object_1_n"][18] = 0.0
    summary["fresh_replay_relative_error"] = 1.0e-4
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "observed_stability_matches_expectation" in result["issues"]
    assert "action_reaction_closes" in result["issues"]
    assert "fresh_force_curve_replay_is_stable" in result["issues"]


def test_mcp_wrapper_reports_invalid_input() -> None:
    result = json.loads(passive_axial_bearing_stiffness_gate("{}"))
    assert result["status"] == "invalid_input"


def test_mcp_wrapper_accepts_force_curve() -> None:
    result = json.loads(passive_axial_bearing_stiffness_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
