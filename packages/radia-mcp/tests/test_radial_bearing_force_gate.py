from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.radial_bearing_force_gate import (
    radial_bearing_force_symmetry_gate as gate,
)
from radia_mcp.radia_ngsolve.server import radial_bearing_force_symmetry_gate


def _summary() -> dict:
    return {
        "excitation_order": ["positive_x", "positive_y", "negative_x", "negative_y"],
        "force_method": "weighted_stress_volume_integral",
        "fresh_replay_relative_error": 0.0,
        "cases": [
            {"role": "balanced", "excitation_a": [6, 6, 6, 6], "force_n": [-0.0049, -0.0155, 0], "force_unit": "N", "coordinate_frame": "cartesian"},
            {"role": "positive_y", "excitation_a": [6, 12, 6, 0], "force_n": [0.1830, 399.8711, 0], "force_unit": "N", "coordinate_frame": "cartesian"},
            {"role": "negative_y", "excitation_a": [6, 0, 6, 12], "force_n": [-0.1735, -399.9312, 0], "force_unit": "N", "coordinate_frame": "cartesian"},
        ],
    }


def test_gate_accepts_null_and_mirrored_magnetic_body_force() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["mirrored_axial_relative_error"] < 2.0e-4


def test_gate_rejects_lorentz_overclaim_and_lost_force_covariance() -> None:
    summary = copy.deepcopy(_summary())
    summary["force_method"] = "lorentz_body_integral"
    summary["cases"][2]["force_n"][1] = 399.9
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "magnetic_body_force_method_recorded" in result["issues"]
    assert "mirrored_axial_forces_reverse_sign" in result["issues"]


def test_mcp_wrapper_accepts_and_rejects_invalid_input() -> None:
    assert json.loads(radial_bearing_force_symmetry_gate(json.dumps(_summary())))["status"] == "ok"
    assert json.loads(radial_bearing_force_symmetry_gate("{}"))["status"] == "invalid_input"
