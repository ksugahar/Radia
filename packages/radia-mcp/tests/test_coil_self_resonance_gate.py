from __future__ import annotations

import copy
import json

import numpy as np

from radia_mcp.radia_ngsolve.coil_self_resonance_gate import (
    coil_self_resonance_sweep_gate as gate,
)
from radia_mcp.radia_ngsolve.server import coil_self_resonance_sweep_gate


def _summary() -> dict:
    frequency = np.arange(1.0e6, 10.0e6 + 0.25e6, 0.25e6)
    resonance = 6.383e6
    resistance = 80.0 + 6.0e4 / (1.0 + ((frequency - 6.5e6) / 0.3e6) ** 2)
    reactance = 1.0e3 * (1.0 - frequency / resonance)
    return {
        "frequency_hz": frequency.tolist(),
        "resistance_ohm": resistance.tolist(),
        "reactance_ohm": reactance.tolist(),
        "dataset_frequency_relative_error": 0.0,
        "dataset_resistance_relative_error": 0.0,
        "dataset_reactance_relative_error": 0.0,
        "replay_frequency_relative_error": 0.0,
        "replay_resistance_relative_error": 1.0e-8,
        "replay_reactance_relative_error": 2.0e-8,
    }


def test_gate_accepts_passive_self_resonance_sweep() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert 6.25e6 < result["metrics"]["interpolated_self_resonance_hz"] < 6.5e6


def test_gate_rejects_missing_capacitive_regime_and_negative_resistance() -> None:
    summary = copy.deepcopy(_summary())
    summary["reactance_ohm"] = [abs(value) for value in summary["reactance_ohm"]]
    summary["resistance_ohm"][0] = -1.0
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "positive_series_resistance_is_passive" in result["issues"]
    assert "single_inductive_to_capacitive_transition" in result["issues"]


def test_mcp_wrapper_reports_invalid_input() -> None:
    assert json.loads(coil_self_resonance_sweep_gate("{}"))["status"] == "invalid_input"


def test_mcp_wrapper_serializes_numpy_backed_checks() -> None:
    result = json.loads(coil_self_resonance_sweep_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
