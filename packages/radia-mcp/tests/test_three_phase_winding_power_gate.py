from __future__ import annotations

import cmath
import json
import math

import pytest

from radia_mcp.radia_ngsolve.server import three_phase_winding_power_balance_gate
from radia_mcp.radia_ngsolve.three_phase_winding_power_gate import (
    three_phase_winding_power_balance_gate as gate,
)


def _triplet(magnitude: float) -> list[list[float]]:
    return [
        [
            (magnitude * cmath.exp(-1j * 2.0 * math.pi * index / 3.0)).real,
            (magnitude * cmath.exp(-1j * 2.0 * math.pi * index / 3.0)).imag,
        ]
        for index in range(3)
    ]


def _summary() -> dict:
    passive_current = math.sqrt(99.9)
    return {
        "phasor_convention": "rms",
        "expected_phase_step_deg": -120.0,
        "voltage_unit": "V",
        "current_unit": "A",
        "resistance_unit": "ohm",
        "source_winding": {
            "label": "source",
            "voltage_phasors": _triplet(100.0),
            "current_phasors": _triplet(1.0),
            "phase_resistance_ohm": 0.1,
        },
        "passive_windings": [
            {
                "label": "shorted pickup",
                "current_phasors": _triplet(passive_current),
                "phase_resistance_ohm": 1.0,
            }
        ],
    }


def test_three_phase_winding_power_gate_accepts_balanced_power_closure():
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["complex_input_power_va"]["real"] == pytest.approx(300.0)
    assert result["metrics"]["copper_loss_w"] == pytest.approx(300.0)
    assert result["metrics"]["active_power_relative_residual"] < 1.0e-14


def test_three_phase_winding_power_gate_rejects_loss_mismatch():
    summary = _summary()
    summary["passive_windings"][0]["phase_resistance_ohm"] = 0.5
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "active_power_closes_to_copper_loss" in result["issues"]


def test_three_phase_winding_power_gate_rejects_wrong_source_direction():
    summary = _summary()
    summary["source_winding"]["current_phasors"] = [
        [-real, -imag]
        for real, imag in summary["source_winding"]["current_phasors"]
    ]
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["active_input_power_positive"] is False


def test_three_phase_winding_power_gate_is_exposed_over_mcp_wrapper():
    payload = json.loads(three_phase_winding_power_balance_gate(json.dumps(_summary())))
    assert payload["status"] == "ok"


def test_three_phase_winding_power_gate_validates_units_and_shape():
    summary = _summary()
    summary["source_winding"]["voltage_phasors"] = [[1.0, 0.0]]
    with pytest.raises(ValueError, match="exactly three"):
        gate(summary)
