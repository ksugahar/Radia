from __future__ import annotations

import copy
import json
import math

from radia_mcp.radia_ngsolve.server import two_winding_frequency_faraday_gate
from radia_mcp.radia_ngsolve.two_winding_frequency_gate import (
    two_winding_frequency_faraday_gate as gate,
)


def _summary() -> dict:
    rows = []
    for frequency in (100.0, 1000.0, 10000.0):
        base = complex(1.0e-4, -2.0e-5) / (1.0 + frequency / 2000.0)
        windings = []
        for turns, resistance, factor in ((40.0, 1.2, 1.0), (10.0, 0.5, 0.98)):
            flux = turns * factor * base
            response = (-1j * 2.0 * math.pi * frequency * flux) / resistance
            windings.append({
                "turns": turns,
                "resistance_ohm": resistance,
                "response": [response.real, response.imag],
                "flux_linkage_Wb_turn": [flux.real, flux.imag],
            })
        rows.append({"frequency_hz": frequency, "windings": windings})
    return {
        "model_contract": {
            "physics": "harmonic_magnetics",
            "two_windings": True,
            "passive_secondary": True,
            "complex_phasors": True,
        },
        "rows": rows,
    }


def test_gate_accepts_two_winding_faraday_identity() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["frequency_count"] == 3
    assert result["metrics"]["maximum_faraday_relative_error"] < 1.0e-12


def test_gate_rejects_phase_and_frequency_order_faults() -> None:
    summary = copy.deepcopy(_summary())
    summary["rows"][1]["windings"][0]["response"][1] *= -1.0
    summary["rows"][2]["frequency_hz"] = 500.0
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "faraday_identity_holds_for_both_windings" in result["issues"]
    assert "positive_strictly_increasing_frequency_axis" in result["issues"]


def test_mcp_wrapper_returns_structured_invalid_input() -> None:
    result = json.loads(two_winding_frequency_faraday_gate('{"rows": []}'))
    assert result["status"] == "invalid_input"
