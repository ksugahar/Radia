from __future__ import annotations

import copy
import json
import math

from radia_mcp.radia_ngsolve.server import single_loop_source_normalized_field_gate
from radia_mcp.radia_ngsolve.single_loop_normalized_field_gate import (
    single_loop_source_normalized_field_gate as gate,
)


def _summary() -> dict:
    routes = []
    for route_index, formulation in enumerate(("surface", "volume")):
        field_rows = []
        power_rows = []
        for frequency in (50.0e3, 125.0e3, 200.0e3):
            current = complex(0.2, -1.0e-3 * route_index)
            transfer = 45.0 * (1.0 + 0.02 * route_index)
            hz = -transfer * abs(current) + 1j * 1.0e-3
            hx = 1.0e-3 * hz
            hy = 5.0e-4 * hz
            magnitude = math.sqrt(abs(hx) ** 2 + abs(hy) ** 2 + abs(hz) ** 2)
            field_rows.append({
                "frequency_hz": frequency,
                "current_a": [current.real, current.imag],
                "h_components_a_per_m": [
                    [hx.real, hx.imag],
                    [hy.real, hy.imag],
                    [hz.real, hz.imag],
                ],
                "h_magnitude_a_per_m": magnitude,
            })
            s11 = complex(0.8, 0.1 * route_index)
            stimulated = 0.5
            power_rows.append({
                "frequency_hz": frequency,
                "s11": [s11.real, s11.imag],
                "stimulated_power_w": stimulated,
                "accepted_power_w": stimulated * (1.0 - abs(s11) ** 2),
            })
        routes.append({
            "name": formulation,
            "formulation": formulation,
            "field_rows": field_rows,
            "power_rows": power_rows,
        })
    return {
        "model_contract": {
            "physics": "harmonic_maxwell",
            "single_turn_loop": True,
            "one_port_per_route": True,
            "same_probe_location": True,
            "raw_port_phase_comparable": False,
        },
        "routes": routes,
    }


def test_gate_accepts_source_normalized_field_transfer() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_source_normalized_transfer_relative_gap"] < 0.03


def test_gate_rejects_raw_port_phase_claim_and_field_drift() -> None:
    summary = copy.deepcopy(_summary())
    summary["model_contract"]["raw_port_phase_comparable"] = True
    summary["routes"][1]["field_rows"][1]["h_magnitude_a_per_m"] *= 1.2
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "single_loop_cross_formulation_contract" in result["issues"]
    assert "source_normalized_field_transfer_agrees" in result["issues"]


def test_mcp_wrapper_returns_structured_invalid_input() -> None:
    result = json.loads(single_loop_source_normalized_field_gate('{"routes": []}'))
    assert result["status"] == "invalid_input"
