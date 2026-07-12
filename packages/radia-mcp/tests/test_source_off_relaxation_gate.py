from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.server import source_off_linear_relaxation_gate as mcp_gate
from radia_mcp.radia_ngsolve.source_off_relaxation_gate import (
    source_off_linear_relaxation_gate,
)


def _summary() -> dict:
    q = 0.44
    current = [1.0, q, q * q, q * q * q]
    return {
        "contract": {
            "resistance_ohm": 8.0,
            "material_response": "linear",
            "response_model": "single_mode_rl_relaxation",
            "source_schedule": "initial_voltage_then_zero",
            "current_semantics": "direct_plus_induced_total_coil_current",
        },
        "rows": [
            {
                "time_s": index * 1.0e-4,
                "source_voltage_v": 8.0 if index == 0 else 0.0,
                "total_coil_current_a": value,
                "field_max_t": 0.43 * value,
            }
            for index, value in enumerate(current)
        ],
    }


def test_source_off_relaxation_accepts_shared_passive_decay_and_dispatches() -> None:
    result = source_off_linear_relaxation_gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["current_and_field_share_one_decay_factor"] is True
    assert result["checks"]["field_scales_linearly_with_total_current"] is True
    assert json.loads(mcp_gate(_summary()))["status"] == "ok"


def test_source_off_relaxation_rejects_lingering_source_and_field_corruption() -> None:
    bad = copy.deepcopy(_summary())
    bad["rows"][2]["source_voltage_v"] = 0.5
    bad["rows"][3]["field_max_t"] *= 1.4
    result = source_off_linear_relaxation_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_active_only_at_initial_sample"] is False
    assert result["checks"]["current_and_field_share_one_decay_factor"] is False
    assert result["checks"]["field_scales_linearly_with_total_current"] is False


def test_source_off_relaxation_rejects_induced_only_current_semantics() -> None:
    bad = _summary()
    bad["contract"]["current_semantics"] = "induced_current_only"
    bad["rows"][0]["total_coil_current_a"] = 0.0
    result = source_off_linear_relaxation_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["total_current_semantics_recorded"] is False
    assert result["checks"]["initial_current_matches_voltage_over_resistance"] is False
