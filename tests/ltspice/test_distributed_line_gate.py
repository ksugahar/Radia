from __future__ import annotations

import copy

from radia.ltspice.distributed_line_gate import distributed_line_delay_loss_gate
from radia.ltspice.mcp_server import distributed_line_delay_loss_gate as mcp_gate


def _summary() -> dict:
    base = {
        "length": 1.0,
        "inductance_h_per_length": 0.5e-6,
        "capacitance_f_per_length": 70e-12,
    }
    return {
        "cases": [
            {
                "label": "low_resistance",
                **base,
                "series_resistance_ohm_per_length": 0.3,
                "measured_one_way_delay_s": 6.0404701218608094e-9,
                "first_pulse_peak_output_v": 0.4816930890083313,
            },
            {
                "label": "high_resistance",
                **base,
                "series_resistance_ohm_per_length": 300.0,
                "measured_one_way_delay_s": 6.24170488610107e-9,
                "first_pulse_peak_output_v": 0.10230796039104462,
            },
        ],
        "replay": {
            "low_resistance": {"delay_relative_gap": 0.0, "peak_relative_gap": 0.0},
            "high_resistance": {"delay_relative_gap": 0.0, "peak_relative_gap": 0.0},
        },
    }


def test_gate_accepts_lc_delay_and_resistance_attenuation() -> None:
    result = distributed_line_delay_loss_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["series_resistance_ratio_high_to_low"] == 1000.0


def test_gate_rejects_milli_suffix_misinterpretation() -> None:
    summary = copy.deepcopy(_summary())
    summary["cases"][0]["series_resistance_ohm_per_length"] = 300.0
    result = distributed_line_delay_loss_gate(summary)
    assert result["status"] == "needs_attention"
    assert "spice_m_suffix_resistance_ratio_is_one_thousand" in result["issues"]


def test_mcp_wrapper_reports_invalid_input() -> None:
    assert mcp_gate({})["status"] == "invalid_input"
