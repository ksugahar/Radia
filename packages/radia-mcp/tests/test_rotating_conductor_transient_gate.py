from __future__ import annotations

import json

from radia_mcp.radia_ngsolve.rotating_conductor_transient_gate import (
    rotating_conductor_transient_gate as gate,
)
from radia_mcp.radia_ngsolve.server import rotating_conductor_transient_gate


def _summary() -> dict:
    time = [0.01 * index for index in range(21)]
    speed = 60.0
    omega = 2.0 * 3.141592653589793
    return {
        "moving_axis_boundary": {
            "axis": "z",
            "before": "open",
            "after": "magnetic",
            "source_modified": False,
        },
        "units": {
            "time": "s",
            "angle": "rad",
            "speed": "r/min",
            "torque": "N*m",
            "loss": "W",
            "current_flux": "A",
        },
        "angle_rows": [[value, omega * value] for value in time],
        "speed_rows": [[value, speed] for value in time[1:]],
        "torque_rows": [[value, 0.2 + 0.01 * index] for index, value in enumerate(time)],
        "loss_rows": [
            {"time_s": value, "total_w": 3.0, "parts_w": [3.0, 0.0]}
            for value in time[1:]
        ],
        "current_flux_rows": [[value, (-1.0) ** index] for index, value in enumerate(time[1:])],
        "energy_balance_contract": {
            "mechanical_power_vs_joule_loss": "diagnostic_only",
            "external_drive_power_available": False,
        },
    }


def test_rotating_conductor_transient_gate_accepts_right_endpoint_tables():
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_kinematic_relative_error"] < 1.0e-13
    assert result["metrics"]["maximum_loss_partition_relative_error"] == 0.0


def test_rotating_conductor_transient_gate_rejects_open_axis_after_migration():
    summary = _summary()
    summary["moving_axis_boundary"]["after"] = "open"
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "moving_axis_boundary_migrated_closed" in result["issues"]


def test_rotating_conductor_transient_gate_rejects_shifted_speed_time_axis():
    summary = _summary()
    summary["speed_rows"][0][0] += 0.001
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "right_endpoint_tables_align" in result["issues"]


def test_rotating_conductor_transient_gate_rejects_false_energy_claim():
    summary = _summary()
    summary["energy_balance_contract"]["mechanical_power_vs_joule_loss"] = "must_equal"
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "mechanical_power_vs_joule_loss_is_diagnostic" in result["issues"]


def test_rotating_conductor_transient_gate_is_exposed_over_mcp_wrapper():
    result = json.loads(rotating_conductor_transient_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
