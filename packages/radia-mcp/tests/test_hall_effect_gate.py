from __future__ import annotations

import copy
import json

import numpy as np

from radia_mcp.radia_ngsolve.hall_effect_gate import (
    hall_effect_transverse_voltage_gate as gate,
)
from radia_mcp.radia_ngsolve.server import hall_effect_transverse_voltage_gate


def _summary() -> dict:
    angle = np.linspace(-90.0, 90.0, 37)
    field = 0.045 * np.exp(-(angle / 28.0) ** 2) - 0.002
    voltage = 4.8 * field
    return {
        "angle_deg": angle.tolist(),
        "hall_voltage_baseline_v": voltage.tolist(),
        "hall_voltage_replay_v": voltage.tolist(),
        "hall_voltage_zero_coefficient_v": (1.0e-6 * voltage).tolist(),
        "hall_voltage_reversed_coefficient_v": (-voltage).tolist(),
        "hall_voltage_scaled_drive_v": (0.5 * voltage).tolist(),
        "drive_scale_ratio": 0.5,
        "magnetic_flux_density_baseline_t": field.tolist(),
        "magnetic_flux_density_replay_t": field.tolist(),
        "magnetic_flux_density_zero_coefficient_t": field.tolist(),
        "magnetic_flux_density_reversed_coefficient_t": field.tolist(),
        "magnetic_flux_density_scaled_drive_t": field.tolist(),
    }


def test_gate_accepts_hall_constitutive_controls() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["first_harmonic_relative_residual"] > 0.1


def test_gate_rejects_false_replay_coefficient_and_drive_claims() -> None:
    summary = copy.deepcopy(_summary())
    center = len(summary["angle_deg"]) // 2
    scale = max(abs(value) for value in summary["hall_voltage_baseline_v"])
    summary["hall_voltage_replay_v"][center] += 0.02 * scale
    summary["hall_voltage_zero_coefficient_v"][center] = 0.1 * scale
    summary["hall_voltage_reversed_coefficient_v"][center] = summary[
        "hall_voltage_baseline_v"
    ][center]
    summary["hall_voltage_scaled_drive_v"][center] = summary[
        "hall_voltage_baseline_v"
    ][center]
    summary["magnetic_flux_density_scaled_drive_t"][center] += 0.01
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "fresh_voltage_replay_is_deterministic" in result["issues"]
    assert "zero_hall_coefficient_suppresses_transverse_voltage" in result["issues"]
    assert "hall_coefficient_reversal_reverses_transverse_voltage" in result["issues"]
    assert "drive_scaling_is_linear" in result["issues"]
    assert "prescribed_magnetic_field_is_case_invariant" in result["issues"]


def test_mcp_wrapper_accepts_hall_sweep() -> None:
    result = json.loads(hall_effect_transverse_voltage_gate(json.dumps(_summary())))
    assert result["status"] == "ok"


def test_mcp_wrapper_reports_invalid_input() -> None:
    result = json.loads(hall_effect_transverse_voltage_gate("{}"))
    assert result["status"] == "invalid_input"
