"""Solver-neutral gate for a loaded two-stage Cockcroft-Walton multiplier."""
from __future__ import annotations

import math


def cockcroft_walton_stage_gate(
    vin_peak_v: float,
    stage1_avg_v: float,
    stage2_avg_v: float,
    stage2_previous_avg_v: float,
    stage1_ripple_vpp: float,
    stage2_ripple_vpp: float,
    load_ohm: float,
    load_avg_a: float,
    source_power_delivered_w: float,
    load_power_w: float,
    max_stage_law_relative_error: float = 0.05,
    max_settling_relative_drift: float = 0.01,
    max_ripple_fraction: float = 0.05,
) -> dict:
    """Gate stage scaling, settling, load balance, ripple, and real-power bounds."""

    values = (
        vin_peak_v,
        stage1_avg_v,
        stage2_avg_v,
        stage2_previous_avg_v,
        stage1_ripple_vpp,
        stage2_ripple_vpp,
        load_ohm,
        load_avg_a,
        source_power_delivered_w,
        load_power_w,
        max_stage_law_relative_error,
        max_settling_relative_drift,
        max_ripple_fraction,
    )
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("all multiplier metrics and tolerances must be finite")
    if any(float(value) <= 0.0 for value in values[:10]):
        raise ValueError("all multiplier measurements must be positive")
    if any(float(value) < 0.0 for value in values[10:]):
        raise ValueError("all tolerances must be nonnegative")

    ideal_stage2_v = 4.0 * vin_peak_v
    stage_law_error = abs(stage2_avg_v - 2.0 * stage1_avg_v) / stage2_avg_v
    settling_drift = abs(stage2_avg_v - stage2_previous_avg_v) / stage2_avg_v
    stage1_ripple_fraction = stage1_ripple_vpp / stage1_avg_v
    stage2_ripple_fraction = stage2_ripple_vpp / stage2_avg_v
    expected_load_current = stage2_avg_v / load_ohm
    load_current_error = abs(load_avg_a - expected_load_current) / expected_load_current
    efficiency = load_power_w / source_power_delivered_w

    checks = {
        "stage_voltage_strictly_increases": vin_peak_v < stage1_avg_v < stage2_avg_v,
        "two_stage_voltage_law": stage_law_error <= max_stage_law_relative_error,
        "loaded_output_below_ideal_ceiling": stage2_avg_v <= 1.01 * ideal_stage2_v,
        "loaded_output_above_70_percent_ideal": stage2_avg_v >= 0.70 * ideal_stage2_v,
        "late_windows_settled": settling_drift <= max_settling_relative_drift,
        "stage1_ripple_bounded": stage1_ripple_fraction <= max_ripple_fraction,
        "stage2_ripple_bounded": stage2_ripple_fraction <= max_ripple_fraction,
        "load_ohm_law_closes": load_current_error <= 0.01,
        "source_power_covers_load": source_power_delivered_w >= load_power_w,
        "efficiency_physical": 0.0 < efficiency <= 1.0,
    }
    return {
        "policy": "cockcroft_walton_stage_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "ideal_stage2_v": ideal_stage2_v,
            "stage_law_relative_error": stage_law_error,
            "settling_relative_drift": settling_drift,
            "stage1_ripple_fraction": stage1_ripple_fraction,
            "stage2_ripple_fraction": stage2_ripple_fraction,
            "expected_load_current_a": expected_load_current,
            "load_current_relative_error": load_current_error,
            "real_power_efficiency": efficiency,
        },
        "notes": [
            "Use late, adjacent integer-cycle windows so startup charging is not mistaken for steady-state gain.",
            "Loaded diode ladders should remain below the ideal 4*Vin peak ceiling; exact diode loss is model-dependent.",
        ],
    }
