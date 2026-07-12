"""No-load/load transient-cycle validation for three-phase rotating machines."""

from __future__ import annotations

import math
from typing import Any


def _finite(value: Any, name: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _vector(record: dict[str, Any], name: str, count: int = 3) -> list[float]:
    values = record.get(name)
    if not isinstance(values, list) or len(values) != count:
        raise ValueError(f"{name} must contain {count} values")
    return [_finite(value, name) for value in values]


def motor_transient_no_load_load_cycle_gate(summary: dict[str, Any]) -> dict[str, Any]:
    """Validate paired no-load and loaded one-cycle motor evidence.

    The gate combines electrical balance, mechanical kinematics, periodic
    endpoint state, and cycle-mean power closure. It also records the common
    post-step convention where the initial speed row is a zero placeholder and
    angle must be integrated with right-endpoint speed values.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    no_load = summary.get("no_load")
    loaded = summary.get("loaded")
    if not isinstance(no_load, dict) or not isinstance(loaded, dict):
        raise ValueError("no_load and loaded must be mappings")

    errors: list[str] = []
    try:
        frequency = _finite(summary["electrical_frequency_hz"], "electrical_frequency_hz")
        angle_span = _finite(summary["mechanical_angle_span_deg"], "mechanical_angle_span_deg")
        speed_rpm = _finite(summary["final_speed_rpm"], "final_speed_rpm")
        no_current = no_load["phase_current_a"]
        no_torque = no_load["torque_nm"]
        no_voltage = no_load["phase_voltage_v"]
        no_flux = no_load["phase_flux_wb"]
        load_current = loaded["phase_current_a"]
        load_torque = loaded["torque_nm"]
        load_voltage = loaded["phase_voltage_v"]
        load_flux = loaded["phase_flux_wb"]
        load_power = loaded["power_w"]
        no_kinematics = no_load["kinematics"]
        load_kinematics = loaded["kinematics"]
        no_duration = _finite(no_load["cycle_duration_s"], "no_load.cycle_duration_s")
        load_duration = _finite(loaded["cycle_duration_s"], "loaded.cycle_duration_s")
        load_current_rms = _vector(load_current, "rms")
        no_voltage_rms = _vector(no_voltage, "rms")
        no_flux_rms = _vector(no_flux, "rms")
        load_voltage_rms = _vector(load_voltage, "rms")
        load_flux_rms = _vector(load_flux, "rms")
    except (KeyError, TypeError, ValueError) as exc:
        return {
            "policy": "motor_transient_no_load_load_cycle_gate_v1",
            "status": "needs_attention",
            "checks": {},
            "errors": [str(exc)],
        }

    mechanical_frequency = speed_rpm / 60.0
    inferred_pole_pairs = frequency / max(mechanical_frequency, 1.0e-30)
    expected_angle = 360.0 * mechanical_frequency * load_duration
    expected_current_rms = _finite(load_current["peak_abs"], "loaded peak current") / math.sqrt(2.0)
    rms_peak_errors = [
        abs(value - expected_current_rms) / max(expected_current_rms, 1.0e-30)
        for value in load_current_rms
    ]
    first_step = _finite(load_kinematics["first_step_angle_deg"], "first_step_angle_deg")
    trapezoid_offset = _finite(
        load_kinematics["trapezoid_speed_angle_max_error_deg"],
        "trapezoid_speed_angle_max_error_deg",
    )
    metrics = {
        "inferred_pole_pairs": inferred_pole_pairs,
        "angle_from_speed_duration_deg": expected_angle,
        "max_loaded_current_rms_peak_relative_error": max(rms_peak_errors),
        "loaded_current_sum_max_abs_a": _finite(load_current["sum_max_abs"], "loaded current sum"),
        "loaded_power_balance_relative_error": _finite(load_power["balance_relative_error"], "power balance"),
        "loaded_torque_mean_nm": _finite(load_torque["mean"], "loaded torque mean"),
        "no_load_torque_mean_abs_nm": abs(_finite(no_torque["mean"], "no-load torque mean")),
        "no_load_current_peak_abs_a": _finite(no_current["peak_abs"], "no-load current peak"),
        "trapezoid_offset_to_half_step_error_deg": abs(trapezoid_offset - 0.5 * first_step),
    }
    checks = {
        "schema_recorded": summary.get("schema") == "motor.transient-no-load-load-cycle.v1",
        "paired_cycle_sample_counts_sufficient": int(no_load.get("sample_count", 0)) >= 100
        and int(loaded.get("sample_count", 0)) >= 100,
        "one_electrical_cycle_recorded": abs(no_duration - 1.0 / frequency) <= 1.0e-12
        and abs(load_duration - 1.0 / frequency) <= 1.0e-12,
        "integer_pole_pair_kinematics": abs(inferred_pole_pairs - round(inferred_pole_pairs)) <= 1.0e-9
        and round(inferred_pole_pairs) >= 1
        and abs(angle_span - expected_angle) <= 1.0e-8,
        "right_endpoint_speed_angle_convention": _finite(
            no_kinematics["right_endpoint_speed_angle_max_error_deg"], "no-load right endpoint"
        ) <= 1.0e-8
        and _finite(load_kinematics["right_endpoint_speed_angle_max_error_deg"], "load right endpoint") <= 1.0e-8
        and metrics["trapezoid_offset_to_half_step_error_deg"] <= 1.0e-9,
        "no_load_open_circuit_behavior": metrics["no_load_current_peak_abs_a"] <= 5.0e-6
        and metrics["no_load_torque_mean_abs_nm"] <= 1.0e-3,
        "no_load_voltage_and_flux_balanced": _finite(no_voltage["rms_spread_relative"], "no-load voltage spread") <= 0.01
        and _finite(no_flux["rms_spread_relative"], "no-load flux spread") <= 0.01
        and all(value > 0.0 for value in no_voltage_rms + no_flux_rms),
        "loaded_four_amp_three_phase_balance": abs(_finite(load_current["peak_abs"], "loaded peak") - 4.0) <= 1.0e-4
        and metrics["max_loaded_current_rms_peak_relative_error"] <= 1.0e-6
        and _finite(load_current["rms_spread_relative"], "loaded current spread") <= 1.0e-5
        and metrics["loaded_current_sum_max_abs_a"] <= 1.0e-5,
        "loaded_voltage_and_flux_balanced": _finite(load_voltage["rms_spread_relative"], "loaded voltage spread") <= 0.01
        and _finite(load_flux["rms_spread_relative"], "loaded flux spread") <= 0.01
        and all(value > 0.0 for value in load_voltage_rms + load_flux_rms),
        "loaded_periodic_endpoint_state": _finite(load_torque["endpoint_relative_error"], "torque endpoint") <= 1.0e-3
        and _finite(load_current["endpoint_relative_error_max"], "current endpoint") <= 1.0e-5
        and _finite(load_flux["endpoint_relative_error_max"], "flux endpoint") <= 1.0e-4,
        "loaded_positive_torque_and_power": metrics["loaded_torque_mean_nm"] > 0.0
        and _finite(load_power["electrical_cycle_mean"], "electrical power") > 0.0
        and _finite(load_power["mechanical_cycle_mean"], "mechanical power") > 0.0
        and _finite(load_power["copper_loss_cycle_mean"], "copper loss") > 0.0,
        "loaded_cycle_power_balance_within_3pct": metrics["loaded_power_balance_relative_error"] <= 0.03,
    }
    failed = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "motor_transient_no_load_load_cycle_gate_v1",
        "status": "ok" if not errors and not failed else "needs_attention",
        "checks": checks,
        "metrics": metrics,
        "errors": errors + failed,
        "lesson": (
            "Validate a rotating-machine cycle as a no-load/load pair. Open-circuit current and mean torque, "
            "loaded three-phase RMS/KCL, periodic endpoint state, speed-angle convention, and cycle power "
            "balance must agree before the waveform is accepted."
        ),
    }
