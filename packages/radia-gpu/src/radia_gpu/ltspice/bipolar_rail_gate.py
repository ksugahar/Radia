"""Solver-neutral bipolar rail power-quality gate."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


_POLICY_LIMITS = {
    "regulation": 0.05,
    "balance": 0.02,
    "ripple": 0.02,
    "efficiency_closure_percent": 1.0e-6,
}


def _finite(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _nonnegative(value: object, name: str) -> float:
    parsed = _finite(value, name)
    if parsed < 0.0:
        raise ValueError(f"{name} must be nonnegative")
    return parsed


def bipolar_rail_power_quality_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate signed bipolar rails, ripple, target regulation, and power closure."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    units = summary.get("units")
    window = summary.get("measure_window_s")
    if not isinstance(units, Mapping):
        raise ValueError("units must be a mapping")
    if not isinstance(window, Sequence) or isinstance(window, (str, bytes)) or len(window) != 2:
        raise ValueError("measure_window_s must contain start and stop")

    start_s = _nonnegative(window[0], "measure_window_s[0]")
    stop_s = _nonnegative(window[1], "measure_window_s[1]")
    target_v = _finite(summary.get("target_rail_voltage_v"), "target_rail_voltage_v")
    positive_v = _finite(summary.get("positive_output_voltage_v"), "positive_output_voltage_v")
    negative_v = _finite(summary.get("negative_output_voltage_v"), "negative_output_voltage_v")
    positive_ripple_v = _nonnegative(
        summary.get("positive_output_ripple_pp_v"), "positive_output_ripple_pp_v"
    )
    negative_ripple_v = _nonnegative(
        summary.get("negative_output_ripple_pp_v"), "negative_output_ripple_pp_v"
    )
    input_power_w = _finite(summary.get("delivered_input_power_w"), "delivered_input_power_w")
    positive_power_w = _finite(summary.get("positive_output_power_w"), "positive_output_power_w")
    negative_power_w = _finite(summary.get("negative_output_power_w"), "negative_output_power_w")
    efficiency_percent = _finite(
        summary.get("reported_efficiency_percent"), "reported_efficiency_percent"
    )
    limits = {
        "regulation": _nonnegative(
            summary.get("max_regulation_relative_error", _POLICY_LIMITS["regulation"]),
            "max_regulation_relative_error",
        ),
        "balance": _nonnegative(
            summary.get("max_rail_imbalance_relative", _POLICY_LIMITS["balance"]),
            "max_rail_imbalance_relative",
        ),
        "ripple": _nonnegative(
            summary.get("max_ripple_fraction", _POLICY_LIMITS["ripple"]),
            "max_ripple_fraction",
        ),
        "efficiency_closure_percent": _nonnegative(
            summary.get(
                "max_efficiency_closure_percent",
                _POLICY_LIMITS["efficiency_closure_percent"],
            ),
            "max_efficiency_closure_percent",
        ),
    }
    if any(limits[name] > _POLICY_LIMITS[name] for name in limits):
        raise ValueError("limits cannot exceed the policy maxima")

    positive_magnitude_v = positive_v
    negative_magnitude_v = -negative_v
    mean_rail_v = (positive_magnitude_v + negative_magnitude_v) / 2.0
    rail_imbalance = (
        abs(positive_magnitude_v - negative_magnitude_v) / mean_rail_v
        if mean_rail_v > 0.0
        else math.inf
    )
    positive_regulation = (
        abs(positive_magnitude_v - target_v) / target_v if target_v > 0.0 else math.inf
    )
    negative_regulation = (
        abs(negative_magnitude_v - target_v) / target_v if target_v > 0.0 else math.inf
    )
    positive_ripple_fraction = (
        positive_ripple_v / positive_magnitude_v if positive_magnitude_v > 0.0 else math.inf
    )
    negative_ripple_fraction = (
        negative_ripple_v / negative_magnitude_v if negative_magnitude_v > 0.0 else math.inf
    )
    output_power_w = positive_power_w + negative_power_w
    power_loss_w = input_power_w - output_power_w
    recomputed_efficiency = (
        output_power_w / input_power_w * 100.0 if input_power_w > 0.0 else math.nan
    )
    efficiency_closure = abs(efficiency_percent - recomputed_efficiency)

    checks = {
        "units_are_explicit": units.get("voltage") == "V"
        and units.get("power") == "W"
        and units.get("efficiency") == "percent"
        and units.get("time") == "s",
        "measure_window_is_ordered": start_s < stop_s,
        "bipolar_polarity_is_signed": positive_v > 0.0 and negative_v < 0.0,
        "target_voltage_is_positive": target_v > 0.0,
        "both_rails_meet_target_regulation": max(positive_regulation, negative_regulation)
        <= limits["regulation"],
        "rail_magnitudes_are_balanced": rail_imbalance <= limits["balance"],
        "both_rail_ripples_are_small": max(
            positive_ripple_fraction, negative_ripple_fraction
        )
        <= limits["ripple"],
        "both_outputs_deliver_positive_power": positive_power_w > 0.0
        and negative_power_w > 0.0,
        "power_balance_is_passive": input_power_w > output_power_w > 0.0
        and power_loss_w >= 0.0,
        "reported_efficiency_is_physical": 0.0 < efficiency_percent <= 100.0,
        "reported_efficiency_recomputes": efficiency_closure
        <= limits["efficiency_closure_percent"],
    }
    return {
        "policy": "bipolar_rail_power_quality_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "limits": limits,
        "metrics": {
            "positive_rail_magnitude_v": positive_magnitude_v,
            "negative_rail_magnitude_v": negative_magnitude_v,
            "positive_regulation_relative_error": positive_regulation,
            "negative_regulation_relative_error": negative_regulation,
            "rail_imbalance_relative": rail_imbalance,
            "positive_ripple_fraction": positive_ripple_fraction,
            "negative_ripple_fraction": negative_ripple_fraction,
            "output_power_w": output_power_w,
            "power_loss_w": power_loss_w,
            "recomputed_efficiency_percent": recomputed_efficiency,
            "efficiency_closure_absolute_percent": efficiency_closure,
        },
        "lesson": (
            "A bipolar converter is operating credibly only when the rails have opposite signs, both magnitudes "
            "meet the same target, ripple and imbalance are small in one late window, output power remains below "
            "delivered input power, and reported efficiency recomputes from that window."
        ),
    }
