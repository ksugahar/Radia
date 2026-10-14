"""Solver-neutral gate for a linear source-off magnetic relaxation."""

from __future__ import annotations

import math
from typing import Any


def _finite(value: Any, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _relative_span(values: list[float]) -> float:
    scale = max((abs(value) for value in values), default=0.0)
    if not values or scale <= 0.0:
        return math.inf
    return (max(values) - min(values)) / scale


def source_off_linear_relaxation_gate(
    summary: dict[str, Any],
    *,
    max_time_grid_relative_spread: float = 1.0e-9,
    max_initial_ohm_relative_error: float = 1.0e-6,
    max_decay_ratio_relative_span: float = 1.0e-3,
    max_field_current_scale_relative_span: float = 1.0e-3,
    source_off_atol_v: float = 1.0e-12,
    minimum_sample_count: int = 3,
) -> dict[str, Any]:
    """Validate a single-mode linear RL relaxation after the source is removed.

    The first sample is voltage driven from a zero-memory state. Later samples
    have zero applied voltage and must decay with one passive geometric factor.
    A linear magnetic field observable must carry the same factor as coil
    current. The inferred time constant is an effective exponential diagnostic;
    it is not converted to inductance without a declared integration scheme.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("contract")
    rows = summary.get("rows")
    if not isinstance(contract, dict) or not isinstance(rows, list):
        raise ValueError("contract must be a mapping and rows must be a list")

    resistance = _finite(contract.get("resistance_ohm"), "resistance_ohm")
    parsed: list[dict[str, float]] = []
    parse_errors: list[str] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            parse_errors.append(f"row {index} is not a mapping")
            continue
        try:
            parsed.append(
                {
                    "time_s": _finite(row["time_s"], f"row {index}.time_s"),
                    "source_voltage_v": _finite(
                        row["source_voltage_v"], f"row {index}.source_voltage_v"
                    ),
                    "total_coil_current_a": _finite(
                        row["total_coil_current_a"],
                        f"row {index}.total_coil_current_a",
                    ),
                    "field_max_t": _finite(
                        row["field_max_t"], f"row {index}.field_max_t"
                    ),
                }
            )
        except (KeyError, TypeError, ValueError) as exc:
            parse_errors.append(f"row {index}: {exc}")

    times = [row["time_s"] for row in parsed]
    source = [row["source_voltage_v"] for row in parsed]
    currents = [row["total_coil_current_a"] for row in parsed]
    current_magnitudes = [abs(value) for value in currents]
    fields = [abs(row["field_max_t"]) for row in parsed]
    time_steps = [right - left for left, right in zip(times, times[1:])]
    mean_time_step = sum(time_steps) / len(time_steps) if time_steps else math.nan
    time_grid_spread = (
        max(abs(step - mean_time_step) for step in time_steps) / mean_time_step
        if time_steps and mean_time_step > 0.0
        else math.inf
    )
    current_ratios = [
        right / left
        for left, right in zip(current_magnitudes, current_magnitudes[1:])
        if left > 0.0
    ]
    field_ratios = [
        right / left for left, right in zip(fields, fields[1:]) if left > 0.0
    ]
    all_ratios = current_ratios + field_ratios
    decay_ratio_span = _relative_span(all_ratios)
    mean_decay_ratio = (
        sum(all_ratios) / len(all_ratios) if all_ratios else math.nan
    )
    field_per_current = [
        field / current
        for field, current in zip(fields, current_magnitudes, strict=True)
        if current > 0.0
    ]
    field_current_span = _relative_span(field_per_current)
    initial_ohm_error = (
        abs(currents[0] - source[0] / resistance)
        / max(abs(currents[0]), abs(source[0] / resistance), 1.0e-300)
        if parsed and resistance > 0.0
        else math.inf
    )
    effective_tau = (
        -mean_time_step / math.log(mean_decay_ratio)
        if mean_time_step > 0.0 and 0.0 < mean_decay_ratio < 1.0
        else math.nan
    )
    nonzero_signs = [math.copysign(1.0, value) for value in currents if value != 0.0]

    checks = {
        "rows_parsed_and_finite": not parse_errors and len(parsed) == len(rows),
        "sample_count_sufficient": len(parsed) >= int(minimum_sample_count),
        "linear_single_mode_contract_recorded": contract.get("material_response")
        == "linear"
        and contract.get("response_model") == "single_mode_rl_relaxation",
        "source_schedule_recorded": contract.get("source_schedule")
        == "initial_voltage_then_zero",
        "total_current_semantics_recorded": contract.get("current_semantics")
        == "direct_plus_induced_total_coil_current",
        "positive_resistance": resistance > 0.0,
        "strictly_increasing_uniform_time_grid": bool(time_steps)
        and all(step > 0.0 for step in time_steps)
        and time_grid_spread <= float(max_time_grid_relative_spread),
        "source_active_only_at_initial_sample": bool(source)
        and abs(source[0]) > float(source_off_atol_v)
        and all(abs(value) <= float(source_off_atol_v) for value in source[1:]),
        "initial_current_matches_voltage_over_resistance": initial_ohm_error
        <= float(max_initial_ohm_relative_error),
        "current_and_field_resolved": bool(current_magnitudes)
        and min(current_magnitudes) > 0.0
        and min(fields) > 0.0,
        "current_polarity_does_not_flip": bool(nonzero_signs)
        and all(sign == nonzero_signs[0] for sign in nonzero_signs),
        "passive_monotone_decay": len(current_ratios) == len(parsed) - 1
        and len(field_ratios) == len(parsed) - 1
        and all(0.0 < ratio < 1.0 for ratio in all_ratios),
        "current_and_field_share_one_decay_factor": bool(all_ratios)
        and decay_ratio_span <= float(max_decay_ratio_relative_span),
        "field_scales_linearly_with_total_current": bool(field_per_current)
        and field_current_span <= float(max_field_current_scale_relative_span),
        "effective_time_constant_positive": math.isfinite(effective_tau)
        and effective_tau > 0.0,
    }
    issues = parse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "schema": "radia-source-off-linear-relaxation/v1",
        "policy": "linear_source_off_total_current_field_decay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(parsed),
            "time_step_s": mean_time_step,
            "time_grid_relative_spread": time_grid_spread,
            "current_decay_ratios": current_ratios,
            "field_decay_ratios": field_ratios,
            "mean_decay_ratio": mean_decay_ratio,
            "decay_ratio_relative_span": decay_ratio_span,
            "field_per_current_relative_span": field_current_span,
            "initial_ohm_relative_error": initial_ohm_error,
            "effective_exponential_time_constant_s": effective_tau,
        },
        "tolerances": {
            "max_time_grid_relative_spread": float(max_time_grid_relative_spread),
            "max_initial_ohm_relative_error": float(max_initial_ohm_relative_error),
            "max_decay_ratio_relative_span": float(max_decay_ratio_relative_span),
            "max_field_current_scale_relative_span": float(
                max_field_current_scale_relative_span
            ),
            "source_off_atol_v": float(source_off_atol_v),
            "minimum_sample_count": int(minimum_sample_count),
        },
        "notes": [
            "For voltage drive, total coil current includes direct V/R and the induced-current contribution.",
            "A linear field observable should decay with the same factor as total current after source removal.",
            "Do not infer inductance from the decay factor until the time-integration scheme is declared.",
            "Do not compare a nonlinear static companion numerically with this linear transient gate.",
        ],
    }
