"""Solver-neutral validation of periodic hysteresis-loss power tables."""

from __future__ import annotations

import math
from typing import Any


def _finite(row: dict[str, Any], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"{key} must be finite")
    return value


def _trapz(values: list[float], dt: float) -> float:
    return sum(0.5 * (left + right) * dt for left, right in zip(values, values[1:]))


def periodic_hysteresis_loss_energy_gate(
    summary: dict[str, Any],
    *,
    max_waveform_repeat_relative_error: float = 1.0e-3,
    max_cycle_energy_relative_error: float = 1.0e-3,
    max_loss_decomposition_relative_error: float = 1.0e-10,
    max_part_total_relative_error: float = 1.0e-10,
) -> dict[str, Any]:
    """Gate two final steady cycles of instantaneous magnetic-loss power.

    Instantaneous hysteresis power may be negative during local energy return.
    Dissipation is established by a positive closed-cycle integral, not by
    requiring every power sample to be nonnegative.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("contract") or {}
    if not isinstance(contract, dict):
        raise ValueError("contract must be a mapping")
    rows = summary.get("rows") or []
    if not isinstance(rows, list):
        raise ValueError("rows must be a list")

    parsed: list[dict[str, float]] = []
    parse_errors: list[str] = []
    keys = (
        "time_s",
        "joule_loss_w",
        "hysteresis_loss_w",
        "iron_loss_w",
        "hysteresis_part_w",
        "hysteresis_total_w",
    )
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            parse_errors.append(f"row {index} is not a mapping")
            continue
        try:
            parsed.append({key: _finite(row, key) for key in keys})
        except (KeyError, TypeError, ValueError) as exc:
            parse_errors.append(f"row {index}: {exc}")

    times = [row["time_s"] for row in parsed]
    increments = [right - left for left, right in zip(times, times[1:])]
    increasing = bool(times) and all(value > 0.0 for value in increments)
    dt = sum(increments) / len(increments) if increments else math.nan
    uniform = bool(increments) and max(abs(value - dt) for value in increments) <= max(
        1.0e-12, abs(dt) * 1.0e-9
    )
    period = float(contract.get("cycle_period_s", math.nan))
    steps_per_cycle = round(period / dt) if uniform and dt > 0.0 and period > 0.0 else 0
    period_on_grid = (
        steps_per_cycle > 0
        and abs(steps_per_cycle * dt - period) <= max(1.0e-12, period * 1.0e-9)
    )
    enough_cycles = len(parsed) >= 2 * steps_per_cycle + 1 if steps_per_cycle else False

    previous: list[dict[str, float]] = []
    latest: list[dict[str, float]] = []
    if enough_cycles:
        previous = parsed[-2 * steps_per_cycle - 1 : -steps_per_cycle]
        latest = parsed[-steps_per_cycle - 1 :]

    previous_power = [row["hysteresis_loss_w"] for row in previous]
    latest_power = [row["hysteresis_loss_w"] for row in latest]
    waveform_scale = max(
        [abs(value) for value in previous_power + latest_power], default=0.0
    )
    waveform_repeat = (
        max(abs(left - right) for left, right in zip(previous_power, latest_power))
        / waveform_scale
        if previous_power and latest_power and waveform_scale > 0.0
        else math.inf
    )
    previous_energy = _trapz(previous_power, dt) if previous_power else math.nan
    latest_energy = _trapz(latest_power, dt) if latest_power else math.nan
    energy_repeat = (
        abs(previous_energy - latest_energy)
        / max(abs(previous_energy), abs(latest_energy), 1.0e-30)
        if previous_power and latest_power
        else math.inf
    )

    loss_scale = max(
        [
            abs(row["iron_loss_w"])
            + abs(row["joule_loss_w"])
            + abs(row["hysteresis_loss_w"])
            for row in parsed
        ],
        default=0.0,
    )
    decomposition_error = max(
        [
            abs(
                row["iron_loss_w"]
                - row["joule_loss_w"]
                - row["hysteresis_loss_w"]
            )
            / max(loss_scale, 1.0e-30)
            for row in parsed
        ],
        default=math.inf,
    )
    part_total_scale = max(
        [
            max(abs(row["hysteresis_part_w"]), abs(row["hysteresis_total_w"]))
            for row in parsed
        ],
        default=0.0,
    )
    part_total_error = max(
        [
            abs(row["hysteresis_part_w"] - row["hysteresis_total_w"])
            / max(part_total_scale, 1.0e-30)
            for row in parsed
        ],
        default=math.inf,
    )
    hysteresis_values = [row["hysteresis_loss_w"] for row in parsed]

    checks = {
        "rows_parsed_and_finite": not parse_errors and len(parsed) == len(rows),
        "time_strictly_increasing": increasing,
        "time_step_uniform": uniform,
        "units_are_seconds_and_watts": contract.get("time_unit") == "s"
        and contract.get("power_unit") == "W",
        "cycle_period_is_on_time_grid": period_on_grid,
        "two_final_complete_cycles_present": enough_cycles,
        "steady_waveform_repeats": waveform_repeat
        <= float(max_waveform_repeat_relative_error),
        "cycle_energies_are_positive": math.isfinite(previous_energy)
        and math.isfinite(latest_energy)
        and previous_energy > 0.0
        and latest_energy > 0.0,
        "cycle_energy_repeats": energy_repeat
        <= float(max_cycle_energy_relative_error),
        "instantaneous_power_has_return_interval": bool(hysteresis_values)
        and min(hysteresis_values) < 0.0
        and max(hysteresis_values) > 0.0,
        "iron_loss_decomposition_closes": decomposition_error
        <= float(max_loss_decomposition_relative_error),
        "single_part_total_closes": contract.get("single_part_total") is True
        and part_total_error <= float(max_part_total_relative_error),
    }
    issues = parse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "policy": "periodic_hysteresis_loss_energy_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(parsed),
            "time_step_s": dt,
            "steps_per_cycle": steps_per_cycle,
            "previous_cycle_energy_j": previous_energy,
            "latest_cycle_energy_j": latest_energy,
            "waveform_repeat_relative_error": waveform_repeat,
            "cycle_energy_relative_error": energy_repeat,
            "loss_decomposition_relative_error": decomposition_error,
            "part_total_relative_error": part_total_error,
            "instantaneous_power_min_w": min(hysteresis_values, default=math.nan),
            "instantaneous_power_max_w": max(hysteresis_values, default=math.nan),
        },
        "tolerances": {
            "max_waveform_repeat_relative_error": float(
                max_waveform_repeat_relative_error
            ),
            "max_cycle_energy_relative_error": float(max_cycle_energy_relative_error),
            "max_loss_decomposition_relative_error": float(
                max_loss_decomposition_relative_error
            ),
            "max_part_total_relative_error": float(max_part_total_relative_error),
        },
        "notes": [
            "Negative instantaneous hysteresis power can represent local energy return; require positive closed-cycle energy instead of clipping it.",
            "Compare phase-aligned complete cycles on a uniform time grid before declaring periodic steady state.",
            "Check both loss decomposition and part-to-total closure so a numerically plausible total does not hide a namespace or domain mismatch.",
        ],
    }
