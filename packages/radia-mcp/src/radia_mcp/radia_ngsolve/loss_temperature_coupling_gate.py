"""Solver-neutral loss-to-temperature coupling validation."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import Any


def _finite(row: Mapping[str, Any], name: str) -> float:
    try:
        value = float(row[name])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"missing or invalid {name}") from exc
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def loss_temperature_coupling_gate(
    magnetic_rows: Iterable[Mapping[str, Any]],
    thermal_rows: Iterable[Mapping[str, Any]],
    *,
    loss_to_heat_scale: float,
    coupling_rtol: float = 2.0e-5,
    decomposition_rtol: float = 1.0e-12,
    minimum_power_coverage: float = 0.90,
    initial_temperature_c: float = 20.0,
    initial_temperature_atol_c: float = 1.0e-9,
) -> dict[str, object]:
    """Gate an induction-heating loss table and its thermal response.

    Thermal row zero is the initial state; each later heat-source row maps to
    one magnetic loss row through ``loss_to_heat_scale``. The scale is explicit
    because duty cycle, symmetry, axial depth, or cycle averaging may make the
    thermal source a reduced version of the electromagnetic loss.
    """

    magnetic = [dict(row) for row in magnetic_rows]
    thermal = [dict(row) for row in thermal_rows]
    if len(magnetic) < 5 or len(thermal) != len(magnetic) + 1:
        raise ValueError("thermal rows must contain initial state plus one row per magnetic row")
    scale = float(loss_to_heat_scale)
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("loss_to_heat_scale must be finite and positive")

    magnetic_values = []
    decomposition_errors = []
    power_coverages = []
    for row in magnetic:
        gear = _finite(row, "target_loss_w")
        auxiliary = _finite(row, "auxiliary_loss_w")
        other = _finite(row, "other_loss_w")
        total = _finite(row, "total_loss_w")
        active = _finite(row, "active_input_power_w")
        components = gear + auxiliary + other
        decomposition_errors.append(abs(total - components) / max(abs(total), 1.0e-300))
        power_coverages.append(total / active if active > 0.0 else math.inf)
        magnetic_values.append(gear)

    times = [_finite(row, "time_s") for row in thermal]
    heat = [_finite(row, "target_heat_source_w") for row in thermal]
    temperatures = [_finite(row, "average_temperature_c") for row in thermal]
    coupling_errors = [
        abs(source - scale * loss) / max(abs(scale * loss), 1.0e-300)
        for loss, source in zip(magnetic_values, heat[1:])
    ]
    time_steps = [b - a for a, b in zip(times, times[1:])]

    checks = {
        "loss_components_nonnegative": all(
            _finite(row, name) >= 0.0
            for row in magnetic
            for name in ("target_loss_w", "auxiliary_loss_w", "other_loss_w", "total_loss_w")
        ),
        "loss_decomposition_closes": max(decomposition_errors) <= float(decomposition_rtol),
        "active_power_positive_and_covers_loss": all(
            math.isfinite(value) and float(minimum_power_coverage) <= value <= 1.0 + float(decomposition_rtol)
            for value in power_coverages
        ),
        "thermal_time_starts_at_zero": abs(times[0]) <= 1.0e-15,
        "thermal_time_strictly_increasing": all(step > 0.0 for step in time_steps),
        "thermal_time_step_uniform": max(time_steps) - min(time_steps) <= 1.0e-12,
        "initial_heat_source_zero": abs(heat[0]) <= 1.0e-12,
        "coupled_heat_source_positive": all(value > 0.0 for value in heat[1:]),
        "loss_to_heat_scale_matches": max(coupling_errors) <= float(coupling_rtol),
        "initial_temperature_matches": abs(temperatures[0] - float(initial_temperature_c)) <= float(initial_temperature_atol_c),
        "average_temperature_strictly_increases": all(a < b for a, b in zip(temperatures, temperatures[1:])),
    }
    return {
        "policy": "loss_temperature_coupling_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "magnetic_row_count": len(magnetic),
            "thermal_row_count": len(thermal),
            "loss_to_heat_scale": scale,
            "maximum_coupling_relative_error": max(coupling_errors),
            "maximum_loss_decomposition_relative_error": max(decomposition_errors),
            "power_coverage_range": [min(power_coverages), max(power_coverages)],
            "time_step_s": time_steps[0],
            "final_temperature_c": temperatures[-1],
            "temperature_rise_c": temperatures[-1] - temperatures[0],
        },
        "lesson": (
            "A magnetic-to-thermal handoff must bind each loss row to one heat-source row with an explicit "
            "scale, preserve power decomposition, and verify the thermal time axis and temperature response. "
            "Do not compare raw loss and heat values before accounting for duty, symmetry, or depth scaling."
        ),
    }
