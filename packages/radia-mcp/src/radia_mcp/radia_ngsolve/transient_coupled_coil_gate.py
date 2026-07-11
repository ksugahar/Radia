"""Solver-neutral transient coupled-coil validation."""

from __future__ import annotations

import math
from collections.abc import Iterable


def _finite_series(values: Iterable[float], name: str) -> list[float]:
    try:
        series = [float(value) for value in values]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a numeric sequence") from exc
    if not series or not all(math.isfinite(value) for value in series):
        raise ValueError(f"{name} must contain only finite values")
    return series


def _fit_first_order_response(primary: list[float], secondary: list[float]) -> tuple[float, float, list[float]]:
    previous = secondary[:-1]
    delta_primary = [later - earlier for earlier, later in zip(primary, primary[1:])]
    target = secondary[1:]

    s_xx = sum(value * value for value in previous)
    s_xy = sum(x * y for x, y in zip(previous, delta_primary))
    s_yy = sum(value * value for value in delta_primary)
    s_xz = sum(x * z for x, z in zip(previous, target))
    s_yz = sum(y * z for y, z in zip(delta_primary, target))
    determinant = s_xx * s_yy - s_xy * s_xy
    scale = max(s_xx * s_yy, 1.0)
    if abs(determinant) <= 1.0e-14 * scale:
        raise ValueError("the primary/secondary history does not resolve a two-parameter response")

    memory = (s_xz * s_yy - s_yz * s_xy) / determinant
    coupling = (s_yz * s_xx - s_xz * s_xy) / determinant
    residual = [
        observed - (memory * old + coupling * delta)
        for old, delta, observed in zip(previous, delta_primary, target)
    ]
    return memory, coupling, residual


def transient_coupled_coil_response_gate(
    times_s: Iterable[float],
    primary_current_a: Iterable[float],
    secondary_current_a: Iterable[float],
    *,
    secondary_resistance_ohm: float,
    secondary_turns: float,
    maximum_relative_residual: float = 1.0e-3,
    time_uniformity_rtol: float = 1.0e-9,
    zero_initial_atol_a: float = 1.0e-12,
    minimum_resolved_current_a: float = 1.0e-6,
    minimum_coupling_gain: float = 1.0e-8,
) -> dict[str, object]:
    """Gate a passive shorted-secondary transient response.

    The response is fitted to ``i2[n] = a*i2[n-1] + b*delta(i1[n])``.
    The sign of ``b`` is intentionally unrestricted because it is set by coil
    orientation. A passive first-order decay requires ``0 < a < 1``.
    """

    times = _finite_series(times_s, "times_s")
    primary = _finite_series(primary_current_a, "primary_current_a")
    secondary = _finite_series(secondary_current_a, "secondary_current_a")
    if len(times) < 8 or len(primary) != len(times) or len(secondary) != len(times):
        raise ValueError("time and current histories must have equal lengths of at least eight")

    resistance = float(secondary_resistance_ohm)
    turns = float(secondary_turns)
    residual_limit = float(maximum_relative_residual)
    if not all(math.isfinite(value) for value in (resistance, turns, residual_limit)):
        raise ValueError("resistance, turns, and residual tolerance must be finite")
    if resistance <= 0.0 or turns <= 0.0 or residual_limit <= 0.0:
        raise ValueError("resistance, turns, and residual tolerance must be positive")

    time_steps = [later - earlier for earlier, later in zip(times, times[1:])]
    if any(step <= 0.0 for step in time_steps):
        raise ValueError("times_s must be strictly increasing")
    mean_step = sum(time_steps) / len(time_steps)
    time_spread = max(abs(step - mean_step) for step in time_steps) / mean_step

    memory, coupling, residual = _fit_first_order_response(primary, secondary)
    response_scale = max(max(abs(value) for value in secondary[1:]), 1.0e-300)
    max_relative_residual = max(abs(value) for value in residual) / response_scale
    rms_residual = math.sqrt(sum(value * value for value in residual) / len(residual))
    inferred_time_constant = -mean_step / math.log(memory) if 0.0 < memory < 1.0 else math.nan

    checks = {
        "uniform_time_grid": time_spread <= float(time_uniformity_rtol),
        "zero_initial_primary_current": abs(primary[0]) <= float(zero_initial_atol_a),
        "zero_initial_secondary_current": abs(secondary[0]) <= float(zero_initial_atol_a),
        "positive_secondary_resistance": resistance > 0.0,
        "positive_secondary_turn_count": turns > 0.0,
        "primary_excitation_resolved": max(abs(value) for value in primary) >= float(minimum_resolved_current_a),
        "secondary_response_resolved": response_scale >= float(minimum_resolved_current_a),
        "passive_first_order_memory": 0.0 < memory < 1.0,
        "nonzero_orientation_aware_coupling": abs(coupling) >= float(minimum_coupling_gain),
        "first_order_response_matches": max_relative_residual <= residual_limit,
    }
    return {
        "policy": "transient_coupled_coil_response_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": len(times),
            "time_step_s": mean_step,
            "time_grid_relative_spread": time_spread,
            "memory_factor": memory,
            "coupling_gain": coupling,
            "coupling_polarity": 1 if coupling > 0.0 else -1,
            "inferred_time_constant_s": inferred_time_constant,
            "maximum_relative_residual": max_relative_residual,
            "rms_residual_a": rms_residual,
        },
        "lesson": (
            "A transient coupled-coil result is an induced-current history, not an electromagnetic-force history. "
            "Validate its time grid, passive memory, orientation-aware coupling, and reduced RL recurrence before "
            "using it as a solver or circuit handoff."
        ),
    }
