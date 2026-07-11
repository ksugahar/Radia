"""Solver-independent history gate for irreversible permanent-magnet demagnetization."""
from __future__ import annotations

import json
import math


def permanent_magnet_demagnetization_history_gate(
    summary_json: str,
    state_tolerance: float = 1.0e-9,
    minimum_damage_fraction: float = 1.0e-3,
) -> dict:
    """Check that a load excursion leaves a bounded, persistent magnet state."""

    if state_tolerance < 0.0 or minimum_damage_fraction <= 0.0:
        raise ValueError("state_tolerance must be nonnegative and minimum_damage_fraction positive")
    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    if not isinstance(summary, dict):
        raise ValueError("summary_json must decode to an object")

    steps = summary.get("steps")
    if not isinstance(steps, list) or len(steps) < 3:
        raise ValueError("steps must contain at least initial, stress, and recovery records")
    stress_index = int(summary.get("stress_step_index", -1))
    recovery_index = int(summary.get("recovery_step_index", -1))
    if not (0 < stress_index < recovery_index < len(steps)):
        raise ValueError("stress_step_index and recovery_step_index must identify an ordered load history")

    vectors: list[list[float]] = []
    peak_fields: list[float] = []
    step_ids: list[int] = []
    expected_size: int | None = None
    for position, step in enumerate(steps):
        if not isinstance(step, dict):
            raise ValueError(f"steps[{position}] must be an object")
        step_ids.append(int(step.get("step", position)))
        raw_state = step.get("magnet_state_fraction")
        if not isinstance(raw_state, list) or not raw_state:
            raise ValueError(f"steps[{position}].magnet_state_fraction must be a nonempty list")
        state = [float(value) for value in raw_state]
        if not all(math.isfinite(value) for value in state):
            raise ValueError("magnet state fractions must be finite")
        if expected_size is None:
            expected_size = len(state)
        elif len(state) != expected_size:
            raise ValueError("all state vectors must use the same magnet partition")
        vectors.append(state)
        peak = float(step.get("peak_flux_density_T", 0.0))
        if not math.isfinite(peak) or peak < 0.0:
            raise ValueError("peak_flux_density_T must be finite and nonnegative")
        peak_fields.append(peak)

    initial = vectors[0]
    stressed = vectors[stress_index]
    recovered = vectors[recovery_index]
    state_increases = [
        later[element] - earlier[element]
        for earlier, later in zip(vectors, vectors[1:])
        for element in range(len(initial))
    ]
    damage = [before - after for before, after in zip(initial, stressed)]
    recovery_drift = [after - during for during, after in zip(stressed, recovered)]
    irreversible_elements = sum(value > minimum_damage_fraction for value in damage)
    stress_to_recovery_drift = max(
        abs(vectors[position][element] - stressed[element])
        for position in range(stress_index + 1, recovery_index + 1)
        for element in range(len(initial))
    )

    checks = {
        "state_unit_is_fraction": summary.get("state_unit") == "fraction_of_reference_remanence",
        "field_unit_is_T": summary.get("field_unit") == "T",
        "step_ids_strictly_increase": all(a < b for a, b in zip(step_ids, step_ids[1:])),
        "state_fractions_bounded": all(
            -state_tolerance <= value <= 1.0 + state_tolerance
            for vector in vectors
            for value in vector
        ),
        "state_never_spontaneously_recovers": max(state_increases, default=0.0) <= state_tolerance,
        "stress_causes_resolved_damage": irreversible_elements > 0,
        "recovery_retains_stressed_state": max((abs(value) for value in recovery_drift), default=0.0) <= state_tolerance,
        "stress_to_recovery_state_is_history_locked": stress_to_recovery_drift <= state_tolerance,
        "field_changes_between_stress_and_recovery": (
            max(peak_fields[stress_index : recovery_index + 1])
            - min(peak_fields[stress_index : recovery_index + 1])
            > 0.0
        ),
    }
    return {
        "policy": "permanent_magnet_demagnetization_history_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "partition_count": len(initial),
            "irreversibly_damaged_partition_count": irreversible_elements,
            "maximum_damage_fraction": max(damage),
            "minimum_stressed_state_fraction": min(stressed),
            "maximum_state_increase": max(state_increases, default=0.0),
            "maximum_recovery_state_drift": max((abs(value) for value in recovery_drift), default=0.0),
            "maximum_stress_to_recovery_state_drift": stress_to_recovery_drift,
            "stress_to_recovery_peak_flux_excursion_T": (
                max(peak_fields[stress_index : recovery_index + 1])
                - min(peak_fields[stress_index : recovery_index + 1])
            ),
        },
        "notes": [
            "demagnetization is a path-dependent state check, not a final-field-only threshold",
            "the state vector must retain the worst excursion even when the magnetic field changes after unloading",
        ],
    }
