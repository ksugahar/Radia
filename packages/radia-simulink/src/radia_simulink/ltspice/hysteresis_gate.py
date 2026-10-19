"""Public-safe cycle-energy validation for hysteretic inductors."""

from __future__ import annotations

import math
from typing import Any


def hysteretic_inductor_cycle_gate(
    cycle_rows: list[dict[str, Any]],
    *,
    expected_current_peak_a: float,
    expected_copper_energy_j: float,
    voltage_thd: float,
    max_current_peak_relative_error: float = 0.01,
    max_copper_energy_relative_error: float = 0.01,
    max_energy_identity_relative_error: float = 0.01,
    max_flux_closure_relative: float = 0.01,
    max_steady_hysteresis_span_relative: float = 0.01,
    min_voltage_thd: float = 0.01,
) -> dict[str, Any]:
    """Gate steady hysteresis loops by energy identity, closure, and repeatability."""

    if not isinstance(cycle_rows, list) or len(cycle_rows) < 2:
        raise ValueError("cycle_rows must contain at least two steady cycles")
    expected_current = float(expected_current_peak_a)
    expected_copper = float(expected_copper_energy_j)
    thd = float(voltage_thd)
    tolerances = [
        float(max_current_peak_relative_error),
        float(max_copper_energy_relative_error),
        float(max_energy_identity_relative_error),
        float(max_flux_closure_relative),
        float(max_steady_hysteresis_span_relative),
        float(min_voltage_thd),
    ]
    if (
        not math.isfinite(expected_current)
        or not math.isfinite(expected_copper)
        or expected_current <= 0.0
        or expected_copper <= 0.0
    ):
        raise ValueError("expected current and copper energy must be finite and positive")
    if not math.isfinite(thd) or thd < 0.0:
        raise ValueError("voltage_thd must be finite and nonnegative")
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    rows = []
    for index, row in enumerate(cycle_rows):
        try:
            normalized = {
                "cycle_index": int(row["cycle_index"]),
                "current_peak_a": float(row["current_peak_a"]),
                "total_energy_j": float(row["total_energy_j"]),
                "copper_energy_j": float(row["copper_energy_j"]),
                "hysteresis_energy_j": float(row["hysteresis_energy_j"]),
                "flux_loop_energy_j": float(row["flux_loop_energy_j"]),
                "flux_closure_relative": float(row["flux_closure_relative"]),
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"cycle row {index} is missing a numeric field") from exc
        if any(not math.isfinite(value) for value in normalized.values()):
            raise ValueError(f"cycle row {index} contains a non-finite value")
        normalized["current_peak_relative_error"] = abs(
            normalized["current_peak_a"] - expected_current
        ) / expected_current
        normalized["copper_energy_relative_error"] = abs(
            normalized["copper_energy_j"] - expected_copper
        ) / expected_copper
        normalized["energy_identity_relative_error"] = abs(
            (normalized["total_energy_j"] - normalized["copper_energy_j"])
            - normalized["flux_loop_energy_j"]
        ) / max(abs(normalized["flux_loop_energy_j"]), math.ulp(1.0))
        rows.append(normalized)

    cycle_indices = [row["cycle_index"] for row in rows]
    hysteresis_energy = [row["hysteresis_energy_j"] for row in rows]
    mean_hysteresis = sum(hysteresis_energy) / len(hysteresis_energy)
    hysteresis_span_relative = (
        (max(hysteresis_energy) - min(hysteresis_energy)) / mean_hysteresis
        if mean_hysteresis > 0.0
        else math.inf
    )
    checks = {
        "cycle_indices_strictly_increase": all(
            a < b for a, b in zip(cycle_indices, cycle_indices[1:])
        ),
        "current_peak_matches_excitation": all(
            row["current_peak_relative_error"] <= max_current_peak_relative_error
            for row in rows
        ),
        "copper_energy_matches_sinusoidal_reference": all(
            row["copper_energy_relative_error"] <= max_copper_energy_relative_error
            for row in rows
        ),
        "total_energy_exceeds_copper_loss": all(
            row["total_energy_j"] > row["copper_energy_j"] > 0.0 for row in rows
        ),
        "hysteresis_energy_is_positive": all(
            row["hysteresis_energy_j"] > 0.0 and row["flux_loop_energy_j"] > 0.0
            for row in rows
        ),
        "terminal_energy_matches_flux_loop_area": all(
            row["energy_identity_relative_error"] <= max_energy_identity_relative_error
            for row in rows
        ),
        "flux_loop_closes_each_steady_cycle": all(
            0.0 <= row["flux_closure_relative"] <= max_flux_closure_relative
            for row in rows
        ),
        "steady_hysteresis_energy_is_repeatable": (
            hysteresis_span_relative <= max_steady_hysteresis_span_relative
        ),
        "voltage_contains_nonlinear_harmonics": thd >= min_voltage_thd,
    }
    return {
        "schema": "radia-spice-lab.hysteretic-inductor-cycle.v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "ok": all(checks.values()),
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "cycle_indices": cycle_indices,
            "hysteresis_energy_j": hysteresis_energy,
            "steady_hysteresis_span_relative": hysteresis_span_relative,
            "max_energy_identity_relative_error": max(
                row["energy_identity_relative_error"] for row in rows
            ),
            "max_flux_closure_relative": max(row["flux_closure_relative"] for row in rows),
            "voltage_thd": thd,
        },
        "rows": rows,
        "notes": [
            "Use settled cycles for closure and repeatability; an initially demagnetized first cycle may establish remanence.",
            "The terminal identity is integral(v*i) - R*integral(i^2) = closed-loop integral(i d lambda).",
        ],
    }
