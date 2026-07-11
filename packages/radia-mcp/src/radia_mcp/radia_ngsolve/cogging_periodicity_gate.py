"""Solver-neutral validation for one mechanical cogging-torque period."""

from __future__ import annotations

import math
from typing import Any


def _as_positive_integer(value: Any, name: str) -> int:
    parsed = int(value)
    if parsed <= 0 or float(value) != parsed:
        raise ValueError(f"{name} must be a positive integer")
    return parsed


def cogging_torque_periodicity_gate(
    summary: dict[str, Any],
    *,
    max_endpoint_relative_mismatch: float = 0.03,
    max_mean_relative_to_peak: float = 0.05,
    min_peak_to_peak_nm: float = 1.0e-9,
    min_sample_count: int = 9,
) -> dict[str, Any]:
    """Validate a zero-current torque sweep over ``360/lcm(Q, poles)`` degrees.

    A weighted-stress observable is accepted only when the selected body is the
    complete rotor. Selecting one magnet or one material subregion can make the
    stress-weighting mask invalid even though the solver returns a number.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    machine = summary.get("machine") or summary.get("machine_contract") or {}
    if not isinstance(machine, dict):
        raise ValueError("machine must be a mapping")
    slots = _as_positive_integer(machine.get("slots"), "slots")
    poles = _as_positive_integer(machine.get("poles"), "poles")
    expected_order = math.lcm(slots, poles)
    expected_period = 360.0 / expected_order

    currents = [float(value) for value in machine.get("phase_currents_a", [])]
    rows = summary.get("rows") or []
    if not isinstance(rows, list):
        raise ValueError("rows must be a list")
    angles: list[float] = []
    torques: list[float] = []
    parse_errors: list[str] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            parse_errors.append(f"row {index} is not a mapping")
            continue
        try:
            angles.append(float(row["angle_mech_deg"]))
            torque_value = row.get("torque_nm", row.get("weighted_stress_torque_nm"))
            torques.append(float(torque_value))
        except (KeyError, TypeError, ValueError):
            parse_errors.append(f"row {index} lacks finite angle/torque fields")

    finite = (
        len(angles) == len(rows)
        and all(math.isfinite(value) for value in angles + torques + currents)
    )
    increasing = finite and all(right > left for left, right in zip(angles, angles[1:]))
    increments = [right - left for left, right in zip(angles, angles[1:])]
    uniform = bool(increments) and max(increments) - min(increments) <= max(
        1.0e-12, abs(sum(increments) / len(increments)) * 1.0e-9
    )
    covers_expected_period = (
        bool(angles)
        and abs(angles[0]) <= 1.0e-12
        and abs(angles[-1] - expected_period) <= max(1.0e-10, expected_period * 1.0e-9)
    )

    peak = max((abs(value) for value in torques), default=0.0)
    peak_to_peak = max(torques) - min(torques) if torques else 0.0
    endpoint_relative = (
        abs(torques[-1] - torques[0]) / peak if torques and peak > 0.0 else math.inf
    )
    periodic = torques[:-1] if len(torques) > 1 else []
    period_mean = sum(periodic) / len(periodic) if periodic else math.nan
    mean_relative = abs(period_mean) / peak if periodic and peak > 0.0 else math.inf

    harmonics: list[dict[str, float | int]] = []
    for order in range(1, len(periodic) // 2 + 1):
        coefficient = sum(
            value
            * complex(
                math.cos(-2.0 * math.pi * order * index / len(periodic)),
                math.sin(-2.0 * math.pi * order * index / len(periodic)),
            )
            for index, value in enumerate(periodic)
        ) / len(periodic)
        harmonics.append(
            {"period_harmonic": order, "amplitude_nm": 2.0 * abs(coefficient)}
        )
    dominant = max(harmonics, key=lambda row: float(row["amplitude_nm"]), default=None)

    observable = summary.get("torque_observable") or {}
    if not isinstance(observable, dict):
        raise ValueError("torque_observable must be a mapping")
    family = str(observable.get("family") or "").strip()
    selected_body = str(observable.get("selected_body") or "").strip()
    supported_families = {
        "weighted_stress_body_torque",
        "air_gap_maxwell_shear_torque",
        "coenergy_derivative_torque",
    }
    body_contract_ok = family != "weighted_stress_body_torque" or selected_body == "complete_rotor"
    expected_dominant = summary.get("expected_dominant_period_harmonic")
    dominant_ok = dominant is not None and float(dominant["amplitude_nm"]) > 0.0
    if expected_dominant is not None and dominant is not None:
        dominant_ok = dominant_ok and int(dominant["period_harmonic"]) == int(expected_dominant)

    checks = {
        "rows_parsed_and_finite": not parse_errors and finite,
        "sample_count_sufficient": len(rows) >= int(min_sample_count),
        "angles_strictly_increase": increasing,
        "angle_step_uniform": uniform,
        "one_lcm_period_covered": covers_expected_period,
        "zero_current_cogging_condition": bool(currents)
        and all(abs(value) <= 1.0e-12 for value in currents),
        "torque_observable_family_supported": family in supported_families,
        "weighted_stress_selects_complete_rotor": body_contract_ok,
        "torque_nontrivial": peak_to_peak >= float(min_peak_to_peak_nm),
        "periodic_endpoint_closure": endpoint_relative
        <= float(max_endpoint_relative_mismatch),
        "near_zero_period_mean": mean_relative <= float(max_mean_relative_to_peak),
        "dominant_period_harmonic_valid": dominant_ok,
    }
    issues = parse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cogging_torque_periodicity_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "cogging_order_per_revolution": expected_order,
            "expected_period_mech_deg": expected_period,
            "sample_count": len(rows),
            "torque_peak_to_peak_nm": peak_to_peak,
            "torque_period_mean_nm": period_mean,
            "endpoint_relative_mismatch": endpoint_relative,
            "mean_relative_to_peak": mean_relative,
            "dominant_period_harmonic": (
                int(dominant["period_harmonic"]) if dominant else None
            ),
            "harmonics": harmonics,
        },
        "tolerances": {
            "max_endpoint_relative_mismatch": float(max_endpoint_relative_mismatch),
            "max_mean_relative_to_peak": float(max_mean_relative_to_peak),
            "min_peak_to_peak_nm": float(min_peak_to_peak_nm),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "Cogging torque is a zero-current observable over 360/lcm(slots,poles) mechanical degrees.",
            "A weighted-stress mask must enclose the complete moving body, not one material subregion.",
            "Record remeshing drift separately; do not silently replace the measured waveform with a fitted sinusoid.",
        ],
    }
