"""Solver-neutral gates for moving-conductor eddy-brake result tables."""

from __future__ import annotations

import math
from typing import Any


def _relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-30)


def moving_conductor_eddy_brake_gate(
    summary: dict[str, Any],
    *,
    max_kinematic_relative_error: float = 1.0e-9,
    max_decomposition_relative_error: float = 1.0e-12,
    min_sample_count: int = 5,
) -> dict[str, Any]:
    """Validate motion, Lorentz-force, and Joule-loss table identities.

    Mechanical work and Joule heat are reported but are not forced equal when
    the magnetic-energy-rate term is absent from the artifact.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    rows = summary.get("rows") or []
    if not isinstance(rows, list):
        raise ValueError("rows must be a list")
    units = summary.get("units") or {}
    observable = summary.get("observable_contract") or {}
    energy_contract = summary.get("energy_balance_contract") or {}
    if not all(isinstance(value, dict) for value in (units, observable, energy_contract)):
        raise ValueError("units and observable/energy contracts must be mappings")

    parsed: list[dict[str, Any]] = []
    parse_errors: list[str] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            parse_errors.append(f"row {index} is not a mapping")
            continue
        try:
            parsed.append(
                {
                    "time_s": float(row["time_s"]),
                    "displacement_m": float(row["displacement_m"]),
                    "velocity_m_s": float(row["velocity_m_s"]),
                    "lorentz_force_n": float(row["lorentz_force_n"]),
                    "lorentz_force_parts_n": [
                        float(value) for value in row["lorentz_force_parts_n"]
                    ],
                    "joule_loss_w": float(row["joule_loss_w"]),
                    "joule_loss_parts_w": [
                        float(value) for value in row["joule_loss_parts_w"]
                    ],
                }
            )
        except (KeyError, TypeError, ValueError):
            parse_errors.append(f"row {index} is incomplete")

    scalars = [
        value
        for row in parsed
        for key, value in row.items()
        if key not in {"lorentz_force_parts_n", "joule_loss_parts_w"}
    ]
    components = [
        value
        for row in parsed
        for key in ("lorentz_force_parts_n", "joule_loss_parts_w")
        for value in row[key]
    ]
    finite = len(parsed) == len(rows) and all(
        math.isfinite(value) for value in scalars + components
    )
    times = [row["time_s"] for row in parsed]
    increasing = finite and all(right > left for left, right in zip(times, times[1:]))

    kinematic_errors = []
    if increasing and len(parsed) >= 3:
        for index in range(1, len(parsed) - 1):
            velocity_from_displacement = (
                parsed[index + 1]["displacement_m"]
                - parsed[index - 1]["displacement_m"]
            ) / (times[index + 1] - times[index - 1])
            kinematic_errors.append(
                _relative_error(velocity_from_displacement, parsed[index]["velocity_m_s"])
            )
    force_errors = [
        _relative_error(row["lorentz_force_n"], sum(row["lorentz_force_parts_n"]))
        for row in parsed
    ]
    loss_errors = [
        _relative_error(row["joule_loss_w"], sum(row["joule_loss_parts_w"]))
        for row in parsed
    ]

    mechanical_work = 0.0
    joule_energy = 0.0
    if increasing:
        for left, right in zip(parsed, parsed[1:]):
            dt = right["time_s"] - left["time_s"]
            mechanical_work += 0.5 * dt * (
                abs(left["lorentz_force_n"] * left["velocity_m_s"])
                + abs(right["lorentz_force_n"] * right["velocity_m_s"])
            )
            joule_energy += 0.5 * dt * (
                left["joule_loss_w"] + right["joule_loss_w"]
            )
    energy_difference_relative = _relative_error(mechanical_work, joule_energy)

    checks = {
        "rows_parsed_and_finite": not parse_errors and finite,
        "sample_count_sufficient": len(rows) >= int(min_sample_count),
        "time_strictly_increases": increasing,
        "si_units_explicit": units
        == {
            "time": "s",
            "displacement": "m",
            "velocity": "m/s",
            "force": "N",
            "loss": "W",
        },
        "observable_families_explicit": observable.get("force_family")
        == "lorentz_force_on_moving_conductor"
        and observable.get("loss_family") == "joule_loss_in_moving_conductor"
        and observable.get("force_value_kind") in {"signed_component", "absolute_magnitude"},
        "displacement_derivative_matches_velocity": bool(kinematic_errors)
        and max(kinematic_errors) <= float(max_kinematic_relative_error),
        "lorentz_force_decomposition_closes": bool(force_errors)
        and max(force_errors) <= float(max_decomposition_relative_error),
        "joule_loss_decomposition_closes": bool(loss_errors)
        and max(loss_errors) <= float(max_decomposition_relative_error),
        "joule_loss_nonnegative": bool(parsed)
        and all(row["joule_loss_w"] >= 0.0 for row in parsed),
        "force_and_loss_nontrivial": bool(parsed)
        and max(row["lorentz_force_n"] for row in parsed) > 0.0
        and max(row["joule_loss_w"] for row in parsed) > 0.0,
        "missing_magnetic_energy_term_acknowledged": energy_contract.get(
            "magnetic_energy_rate_available"
        )
        is False
        and energy_contract.get("mechanical_work_vs_joule_heat") == "diagnostic_only",
    }
    issues = parse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "policy": "moving_conductor_eddy_brake_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(rows),
            "maximum_kinematic_relative_error": max(kinematic_errors, default=math.inf),
            "maximum_force_decomposition_relative_error": max(force_errors, default=math.inf),
            "maximum_loss_decomposition_relative_error": max(loss_errors, default=math.inf),
            "mechanical_work_magnitude_j": mechanical_work,
            "joule_energy_j": joule_energy,
            "mechanical_joule_relative_difference_diagnostic": energy_difference_relative,
            "peak_lorentz_force_n": max(
                (row["lorentz_force_n"] for row in parsed), default=0.0
            ),
            "peak_joule_loss_w": max(
                (row["joule_loss_w"] for row in parsed), default=0.0
            ),
        },
        "tolerances": {
            "max_kinematic_relative_error": float(max_kinematic_relative_error),
            "max_decomposition_relative_error": float(max_decomposition_relative_error),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "Different force observables can share units but represent different bodies or extraction methods.",
            "Do not assert mechanical work equals Joule heat unless magnetic-energy change and sign conventions are available.",
        ],
    }
