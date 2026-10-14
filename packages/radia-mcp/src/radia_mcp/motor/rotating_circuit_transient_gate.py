"""Rotating-circuit identities and endpoint-state classification."""

from __future__ import annotations

import math
from typing import Any


def _finite_float(value: Any, name: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _finite_vector(value: Any, name: str, expected: int | None = None) -> list[float]:
    if not isinstance(value, list):
        raise ValueError(f"{name} must be a list")
    parsed = [_finite_float(item, name) for item in value]
    if expected is not None and len(parsed) != expected:
        raise ValueError(f"{name} must contain {expected} values")
    return parsed


def _finite_pairs(value: Any, name: str, expected: int) -> list[list[float]]:
    if not isinstance(value, list) or len(value) != expected:
        raise ValueError(f"{name} must contain {expected} pairs")
    return [_finite_vector(pair, f"{name}[{index}]", 2) for index, pair in enumerate(value)]


def _norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-30)


def rotating_circuit_transient_gate(
    summary: dict[str, Any],
    *,
    max_phase_kcl_relative_error: float = 2.0e-3,
    max_paired_branch_global_relative_error: float = 2.0e-3,
    max_total_power_relative_error: float = 1.0e-12,
    max_kinematic_error_deg: float = 1.0e-8,
    max_cycle_span_error_deg: float = 1.0e-8,
    periodic_state_tolerance: float = 1.0e-2,
    min_sample_count: int = 3,
) -> dict[str, Any]:
    """Validate circuit topology and classify a geometric-cycle endpoint.

    A 0-to-360 degree table is not automatically periodic. The electrical and
    mechanical state must match before a repeated endpoint may be removed for
    an FFT. Transient tables that fail this state check remain valid results,
    but are explicitly classified as not FFT-ready.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("contract") or {}
    rows = summary.get("rows") or []
    if not isinstance(contract, dict):
        raise ValueError("contract must be a mapping")
    if not isinstance(rows, list):
        raise ValueError("rows must be a list")

    phase_count = int(contract.get("phase_count", 0))
    branch_count = int(contract.get("paired_branch_count", 0))
    cycle_deg = _finite_float(contract.get("geometric_cycle_deg", math.nan), "cycle")
    parsed: list[dict[str, Any]] = []
    parse_errors: list[str] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            parse_errors.append(f"row {index} is not a mapping")
            continue
        try:
            parsed.append(
                {
                    "time_s": _finite_float(row["time_s"], f"row {index}.time_s"),
                    "angle_deg": _finite_float(
                        row["angle_deg"], f"row {index}.angle_deg"
                    ),
                    "speed_rpm": _finite_float(
                        row["speed_rpm"], f"row {index}.speed_rpm"
                    ),
                    "torque_nm": _finite_float(
                        row["torque_nm"], f"row {index}.torque_nm"
                    ),
                    "phase_currents": _finite_vector(
                        row["phase_currents_a"],
                        f"row {index}.phase_currents_a",
                        phase_count,
                    ),
                    "phase_flux": _finite_vector(
                        row["phase_flux_linkages_wb"],
                        f"row {index}.phase_flux_linkages_wb",
                        phase_count,
                    ),
                    "branch_currents": _finite_pairs(
                        row["paired_branch_currents_a"],
                        f"row {index}.paired_branch_currents_a",
                        branch_count,
                    ),
                    "branch_powers": _finite_pairs(
                        row["paired_branch_powers_w"],
                        f"row {index}.paired_branch_powers_w",
                        branch_count,
                    ),
                    "power_components": _finite_vector(
                        row["circuit_power_components_w"],
                        f"row {index}.circuit_power_components_w",
                    ),
                    "reported_total_power": _finite_float(
                        row["reported_total_circuit_power_w"],
                        f"row {index}.reported_total_circuit_power_w",
                    ),
                }
            )
        except (KeyError, TypeError, ValueError) as exc:
            parse_errors.append(f"row {index}: {exc}")

    phase_kcl_errors: list[float] = []
    branch_current_residuals: list[float] = []
    branch_current_magnitudes: list[float] = []
    branch_power_residuals: list[float] = []
    branch_power_magnitudes: list[float] = []
    total_power_errors: list[float] = []
    kinematic_errors: list[float] = []
    for row in parsed:
        phase_kcl_errors.append(
            abs(sum(row["phase_currents"]))
            / max(max((abs(value) for value in row["phase_currents"]), default=0.0), 1.0e-30)
        )
        for left, right in row["branch_currents"]:
            branch_current_residuals.append(abs(left + right))
            branch_current_magnitudes.extend((abs(left), abs(right)))
        for left, right in row["branch_powers"]:
            branch_power_residuals.append(abs(left - right))
            branch_power_magnitudes.extend((abs(left), abs(right)))
        total_power_errors.append(
            _relative_difference(row["reported_total_power"], sum(row["power_components"]))
        )
        kinematic_errors.append(
            abs(row["angle_deg"] - 6.0 * row["speed_rpm"] * row["time_s"])
        )

    times = [row["time_s"] for row in parsed]
    angles = [row["angle_deg"] for row in parsed]
    angle_span = angles[-1] - angles[0] if parsed else math.nan
    current_endpoint_mismatch = math.inf
    flux_endpoint_mismatch = math.inf
    torque_endpoint_mismatch = math.inf
    endpoint_state_mismatch = math.inf
    if parsed:
        first = parsed[0]
        last = parsed[-1]
        current_endpoint_mismatch = _norm(
            [right - left for left, right in zip(first["phase_currents"], last["phase_currents"], strict=True)]
        ) / max(_norm(first["phase_currents"]), _norm(last["phase_currents"]), 1.0e-30)
        flux_endpoint_mismatch = _norm(
            [right - left for left, right in zip(first["phase_flux"], last["phase_flux"], strict=True)]
        ) / max(_norm(first["phase_flux"]), _norm(last["phase_flux"]), 1.0e-30)
        torque_endpoint_mismatch = _relative_difference(
            first["torque_nm"], last["torque_nm"]
        )
        endpoint_state_mismatch = max(
            current_endpoint_mismatch,
            flux_endpoint_mismatch,
            torque_endpoint_mismatch,
        )

    endpoint_classification = (
        "periodic_state"
        if endpoint_state_mismatch <= float(periodic_state_tolerance)
        else "nonperiodic_transient"
    )
    current_global_error = max(branch_current_residuals, default=math.inf) / max(
        max(branch_current_magnitudes, default=0.0), 1.0e-30
    )
    power_global_error = max(branch_power_residuals, default=math.inf) / max(
        max(branch_power_magnitudes, default=0.0), 1.0e-30
    )
    checks = {
        "rows_parsed_and_finite": not parse_errors and len(parsed) == len(rows),
        "sample_count_sufficient": len(parsed) >= int(min_sample_count),
        "three_phase_contract": phase_count == 3,
        "paired_branch_contract": branch_count > 0,
        "time_strictly_increases": bool(times)
        and all(right > left for left, right in zip(times, times[1:])),
        "angle_monotone": bool(angles)
        and all(right >= left for left, right in zip(angles, angles[1:])),
        "geometric_cycle_complete": math.isfinite(angle_span)
        and abs(angle_span - cycle_deg) <= float(max_cycle_span_error_deg),
        "angle_speed_time_kinematics": bool(kinematic_errors)
        and max(kinematic_errors) <= float(max_kinematic_error_deg),
        "three_phase_kcl": bool(phase_kcl_errors)
        and max(phase_kcl_errors) <= float(max_phase_kcl_relative_error),
        "paired_branch_currents_antisymmetric": current_global_error
        <= float(max_paired_branch_global_relative_error),
        "paired_branch_powers_symmetric": power_global_error
        <= float(max_paired_branch_global_relative_error),
        "reported_total_power_matches_component_sum": bool(total_power_errors)
        and max(total_power_errors) <= float(max_total_power_relative_error),
        "endpoint_policy_recorded": contract.get("endpoint_periodicity_policy")
        == "require_state_match_before_fft",
        "endpoint_classification_matches": contract.get("expected_endpoint_state")
        == endpoint_classification,
        "fft_readiness_matches_state": bool(contract.get("fft_ready"))
        == (endpoint_classification == "periodic_state"),
    }
    issues = parse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "policy": "rotating_circuit_transient_endpoint_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "endpoint_classification": endpoint_classification,
        "fft_ready": endpoint_classification == "periodic_state",
        "metrics": {
            "sample_count": len(parsed),
            "angle_span_deg": angle_span,
            "max_phase_kcl_relative_error": max(phase_kcl_errors, default=math.inf),
            "paired_branch_current_global_relative_error": current_global_error,
            "paired_branch_power_global_relative_error": power_global_error,
            "max_total_power_relative_error": max(total_power_errors, default=math.inf),
            "max_kinematic_error_deg": max(kinematic_errors, default=math.inf),
            "current_endpoint_relative_mismatch": current_endpoint_mismatch,
            "flux_endpoint_relative_mismatch": flux_endpoint_mismatch,
            "torque_endpoint_relative_mismatch": torque_endpoint_mismatch,
            "endpoint_state_relative_mismatch": endpoint_state_mismatch,
        },
        "tolerances": {
            "max_phase_kcl_relative_error": float(max_phase_kcl_relative_error),
            "max_paired_branch_global_relative_error": float(
                max_paired_branch_global_relative_error
            ),
            "max_total_power_relative_error": float(max_total_power_relative_error),
            "max_kinematic_error_deg": float(max_kinematic_error_deg),
            "max_cycle_span_error_deg": float(max_cycle_span_error_deg),
            "periodic_state_tolerance": float(periodic_state_tolerance),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "Normalize paired-branch residuals by a global branch scale; local relative errors explode at zero crossings.",
            "A geometric 0-to-360-degree endpoint is not a repeated state until currents, flux linkages, and torque agree.",
            "Only periodic_state artifacts are FFT-ready; a nonperiodic transient remains valid time-domain evidence.",
        ],
    }
