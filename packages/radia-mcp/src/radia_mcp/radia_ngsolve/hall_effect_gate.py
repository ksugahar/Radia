"""Solver-neutral constitutive-control gate for Hall-effect sensor sweeps."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _array(value: object, name: str):
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    try:
        result = [float(item) for item in value]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain finite numeric values") from exc
    if (
        len(result) < 9
        or len(result) % 2 == 0
        or not all(math.isfinite(item) for item in result)
    ):
        raise ValueError(f"{name} must contain an odd number of at least nine finite values")
    return result


def _positive_fraction(value: object, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(result) or not 0.0 < result < 1.0:
        raise ValueError(f"{name} must be finite and between zero and one")
    return result


def _l2(values: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def _max_abs(values: Sequence[float]) -> float:
    return max(abs(value) for value in values)


def _diff(values: Sequence[float]) -> list[float]:
    return [right - left for left, right in zip(values, values[1:])]


def _median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return 0.5 * (ordered[middle - 1] + ordered[middle])


def _subtract(left: Sequence[float], right: Sequence[float]) -> list[float]:
    return [a - b for a, b in zip(left, right)]


def _add(left: Sequence[float], right: Sequence[float]) -> list[float]:
    return [a + b for a, b in zip(left, right)]


def _scale(value: float, items: Sequence[float]) -> list[float]:
    return [value * item for item in items]


def _solve_3x3(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    augmented = [row[:] + [rhs[index]] for index, row in enumerate(matrix)]
    for column in range(3):
        pivot = max(range(column, 3), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= 1.0e-30:
            raise ValueError("first-harmonic basis is singular")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(3):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [augmented[row][3] for row in range(3)]


def _first_harmonic_fit(
    angle_deg: Sequence[float], values: Sequence[float]
) -> tuple[list[float], list[float]]:
    rows = [
        [1.0, math.sin(math.radians(angle)), math.cos(math.radians(angle))]
        for angle in angle_deg
    ]
    normal = [[0.0, 0.0, 0.0] for _ in range(3)]
    rhs = [0.0, 0.0, 0.0]
    for row, value in zip(rows, values):
        for i in range(3):
            rhs[i] += row[i] * value
            for j in range(3):
                normal[i][j] += row[i] * row[j]
    coefficients = _solve_3x3(normal, rhs)
    fitted = [
        coefficients[0] + coefficients[1] * row[1] + coefficients[2] * row[2]
        for row in rows
    ]
    return coefficients, fitted


def hall_effect_transverse_voltage_gate(
    summary: Mapping[str, object],
) -> dict[str, Any]:
    """Gate Hall voltage by coefficient, drive, field, and replay controls."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")

    names = (
        "angle_deg",
        "hall_voltage_baseline_v",
        "hall_voltage_replay_v",
        "hall_voltage_zero_coefficient_v",
        "hall_voltage_reversed_coefficient_v",
        "hall_voltage_scaled_drive_v",
        "magnetic_flux_density_baseline_t",
        "magnetic_flux_density_replay_t",
        "magnetic_flux_density_zero_coefficient_t",
        "magnetic_flux_density_reversed_coefficient_t",
        "magnetic_flux_density_scaled_drive_t",
    )
    arrays = {name: _array(summary.get(name), name) for name in names}
    sizes = {len(value) for value in arrays.values()}
    if len(sizes) != 1:
        raise ValueError("all Hall sweep arrays must have equal length")

    angle = arrays["angle_deg"]
    if not all(value > 0.0 for value in _diff(angle)):
        raise ValueError("angle_deg must be strictly increasing")
    drive_ratio = _positive_fraction(summary.get("drive_scale_ratio"), "drive_scale_ratio")
    baseline = arrays["hall_voltage_baseline_v"]
    voltage_norm = max(_l2(baseline), 1.0e-30)
    voltage_scale = max(_max_abs(baseline), 1.0e-30)
    field_baseline = arrays["magnetic_flux_density_baseline_t"]
    field_scale = max(_max_abs(field_baseline), 1.0e-30)

    def relative_l2(residual: Sequence[float]) -> float:
        return _l2(residual) / voltage_norm

    def relative_max(residual: Sequence[float]) -> float:
        return _max_abs(residual) / voltage_scale

    angle_scale = max(_max_abs(angle), 1.0)
    axis_symmetry = (
        _max_abs([a + b for a, b in zip(angle, reversed(angle))]) / angle_scale
    )
    spacing = _diff(angle)
    median_spacing = _median(spacing)
    spacing_error = _max_abs([value - median_spacing for value in spacing]) / max(
        abs(median_spacing), 1.0
    )
    replay_error = relative_l2(_subtract(arrays["hall_voltage_replay_v"], baseline))
    zero_error = relative_max(arrays["hall_voltage_zero_coefficient_v"])
    reversal_error = relative_l2(
        _add(arrays["hall_voltage_reversed_coefficient_v"], baseline)
    )
    drive_error = relative_l2(
        _subtract(arrays["hall_voltage_scaled_drive_v"], _scale(drive_ratio, baseline))
    )
    field_drift = max(
        _max_abs(_subtract(arrays[name], field_baseline)) / field_scale
        for name in (
            "magnetic_flux_density_replay_t",
            "magnetic_flux_density_zero_coefficient_t",
            "magnetic_flux_density_reversed_coefficient_t",
            "magnetic_flux_density_scaled_drive_t",
        )
    )

    coefficients, first_harmonic = _first_harmonic_fit(angle, baseline)
    first_harmonic_residual = _l2(_subtract(baseline, first_harmonic)) / voltage_norm
    dynamic_range = max(baseline) - min(baseline)

    checks = {
        "angle_sweep_is_symmetric_and_uniform": axis_symmetry <= 1.0e-12
        and spacing_error <= 1.0e-12,
        "baseline_has_nonzero_transverse_response": voltage_scale > 1.0e-15
        and dynamic_range / voltage_scale >= 1.0e-3,
        "fresh_voltage_replay_is_deterministic": replay_error <= 1.0e-10,
        "zero_hall_coefficient_suppresses_transverse_voltage": zero_error <= 5.0e-5,
        "hall_coefficient_reversal_reverses_transverse_voltage": reversal_error
        <= 5.0e-4,
        "drive_scaling_is_linear": drive_error <= 1.0e-10,
        "prescribed_magnetic_field_is_case_invariant": field_drift <= 1.0e-12,
    }
    checks = {name: bool(value) for name, value in checks.items()}
    return {
        "policy": "hall_effect_transverse_voltage_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "point_count": len(angle),
            "angle_start_deg": angle[0],
            "angle_stop_deg": angle[-1],
            "angle_axis_symmetry_relative": axis_symmetry,
            "angle_spacing_relative_error": spacing_error,
            "hall_voltage_min_v": min(baseline),
            "hall_voltage_max_v": max(baseline),
            "hall_voltage_dynamic_range_v": dynamic_range,
            "fresh_replay_relative_l2": replay_error,
            "zero_coefficient_relative_max": zero_error,
            "coefficient_sign_reversal_relative_l2": reversal_error,
            "drive_scale_ratio": drive_ratio,
            "drive_scaling_relative_l2": drive_error,
            "magnetic_field_case_relative_max": field_drift,
            "first_harmonic_relative_residual": first_harmonic_residual,
            "first_harmonic_coefficients_v": coefficients,
        },
        "lesson": (
            "Validate a Hall sensor with constitutive controls rather than assuming "
            "that an angle sweep is sinusoidal. Zero Hall coefficient must suppress "
            "the transverse voltage, coefficient reversal must reverse it, drive "
            "scaling must remain linear, and a prescribed one-way magnetic field "
            "must not drift. Report the first-harmonic residual because a finite "
            "near-field source can produce a strongly nonsinusoidal response."
        ),
    }
