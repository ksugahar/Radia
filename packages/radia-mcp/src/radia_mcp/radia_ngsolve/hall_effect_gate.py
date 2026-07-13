"""Solver-neutral constitutive-control gate for Hall-effect sensor sweeps."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np


def _array(value: object, name: str) -> np.ndarray:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    result = np.asarray(value, dtype=float).ravel()
    if result.size < 9 or result.size % 2 == 0 or not np.all(np.isfinite(result)):
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
    sizes = {value.size for value in arrays.values()}
    if len(sizes) != 1:
        raise ValueError("all Hall sweep arrays must have equal length")

    angle = arrays["angle_deg"]
    if not np.all(np.diff(angle) > 0.0):
        raise ValueError("angle_deg must be strictly increasing")
    drive_ratio = _positive_fraction(summary.get("drive_scale_ratio"), "drive_scale_ratio")
    baseline = arrays["hall_voltage_baseline_v"]
    voltage_norm = max(float(np.linalg.norm(baseline)), 1.0e-30)
    voltage_scale = max(float(np.max(np.abs(baseline))), 1.0e-30)
    field_baseline = arrays["magnetic_flux_density_baseline_t"]
    field_scale = max(float(np.max(np.abs(field_baseline))), 1.0e-30)

    def relative_l2(residual: np.ndarray) -> float:
        return float(np.linalg.norm(residual) / voltage_norm)

    def relative_max(residual: np.ndarray) -> float:
        return float(np.max(np.abs(residual)) / voltage_scale)

    angle_scale = max(float(np.max(np.abs(angle))), 1.0)
    axis_symmetry = float(np.max(np.abs(angle + angle[::-1])) / angle_scale)
    spacing = np.diff(angle)
    spacing_error = float(
        np.max(np.abs(spacing - np.median(spacing)))
        / max(abs(float(np.median(spacing))), 1.0)
    )
    replay_error = relative_l2(arrays["hall_voltage_replay_v"] - baseline)
    zero_error = relative_max(arrays["hall_voltage_zero_coefficient_v"])
    reversal_error = relative_l2(
        arrays["hall_voltage_reversed_coefficient_v"] + baseline
    )
    drive_error = relative_l2(
        arrays["hall_voltage_scaled_drive_v"] - drive_ratio * baseline
    )
    field_drift = max(
        float(np.max(np.abs(arrays[name] - field_baseline)) / field_scale)
        for name in (
            "magnetic_flux_density_replay_t",
            "magnetic_flux_density_zero_coefficient_t",
            "magnetic_flux_density_reversed_coefficient_t",
            "magnetic_flux_density_scaled_drive_t",
        )
    )

    theta = np.deg2rad(angle)
    first_harmonic_basis = np.column_stack(
        (np.ones_like(theta), np.sin(theta), np.cos(theta))
    )
    coefficients, *_ = np.linalg.lstsq(first_harmonic_basis, baseline, rcond=None)
    first_harmonic = first_harmonic_basis @ coefficients
    first_harmonic_residual = float(
        np.linalg.norm(baseline - first_harmonic) / voltage_norm
    )
    dynamic_range = float(np.ptp(baseline))

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
            "point_count": int(angle.size),
            "angle_start_deg": float(angle[0]),
            "angle_stop_deg": float(angle[-1]),
            "angle_axis_symmetry_relative": axis_symmetry,
            "angle_spacing_relative_error": spacing_error,
            "hall_voltage_min_v": float(np.min(baseline)),
            "hall_voltage_max_v": float(np.max(baseline)),
            "hall_voltage_dynamic_range_v": dynamic_range,
            "fresh_replay_relative_l2": replay_error,
            "zero_coefficient_relative_max": zero_error,
            "coefficient_sign_reversal_relative_l2": reversal_error,
            "drive_scale_ratio": drive_ratio,
            "drive_scaling_relative_l2": drive_error,
            "magnetic_field_case_relative_max": field_drift,
            "first_harmonic_relative_residual": first_harmonic_residual,
            "first_harmonic_coefficients_v": coefficients.tolist(),
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
