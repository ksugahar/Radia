"""Validate gauge-invariant fields and losses in frequency-domain magnetics."""

from __future__ import annotations

import cmath
import json
import math
from collections.abc import Mapping, Sequence


def _number(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _complex_matrix(real: object, imag: object, name: str) -> list[list[complex]]:
    if not isinstance(real, Sequence) or isinstance(real, (str, bytes)):
        raise ValueError(f"{name}.real must be a matrix")
    if not isinstance(imag, Sequence) or isinstance(imag, (str, bytes)):
        raise ValueError(f"{name}.imag must be a matrix")
    if len(real) != 3 or len(imag) != 3:
        raise ValueError(f"{name} must contain three vector components")
    rows: list[list[complex]] = []
    for component, (real_row, imag_row) in enumerate(zip(real, imag)):
        if not isinstance(real_row, Sequence) or isinstance(real_row, (str, bytes)):
            raise ValueError(f"{name}.real[{component}] must be an array")
        if not isinstance(imag_row, Sequence) or isinstance(imag_row, (str, bytes)):
            raise ValueError(f"{name}.imag[{component}] must be an array")
        if len(real_row) != len(imag_row) or not real_row:
            raise ValueError(f"{name} component lengths must match and be non-empty")
        rows.append([
            complex(_number(r, name), _number(i, name))
            for r, i in zip(real_row, imag_row)
        ])
    if len({len(row) for row in rows}) != 1:
        raise ValueError(f"{name} component sample counts must match")
    return rows


def evaluate_gauge_invariance(summary: Mapping[str, object]) -> dict:
    """Gate physical-field/loss invariance while leaving potential ungated."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    units = summary.get("units")
    tolerances = summary.get("tolerances")
    if not isinstance(units, Mapping) or not isinstance(tolerances, Mapping):
        raise ValueError("units and tolerances must be objects")
    without = summary.get("b_without_gauge")
    with_gauge = summary.get("b_with_gauge")
    if not isinstance(without, Mapping) or not isinstance(with_gauge, Mapping):
        raise ValueError("both magnetic-field records must be objects")
    field_without = _complex_matrix(without.get("real"), without.get("imag"), "b_without_gauge")
    field_with = _complex_matrix(with_gauge.get("real"), with_gauge.get("imag"), "b_with_gauge")
    if len(field_without[0]) != len(field_with[0]):
        raise ValueError("gauge variants must use the same sample count")

    frequency = _number(summary.get("frequency_hz"), "frequency_hz")
    radius = _number(summary.get("radius_m"), "radius_m")
    skin_depth = _number(summary.get("skin_depth_m"), "skin_depth_m")
    applied_field = _number(summary.get("applied_flux_density_T"), "applied_flux_density_T")
    sphere_loss_without = _number(
        summary.get("sphere_loss_without_gauge_W"), "sphere_loss_without_gauge_W"
    )
    sphere_loss_with = _number(
        summary.get("sphere_loss_with_gauge_W"), "sphere_loss_with_gauge_W"
    )
    air_loss_without = abs(
        _number(summary.get("air_loss_without_gauge_W"), "air_loss_without_gauge_W")
    )
    air_loss_with = abs(
        _number(summary.get("air_loss_with_gauge_W"), "air_loss_with_gauge_W")
    )
    air_conductivity = _number(
        summary.get("air_conductivity_S_per_m"), "air_conductivity_S_per_m"
    )
    if min(frequency, radius, skin_depth, applied_field, sphere_loss_without, sphere_loss_with) <= 0:
        raise ValueError("frequency, lengths, field, and conductor losses must be positive")
    if air_conductivity < 0:
        raise ValueError("air conductivity must be nonnegative")

    sample_errors = []
    sample_magnitudes = []
    for sample in range(len(field_without[0])):
        left = [field_without[component][sample] for component in range(3)]
        right = [field_with[component][sample] for component in range(3)]
        left_norm = math.sqrt(sum(abs(value) ** 2 for value in left))
        right_norm = math.sqrt(sum(abs(value) ** 2 for value in right))
        difference = math.sqrt(sum(abs(a - b) ** 2 for a, b in zip(left, right)))
        sample_errors.append(difference / max(left_norm, right_norm, 1.0e-300))
        sample_magnitudes.append(0.5 * (left_norm + right_norm))
    field_relative_error = max(sample_errors)
    loss_relative_error = abs(sphere_loss_without - sphere_loss_with) / max(
        sphere_loss_without, sphere_loss_with
    )
    air_loss_reduction_ratio = air_loss_with / max(air_loss_without, 1.0e-300)
    skin_depth_to_radius = skin_depth / radius
    internal_field_amplification = sample_magnitudes[0] / applied_field
    dominant_phase_lag_rad = abs(cmath.phase(field_with[0][0]))

    field_tolerance = _number(
        tolerances.get("magnetic_field_relative_error"),
        "tolerances.magnetic_field_relative_error",
    )
    loss_tolerance = _number(
        tolerances.get("conductor_loss_relative_error"),
        "tolerances.conductor_loss_relative_error",
    )
    air_reduction_tolerance = _number(
        tolerances.get("air_loss_reduction_ratio"),
        "tolerances.air_loss_reduction_ratio",
    )
    phase_tolerance = _number(
        tolerances.get("weak_skin_phase_lag_rad"),
        "tolerances.weak_skin_phase_lag_rad",
    )
    checks = {
        "units_explicit": units.get("magnetic_flux_density") == "T"
        and units.get("length") == "m"
        and units.get("power") == "W"
        and units.get("frequency") == "Hz",
        "frequency_domain_declared": summary.get("analysis") == "frequency_domain",
        "same_physical_samples_compared": len(sample_errors) >= 3,
        "magnetic_field_is_gauge_invariant": field_relative_error <= field_tolerance,
        "conductor_loss_is_gauge_invariant": loss_relative_error <= loss_tolerance,
        "zero_conductivity_air_declared": air_conductivity == 0.0,
        "gauge_fixing_suppresses_air_loss_artifact": air_loss_reduction_ratio
        <= air_reduction_tolerance,
        "weak_skin_regime_identified": skin_depth_to_radius >= 1.0,
        "internal_field_amplification_is_physical": 1.0 < internal_field_amplification <= 3.1,
        "weak_skin_phase_lag_is_small": dominant_phase_lag_rad <= phase_tolerance,
        "vector_potential_is_not_used_as_invariant": summary.get(
            "potential_comparison_policy"
        )
        == "not_gated_gauge_dependent",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "schema": "radia-differential-forms-gauge-invariance/v1",
        "policy": "physical_fields_and_losses_invariant_potential_not_gated_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(sample_errors),
            "magnetic_field_max_relative_error": field_relative_error,
            "conductor_loss_relative_error": loss_relative_error,
            "air_loss_reduction_ratio": air_loss_reduction_ratio,
            "skin_depth_to_radius": skin_depth_to_radius,
            "internal_field_amplification": internal_field_amplification,
            "dominant_phase_lag_rad": dominant_phase_lag_rad,
        },
        "notes": [
            "Gauge fixing may change the vector potential but must not change magnetic flux density or physical conductor loss.",
            "For a high-permeability sphere in the weak-skin regime, the internal field approaches three times the applied field with a small phase lag.",
            "A zero-conductivity exterior should not be credited with physical Joule loss; use its residual only as a numerical diagnostic.",
        ],
    }


def gauge_invariance_gate(summary_json: str) -> str:
    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    return json.dumps(evaluate_gauge_invariance(summary), indent=2, sort_keys=True)
