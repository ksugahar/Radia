"""Solver-neutral Hartmann channel-flow profile gate."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


_POLICY_MAXIMA = {
    "profile_linf": 2.5e-2,
    "profile_rms": 3.0e-3,
    "symmetry": 1.0e-8,
    "average": 1.0e-4,
    "center": 5.0e-3,
    "field_scaling": 1.0e-10,
}
_POLICY_MIN_SAMPLES = 101


def _finite(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _array(value: object, name: str, count: int | None = None) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    parsed = [_finite(item, f"{name}[{index}]") for index, item in enumerate(value)]
    if count is not None and len(parsed) != count:
        raise ValueError(f"{name} must contain {count} values")
    return parsed


def _limit(source: Mapping[str, object], name: str) -> float:
    value = _finite(source.get(name, _POLICY_MAXIMA[name]), name)
    if value < 0.0 or value > _POLICY_MAXIMA[name]:
        raise ValueError(f"{name} must be between zero and the policy maximum")
    return value


def hartmann_profile_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate a Hartmann-number sweep against the independent channel profile."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    units = summary.get("units")
    tolerances = summary.get("gate_tolerances") or {}
    if not isinstance(units, Mapping) or not isinstance(tolerances, Mapping):
        raise ValueError("units and gate_tolerances must be objects")

    hartmann = _array(summary.get("hartmann_numbers"), "hartmann_numbers")
    if len(hartmann) < 3:
        raise ValueError("hartmann_numbers must contain at least three values")
    count = len(hartmann)
    samples = _array(summary.get("profile_sample_counts"), "profile_sample_counts", count)
    linf = _array(summary.get("profile_max_abs_errors"), "profile_max_abs_errors", count)
    rms = _array(summary.get("profile_rms_errors"), "profile_rms_errors", count)
    symmetry = _array(summary.get("profile_symmetry_errors"), "profile_symmetry_errors", count)
    averages = _array(
        summary.get("normalized_analytic_averages"), "normalized_analytic_averages", count
    )
    center_fem = _array(
        summary.get("normalized_center_velocity_fem"),
        "normalized_center_velocity_fem",
        count,
    )
    center_analytic = _array(
        summary.get("normalized_center_velocity_analytic"),
        "normalized_center_velocity_analytic",
        count,
    )
    boundary_layers = _array(
        summary.get("boundary_layer_fractions"), "boundary_layer_fractions", count
    )
    b0 = _array(summary.get("magnetic_flux_density_T"), "magnetic_flux_density_T", count)

    limits = {name: _limit(tolerances, name) for name in _POLICY_MAXIMA}
    minimum_samples = int(
        _finite(
            tolerances.get("minimum_profile_sample_count", _POLICY_MIN_SAMPLES),
            "minimum_profile_sample_count",
        )
    )
    if minimum_samples < _POLICY_MIN_SAMPLES:
        raise ValueError("minimum_profile_sample_count cannot relax the policy minimum")
    if any(value <= 0.0 for value in hartmann) or any(
        right <= left for left, right in zip(hartmann, hartmann[1:])
    ):
        raise ValueError("hartmann_numbers must be positive and strictly increasing")
    if any(value < 0.0 for values in (linf, rms, symmetry) for value in values):
        raise ValueError("error values must be nonnegative")

    center_relative = [
        abs(fem - analytic) / max(abs(analytic), 1.0e-300)
        for fem, analytic in zip(center_fem, center_analytic)
    ]
    b0_over_ha = [field / ha for field, ha in zip(b0, hartmann)]
    field_scale = sum(b0_over_ha) / count
    field_spread = (
        (max(b0_over_ha) - min(b0_over_ha)) / abs(field_scale)
        if field_scale
        else math.inf
    )
    checks = {
        "units_are_explicit": units.get("hartmann_number") == "1"
        and units.get("normalized_profile") == "1"
        and units.get("magnetic_flux_density") == "T",
        "sweep_has_at_least_three_strictly_increasing_values": count >= 3,
        "profiles_are_well_sampled": all(value >= minimum_samples for value in samples),
        "independent_profile_linf_errors_pass": max(linf) <= limits["profile_linf"],
        "independent_profile_rms_errors_pass": max(rms) <= limits["profile_rms"],
        "profiles_are_mirror_symmetric": max(symmetry) <= limits["symmetry"],
        "analytic_profiles_preserve_unit_average": max(abs(value - 1.0) for value in averages)
        <= limits["average"],
        "center_velocity_matches_independent_profile": max(center_relative)
        <= limits["center"],
        "center_velocity_tends_monotonically_to_plug_flow": all(value > 1.0 for value in center_fem)
        and all(right < left for left, right in zip(center_fem, center_fem[1:])),
        "boundary_layer_thins_monotonically": all(0.0 < value < 1.0 for value in boundary_layers)
        and all(right < left for left, right in zip(boundary_layers, boundary_layers[1:])),
        "magnetic_field_scales_with_hartmann_number": all(value > 0.0 for value in b0)
        and field_spread <= limits["field_scaling"],
    }
    return {
        "policy": "hartmann_profile_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "limits": {**limits, "minimum_profile_sample_count": minimum_samples},
        "metrics": {
            "sweep_count": count,
            "max_profile_linf_error": max(linf),
            "max_profile_rms_error": max(rms),
            "max_symmetry_error": max(symmetry),
            "max_average_error": max(abs(value - 1.0) for value in averages),
            "max_center_velocity_relative_error": max(center_relative),
            "boundary_layer_fractions": boundary_layers,
            "magnetic_flux_density_per_hartmann_T": b0_over_ha,
            "magnetic_field_scaling_relative_spread": field_spread,
        },
        "lesson": (
            "A magnetohydrodynamic channel result should be checked as a parameterized profile: the independent "
            "Hartmann solution keeps unit mean flow, remains symmetric, approaches plug flow, develops a thinner "
            "wall layer as the Hartmann number grows, and uses a magnetic field proportional to that number."
        ),
    }
