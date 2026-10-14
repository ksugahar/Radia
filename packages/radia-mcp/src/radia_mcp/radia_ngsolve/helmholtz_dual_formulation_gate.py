"""Solver-neutral Helmholtz-coil dual-formulation axis gate."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    if positive and parsed <= 0.0:
        raise ValueError(f"{name} must be positive")
    return parsed


def _array(value: object, name: str, expected: int | None = None) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    parsed = [_finite(item, f"{name}[{index}]") for index, item in enumerate(value)]
    if expected is not None and len(parsed) != expected:
        raise ValueError(f"{name} must contain {expected} values")
    return parsed


def helmholtz_dual_formulation_axis_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate symmetry, central flatness, and agreement of two field formulations."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    units = summary.get("units")
    tolerances = summary.get("gate_tolerances")
    if not isinstance(units, Mapping) or not isinstance(tolerances, Mapping):
        raise ValueError("units and gate_tolerances must be objects")

    axis = _array(summary.get("axis_m"), "axis_m")
    if len(axis) < 9 or len(axis) % 2 != 1:
        raise ValueError("axis_m must contain an odd number of at least nine samples")
    if any(right <= left for left, right in zip(axis, axis[1:])):
        raise ValueError("axis_m must be strictly increasing")
    count = len(axis)
    primary = _array(summary.get("primary_field_T"), "primary_field_T", count)
    secondary = _array(summary.get("secondary_field_T"), "secondary_field_T", count)
    primary_gradient = _array(
        summary.get("primary_gradient_T_per_m"), "primary_gradient_T_per_m", count
    )
    secondary_gradient = _array(
        summary.get("secondary_gradient_T_per_m"), "secondary_gradient_T_per_m", count
    )

    center = count // 2
    spacing = axis[center + 1] - axis[center]
    spacing_error = max(
        abs((right - left) - spacing) for left, right in zip(axis, axis[1:])
    ) / spacing
    axis_symmetry = max(abs(left + right) for left, right in zip(axis, reversed(axis))) / max(
        axis[-1] - axis[0], spacing
    )
    primary_scale = abs(primary[center])
    secondary_scale = abs(secondary[center])
    if primary_scale == 0.0 or secondary_scale == 0.0:
        raise ValueError("both center fields must be nonzero")
    field_scale = max(max(abs(value) for value in primary), primary_scale)
    gradient_scale = max(max(abs(value) for value in primary_gradient), 1.0e-300)

    primary_symmetry = max(
        abs(left - right) for left, right in zip(primary, reversed(primary))
    ) / primary_scale
    secondary_symmetry = max(
        abs(left - right) for left, right in zip(secondary, reversed(secondary))
    ) / secondary_scale
    field_agreement = max(
        abs(left - right) for left, right in zip(primary, secondary)
    ) / field_scale
    gradient_agreement = max(
        abs(left - right) for left, right in zip(primary_gradient, secondary_gradient)
    ) / gradient_scale
    primary_gradient_odd_error = max(
        abs(left + right)
        for left, right in zip(primary_gradient, reversed(primary_gradient))
    ) / gradient_scale
    secondary_gradient_scale = max(
        max(abs(value) for value in secondary_gradient), 1.0e-300
    )
    secondary_gradient_odd_error = max(
        abs(left + right)
        for left, right in zip(secondary_gradient, reversed(secondary_gradient))
    ) / secondary_gradient_scale
    center_gradient = max(
        abs(primary_gradient[center]) * spacing / primary_scale,
        abs(secondary_gradient[center]) * spacing / secondary_scale,
    )
    primary_curvature = abs(
        primary[center + 1] - 2.0 * primary[center] + primary[center - 1]
    ) / primary_scale
    secondary_curvature = abs(
        secondary[center + 1] - 2.0 * secondary[center] + secondary[center - 1]
    ) / secondary_scale
    edge_to_center = max(abs(primary[0]), abs(primary[-1])) / primary_scale

    limits = {
        "symmetry": _finite(
            tolerances.get("maximum_field_symmetry_relative"),
            "maximum_field_symmetry_relative",
            positive=True,
        ),
        "field_agreement": _finite(
            tolerances.get("maximum_formulation_field_relative_error"),
            "maximum_formulation_field_relative_error",
            positive=True,
        ),
        "gradient_agreement": _finite(
            tolerances.get("maximum_formulation_gradient_relative_error"),
            "maximum_formulation_gradient_relative_error",
            positive=True,
        ),
        "center_gradient": _finite(
            tolerances.get("maximum_center_gradient_normalized"),
            "maximum_center_gradient_normalized",
            positive=True,
        ),
        "center_curvature": _finite(
            tolerances.get("maximum_center_curvature_normalized"),
            "maximum_center_curvature_normalized",
            positive=True,
        ),
        "gradient_odd": _finite(
            tolerances.get("maximum_gradient_odd_symmetry_relative"),
            "maximum_gradient_odd_symmetry_relative",
            positive=True,
        ),
        "edge_ratio": _finite(
            tolerances.get("maximum_edge_to_center_abs_ratio"),
            "maximum_edge_to_center_abs_ratio",
            positive=True,
        ),
    }
    checks = {
        "units_are_explicit": units.get("axis") == "m"
        and units.get("field") == "T"
        and units.get("gradient") == "T/m",
        "axis_is_uniform_and_centered": spacing_error <= 1.0e-12
        and axis_symmetry <= 1.0e-12
        and abs(axis[center]) <= 1.0e-15,
        "center_fields_have_same_orientation": primary[center] * secondary[center] > 0.0,
        "both_fields_are_mirror_symmetric": max(primary_symmetry, secondary_symmetry)
        <= limits["symmetry"],
        "field_formulations_agree": field_agreement <= limits["field_agreement"],
        "gradient_formulations_agree": gradient_agreement
        <= limits["gradient_agreement"],
        "gradients_are_approximately_odd": max(
            primary_gradient_odd_error, secondary_gradient_odd_error
        )
        <= limits["gradient_odd"],
        "center_gradient_is_small": center_gradient <= limits["center_gradient"],
        "center_curvature_is_small": max(primary_curvature, secondary_curvature)
        <= limits["center_curvature"],
        "axis_profile_is_nontrivial_and_decays": edge_to_center
        <= limits["edge_ratio"],
    }
    return {
        "policy": "helmholtz_dual_formulation_axis_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": count,
            "axis_spacing_m": spacing,
            "primary_center_field_T": primary[center],
            "secondary_center_field_T": secondary[center],
            "primary_field_symmetry_relative": primary_symmetry,
            "secondary_field_symmetry_relative": secondary_symmetry,
            "formulation_field_relative_error": field_agreement,
            "formulation_gradient_relative_error": gradient_agreement,
            "primary_gradient_odd_symmetry_relative": primary_gradient_odd_error,
            "secondary_gradient_odd_symmetry_relative": secondary_gradient_odd_error,
            "center_gradient_normalized": center_gradient,
            "primary_center_curvature_normalized": primary_curvature,
            "secondary_center_curvature_normalized": secondary_curvature,
            "edge_to_center_abs_ratio": edge_to_center,
        },
        "lesson": (
            "A Helmholtz-coil axis result should be tested as a field profile: two formulations must agree, "
            "the field must be even, its gradient approximately odd, and the central gradient and curvature "
            "small without collapsing to a constant or zero profile."
        ),
    }
