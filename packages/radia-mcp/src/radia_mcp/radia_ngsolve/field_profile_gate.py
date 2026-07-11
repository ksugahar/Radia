"""Solver-neutral gates for symmetric field-profile cross-validation."""

from __future__ import annotations

import math
from typing import Any


def symmetric_axial_field_profile_gate(
    axis_positions: list[float],
    axial_field: list[float],
    *,
    expected_center_field: float,
    transverse_field_1: list[float] | None = None,
    transverse_field_2: list[float] | None = None,
    min_sample_count: int = 5,
    max_center_relative_error: float = 1.0e-6,
    max_symmetry_relative: float = 1.0e-9,
    max_transverse_relative: float = 1.0e-9,
    max_axis_symmetry_relative: float = 1.0e-9,
    monotonic_relative_slack: float = 1.0e-12,
) -> dict[str, Any]:
    """Gate an odd, origin-centered axial field profile against an analytic value.

    The gate is solver-neutral. It checks the whole sampled profile rather than
    accepting a center-point match alone.
    """

    positions = [float(value) for value in axis_positions]
    axial = [float(value) for value in axial_field]
    transverse_1 = (
        [float(value) for value in transverse_field_1]
        if transverse_field_1 is not None
        else [0.0] * len(positions)
    )
    transverse_2 = (
        [float(value) for value in transverse_field_2]
        if transverse_field_2 is not None
        else [0.0] * len(positions)
    )
    expected = float(expected_center_field)
    tolerances = {
        "max_center_relative_error": float(max_center_relative_error),
        "max_symmetry_relative": float(max_symmetry_relative),
        "max_transverse_relative": float(max_transverse_relative),
        "max_axis_symmetry_relative": float(max_axis_symmetry_relative),
        "monotonic_relative_slack": float(monotonic_relative_slack),
    }
    if int(min_sample_count) < 5:
        raise ValueError("min_sample_count must be at least 5")
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances.values()):
        raise ValueError("tolerances must be finite and nonnegative")
    if len(positions) != len(axial):
        raise ValueError("axis_positions and axial_field must have equal length")
    if len(transverse_1) != len(positions) or len(transverse_2) != len(positions):
        raise ValueError("transverse field arrays must match axis_positions")
    if not all(math.isfinite(value) for values in (positions, axial, transverse_1, transverse_2) for value in values):
        raise ValueError("profile values must be finite")
    if not math.isfinite(expected) or expected == 0.0:
        raise ValueError("expected_center_field must be finite and nonzero")

    count = len(positions)
    center_index = count // 2
    field_scale = max((abs(value) for value in axial), default=0.0)
    axis_scale = max((abs(value) for value in positions), default=0.0)
    center_value = axial[center_index] if count else math.nan
    center_relative_error = abs(center_value - expected) / abs(expected) if count else math.inf
    symmetry_relative = (
        max((abs(axial[i] - axial[-1 - i]) for i in range(count // 2)), default=0.0)
        / field_scale
        if field_scale > 0.0
        else math.inf
    )
    axis_symmetry_relative = (
        max((abs(positions[i] + positions[-1 - i]) for i in range(count // 2)), default=0.0)
        / axis_scale
        if axis_scale > 0.0
        else math.inf
    )
    transverse_relative = (
        max((abs(value) for value in transverse_1 + transverse_2), default=0.0) / field_scale
        if field_scale > 0.0
        else math.inf
    )
    monotonic_slack = field_scale * tolerances["monotonic_relative_slack"]
    increasing_to_center = all(
        axial[index + 1] + monotonic_slack >= axial[index]
        for index in range(center_index)
    )
    decreasing_from_center = all(
        axial[index + 1] <= axial[index] + monotonic_slack
        for index in range(center_index, max(center_index, count - 1))
    )
    strictly_increasing_axis = all(
        positions[index + 1] > positions[index] for index in range(max(0, count - 1))
    )

    checks = {
        "odd_sample_count_sufficient": count >= int(min_sample_count) and count % 2 == 1,
        "axis_strictly_increasing": strictly_increasing_axis,
        "axis_straddles_origin": count > 0 and positions[0] < 0.0 < positions[-1],
        "center_sample_at_origin": (
            count > 0
            and axis_scale > 0.0
            and abs(positions[center_index]) <= tolerances["max_axis_symmetry_relative"] * axis_scale
        ),
        "axis_is_antisymmetric": axis_symmetry_relative <= tolerances["max_axis_symmetry_relative"],
        "axial_field_nonzero": field_scale > 0.0,
        "center_matches_analytic_value": center_relative_error <= tolerances["max_center_relative_error"],
        "axial_field_is_symmetric": symmetry_relative <= tolerances["max_symmetry_relative"],
        "transverse_field_is_negligible": transverse_relative <= tolerances["max_transverse_relative"],
        "field_increases_toward_center": increasing_to_center,
        "field_decreases_from_center": decreasing_from_center,
    }
    return {
        "policy": "symmetric_axial_field_profile_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": count,
            "center_index": center_index,
            "center_axis": positions[center_index] if count else None,
            "center_field": center_value if count else None,
            "expected_center_field": expected,
            "center_relative_error": center_relative_error,
            "axis_symmetry_relative": axis_symmetry_relative,
            "field_symmetry_relative": symmetry_relative,
            "transverse_relative": transverse_relative,
        },
        "tolerances": {"min_sample_count": int(min_sample_count), **tolerances},
        "notes": [
            "center agreement alone is insufficient; symmetry and both half-profiles are gated",
            "transverse components must be negligible for an axial symmetry-line comparison",
        ],
    }


def dual_formulation_symmetric_field_profile_gate(
    summary: dict[str, Any],
    *,
    max_profile_relative_difference: float = 0.01,
    max_center_relative_difference: float = 0.01,
    max_symmetry_relative: float = 0.01,
    min_sample_count: int = 21,
) -> dict[str, Any]:
    """Require two formulations to agree on a nonzero symmetric profile."""

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    tolerances = (
        max_profile_relative_difference,
        max_center_relative_difference,
        max_symmetry_relative,
    )
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("relative tolerances must be finite and nonnegative")
    if int(min_sample_count) < 3:
        raise ValueError("min_sample_count must be at least 3")

    def finite_number(name: str) -> float | None:
        try:
            value = float(summary.get(name))
        except (TypeError, ValueError):
            return None
        return value if math.isfinite(value) else None

    sample_count = finite_number("sample_count")
    axis_min = finite_number("axis_min")
    axis_max = finite_number("axis_max")
    center_axis = finite_number("center_axis")
    center_a = finite_number("center_value_a")
    center_b = finite_number("center_value_b")
    profile_difference = finite_number("profile_relative_l2_difference")
    center_difference = finite_number("center_relative_difference")
    symmetry_a = finite_number("symmetry_relative_a")
    symmetry_b = finite_number("symmetry_relative_b")
    field_scale = max(abs(center_a or 0.0), abs(center_b or 0.0))
    axis_span = (axis_max - axis_min) if axis_min is not None and axis_max is not None else None

    checks = {
        "component_id_recorded": bool(str(summary.get("component_id") or "").strip()),
        "field_unit_recorded": bool(str(summary.get("field_unit") or "").strip()),
        "axis_unit_recorded": bool(str(summary.get("axis_unit") or "").strip()),
        "sample_count_sufficient": sample_count is not None and sample_count >= int(min_sample_count),
        "axis_straddles_center": axis_min is not None and axis_max is not None and axis_min < 0.0 < axis_max,
        "center_sample_near_origin": (
            center_axis is not None and axis_span is not None and axis_span > 0.0
            and abs(center_axis) <= max(1e-12, 1e-6 * axis_span)
        ),
        "center_field_nonzero": field_scale > 0.0,
        "profile_formulations_agree": (
            profile_difference is not None
            and profile_difference <= float(max_profile_relative_difference)
        ),
        "center_formulations_agree": (
            center_difference is not None
            and center_difference <= float(max_center_relative_difference)
        ),
        "formulation_a_is_symmetric": (
            symmetry_a is not None and symmetry_a <= float(max_symmetry_relative)
        ),
        "formulation_b_is_symmetric": (
            symmetry_b is not None and symmetry_b <= float(max_symmetry_relative)
        ),
    }
    return {
        "policy": "dual_formulation_symmetric_field_profile_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": sample_count,
            "axis_min": axis_min,
            "axis_max": axis_max,
            "center_axis": center_axis,
            "center_value_a": center_a,
            "center_value_b": center_b,
            "profile_relative_l2_difference": profile_difference,
            "center_relative_difference": center_difference,
            "symmetry_relative_a": symmetry_a,
            "symmetry_relative_b": symmetry_b,
        },
        "tolerances": {
            "max_profile_relative_difference": float(max_profile_relative_difference),
            "max_center_relative_difference": float(max_center_relative_difference),
            "max_symmetry_relative": float(max_symmetry_relative),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "agreement at one point is insufficient; compare the full profile and center",
            "each formulation must independently satisfy the expected reflection symmetry",
            "a zero field cannot pass a relative-agreement gate",
        ],
    }
