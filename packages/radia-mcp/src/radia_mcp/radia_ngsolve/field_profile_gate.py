"""Solver-neutral gates for symmetric field-profile cross-validation."""

from __future__ import annotations

import math
from typing import Any


def nonlinear_magnetic_spatial_evidence_gate(
    summary: dict[str, Any],
    *,
    max_average_vector_relative_difference: float = 0.07,
    max_rms_magnitude_relative_difference: float = 0.10,
    min_tensor_gauss_samples: int = 27,
) -> dict[str, Any]:
    """Gate nonlinear magnetic evidence using a matched volume observable.

    A converged fixed-point residual or one field sample is not sufficient. The
    response and material-update polynomial orders must be compatible, and the
    candidate must be compared with an independent result carrying the same
    observable identity, field unit, and coordinate system.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    average_limit = float(max_average_vector_relative_difference)
    rms_limit = float(max_rms_magnitude_relative_difference)
    min_samples = int(min_tensor_gauss_samples)
    if any(
        not math.isfinite(value) or value < 0.0
        for value in (average_limit, rms_limit)
    ):
        raise ValueError("relative tolerances must be finite and nonnegative")
    if min_samples < 1:
        raise ValueError("min_tensor_gauss_samples must be positive")

    def finite_number(payload: dict[str, Any], name: str) -> float | None:
        try:
            value = float(payload.get(name))
        except (TypeError, ValueError):
            return None
        return value if math.isfinite(value) else None

    def finite_vector(payload: dict[str, Any], name: str) -> list[float] | None:
        raw = payload.get(name)
        if not isinstance(raw, (list, tuple)) or len(raw) != 3:
            return None
        try:
            values = [float(value) for value in raw]
        except (TypeError, ValueError):
            return None
        return values if all(math.isfinite(value) for value in values) else None

    try:
        response_order = int(summary.get("response_order"))
        material_order = int(summary.get("material_update_order"))
        sample_count = int(summary.get("sample_count", 0))
        integration_order = int(summary.get("integration_order", 0))
    except (TypeError, ValueError):
        response_order = material_order = sample_count = integration_order = -1
    observable = str(summary.get("spatial_observable") or "").strip()
    average = finite_vector(summary, "average_field_T")
    volume = finite_number(summary, "volume_m3")
    rms = finite_number(summary, "rms_magnitude_T")
    average_magnitude = (
        math.sqrt(sum(value * value for value in average)) if average is not None else None
    )

    reference = summary.get("reference")
    reference = reference if isinstance(reference, dict) else {}
    reference_average = finite_vector(reference, "average_field_T")
    reference_rms = finite_number(reference, "rms_magnitude_T")
    reference_average_magnitude = (
        math.sqrt(sum(value * value for value in reference_average))
        if reference_average is not None
        else None
    )
    average_difference = (
        math.sqrt(
            sum(
                (candidate - expected) ** 2
                for candidate, expected in zip(average, reference_average)
            )
        )
        / reference_average_magnitude
        if average is not None
        and reference_average is not None
        and reference_average_magnitude is not None
        and reference_average_magnitude > 0.0
        else math.inf
    )
    rms_difference = (
        abs(rms - reference_rms) / reference_rms
        if rms is not None and reference_rms is not None and reference_rms > 0.0
        else math.inf
    )
    identity_fields = ("observable_id", "field_unit", "coordinate_system")
    identity_matches = all(
        bool(str(summary.get(name) or "").strip())
        and str(summary.get(name)).strip() == str(reference.get(name) or "").strip()
        for name in identity_fields
    )
    sampling_sufficient = (
        observable == "volume_integral"
        and integration_order >= max(2, 2 * max(response_order, 1))
    ) or (observable == "tensor_gauss" and sample_count >= min_samples)

    checks = {
        "nonlinear_constitutive_result": summary.get("nonlinear") is True,
        "solver_converged": summary.get("solver_converged") is True,
        "not_linear_only_evidence": summary.get("linear_reference_only") is not True,
        "response_material_orders_compatible": (
            response_order >= 1 and material_order >= max(response_order - 1, 0)
        ),
        "spatial_observable_supported": observable in {"volume_integral", "tensor_gauss"},
        "spatial_sampling_sufficient": sampling_sufficient,
        "volume_positive": volume is not None and volume > 0.0,
        "field_finite_nonzero": average_magnitude is not None and average_magnitude > 0.0,
        "rms_consistent_with_average": (
            rms is not None
            and average_magnitude is not None
            and rms + 1.0e-12 * max(rms, average_magnitude, 1.0) >= average_magnitude
        ),
        "reference_identity_matches": identity_matches,
        "average_vector_matches_reference": average_difference <= average_limit,
        "rms_magnitude_matches_reference": rms_difference <= rms_limit,
    }
    return {
        "policy": "nonlinear_magnetic_spatial_evidence_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
        "metrics": {
            "response_order": response_order,
            "material_update_order": material_order,
            "sample_count": sample_count,
            "integration_order": integration_order,
            "volume_m3": volume,
            "average_magnitude_T": average_magnitude,
            "rms_magnitude_T": rms,
            "average_vector_relative_difference": average_difference,
            "rms_magnitude_relative_difference": rms_difference,
        },
        "tolerances": {
            "max_average_vector_relative_difference": average_limit,
            "max_rms_magnitude_relative_difference": rms_limit,
            "min_tensor_gauss_samples": min_samples,
        },
        "notes": [
            "a center-point match and a converged nonlinear iteration are necessary but insufficient",
            "linear-source agreement cannot establish nonlinear constitutive parity",
            "response p-refinement is not evidence when the nonlinear material update remains lower order",
        ],
    }


def symmetric_complex_field_curve_gate(
    axis_positions: list[float],
    field_real: list[float],
    field_imag: list[float] | None = None,
    *,
    axis_unit: str = "m",
    field_unit: str = "A/m",
    log10_relative_residual: float,
    min_sample_count: int = 9,
    max_axis_symmetry_relative: float = 1.0e-9,
    max_field_symmetry_relative: float = 2.0e-3,
    max_log10_relative_residual: float = -8.0,
) -> dict[str, Any]:
    """Gate an origin-centered complex field curve by mirror symmetry.

    Unlike :func:`symmetric_axial_field_profile_gate`, this gate accepts even
    sample counts and does not require an analytic center value.  It is useful
    for result readers that sample a full line without placing a node exactly
    at the origin.
    """

    positions = [float(value) for value in axis_positions]
    real = [float(value) for value in field_real]
    imag = (
        [float(value) for value in field_imag]
        if field_imag is not None
        else [0.0] * len(real)
    )
    count = len(positions)
    if int(min_sample_count) < 5:
        raise ValueError("min_sample_count must be at least 5")
    if len(real) != count or len(imag) != count:
        raise ValueError("axis and field arrays must have equal length")
    if not str(axis_unit).strip() or not str(field_unit).strip():
        raise ValueError("axis_unit and field_unit must be non-empty")
    if not all(
        math.isfinite(value)
        for values in (positions, real, imag)
        for value in values
    ):
        raise ValueError("axis and field values must be finite")
    residual = float(log10_relative_residual)
    axis_tol = float(max_axis_symmetry_relative)
    field_tol = float(max_field_symmetry_relative)
    residual_limit = float(max_log10_relative_residual)
    if not math.isfinite(residual):
        raise ValueError("log10_relative_residual must be finite")
    if any(not math.isfinite(value) or value < 0.0 for value in (axis_tol, field_tol)):
        raise ValueError("relative tolerances must be finite and nonnegative")
    if not math.isfinite(residual_limit) or residual_limit >= 0.0:
        raise ValueError("max_log10_relative_residual must be finite and negative")

    axis_scale = max((abs(value) for value in positions), default=0.0)
    field_scale = max(
        (abs(complex(real[index], imag[index])) for index in range(count)),
        default=0.0,
    )
    pair_count = count // 2
    axis_symmetry_relative = (
        max(
            (abs(positions[index] + positions[-1 - index]) for index in range(pair_count)),
            default=0.0,
        )
        / axis_scale
        if axis_scale > 0.0
        else math.inf
    )
    field_symmetry_relative = (
        max(
            (
                abs(
                    complex(real[index], imag[index])
                    - complex(real[-1 - index], imag[-1 - index])
                )
                for index in range(pair_count)
            ),
            default=0.0,
        )
        / field_scale
        if field_scale > 0.0
        else math.inf
    )
    strictly_increasing = all(
        positions[index + 1] > positions[index]
        for index in range(max(0, count - 1))
    )
    if count % 2:
        center_bracketed = count > 0 and abs(positions[count // 2]) <= axis_tol * axis_scale
    else:
        center_bracketed = count >= 2 and positions[count // 2 - 1] < 0.0 < positions[count // 2]

    checks = {
        "sample_count_sufficient": count >= int(min_sample_count),
        "axis_strictly_increasing": strictly_increasing,
        "axis_straddles_origin": count > 1 and positions[0] < 0.0 < positions[-1],
        "origin_sampled_or_bracketed": center_bracketed,
        "axis_is_antisymmetric": axis_symmetry_relative <= axis_tol,
        "complex_field_nonzero": field_scale > 0.0,
        "complex_field_is_mirror_symmetric": field_symmetry_relative <= field_tol,
        "solver_residual_converged": residual <= residual_limit,
    }
    return {
        "policy": "symmetric_complex_field_curve_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": count,
            "pair_count": pair_count,
            "axis_min": positions[0] if count else None,
            "axis_max": positions[-1] if count else None,
            "field_scale": field_scale,
            "axis_symmetry_relative": axis_symmetry_relative,
            "field_symmetry_relative": field_symmetry_relative,
            "log10_relative_residual": residual,
        },
        "units": {"axis": str(axis_unit), "field": str(field_unit)},
        "tolerances": {
            "min_sample_count": int(min_sample_count),
            "max_axis_symmetry_relative": axis_tol,
            "max_field_symmetry_relative": field_tol,
            "max_log10_relative_residual": residual_limit,
        },
        "notes": [
            "even sample counts are valid when the origin is bracketed by the central pair",
            "the complex field is compared directly, so magnitude-only phase errors cannot pass",
            "mirror symmetry is a validation identity, not an independent absolute-field reference",
        ],
    }


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
