"""Solver-neutral profile gate for two magnetic-force formulations."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from datetime import datetime


_DIMENSION_UNITS = {
    ("axisymmetric_total", "N"),
    ("3d_total", "N"),
    ("2d_per_length", "N/m"),
}


def _profile(value: object, name: str) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    result = [float(item) for item in value]
    if not result or not all(math.isfinite(item) for item in result):
        raise ValueError(f"{name} must contain finite values")
    return result


def _maximum_relative_difference(left: Sequence[float], right: Sequence[float]) -> float:
    return max(
        abs(a - b) / max(abs(a), abs(b), 1.0e-300)
        for a, b in zip(left, right)
    )


def _trapezoid_integral(positions: Sequence[float], values: Sequence[float]) -> float:
    return sum(
        0.5 * (left_value + right_value) * (right_x - left_x)
        for left_x, right_x, left_value, right_value in zip(
            positions, positions[1:], values, values[1:]
        )
    )


def magnetic_force_method_profile_gate(
    summary: Mapping[str, object],
    *,
    maximum_method_relative_difference: float = 0.05,
    maximum_independent_stress_relative_difference: float = 0.02,
    minimum_selection_scope_relative_difference: float = 0.25,
    maximum_all_body_to_target_magnitude_ratio: float = 0.75,
    maximum_work_relative_difference: float = 0.05,
    maximum_parsed_replay_absolute_difference: float = 1.0e-12,
    minimum_sample_count: int = 5,
) -> dict[str, object]:
    """Compare target-body element force with closed-surface stress force.

    The gate deliberately requires an all-body element-force control. This catches
    a common false comparison where two force formulations are evaluated over
    different bodies or surfaces while being presented as method disagreement.
    """
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    if minimum_sample_count < 3:
        raise ValueError("minimum_sample_count must be at least 3")

    tolerances = {
        "maximum_method_relative_difference": float(maximum_method_relative_difference),
        "maximum_independent_stress_relative_difference": float(
            maximum_independent_stress_relative_difference
        ),
        "minimum_selection_scope_relative_difference": float(
            minimum_selection_scope_relative_difference
        ),
        "maximum_all_body_to_target_magnitude_ratio": float(
            maximum_all_body_to_target_magnitude_ratio
        ),
        "maximum_work_relative_difference": float(maximum_work_relative_difference),
        "maximum_parsed_replay_absolute_difference": float(
            maximum_parsed_replay_absolute_difference
        ),
    }
    if not all(math.isfinite(value) and value >= 0.0 for value in tolerances.values()):
        raise ValueError("all tolerances must be finite and nonnegative")

    profiles = {
        name: _profile(summary.get(name), name)
        for name in (
            "positions",
            "moving_body_element_force",
            "closed_surface_maxwell_stress_force",
            "independent_closed_surface_force",
            "all_body_element_force",
        )
    }
    lengths = {len(values) for values in profiles.values()}
    if len(lengths) != 1:
        raise ValueError("positions and all force profiles must have the same length")

    positions = profiles["positions"]
    target = profiles["moving_body_element_force"]
    stress = profiles["closed_surface_maxwell_stress_force"]
    independent_stress = profiles["independent_closed_surface_force"]
    all_body = profiles["all_body_element_force"]
    increasing_positions = all(right > left for left, right in zip(positions, positions[1:]))

    identity_value = summary.get("artifact_identity")
    identity_present = isinstance(identity_value, Mapping)
    one_sweep_generation_ok = True
    demag_reference_ok = True
    if identity_value is not None and not identity_present:
        one_sweep_generation_ok = False
        demag_reference_ok = False
    elif identity_present:
        generations = identity_value.get("position_force_sample_generations")
        timestamps = identity_value.get("sample_acquired_at_utc")

        def timestamp(value: object) -> float:
            try:
                return datetime.fromisoformat(
                    str(value).replace("Z", "+00:00")
                ).timestamp()
            except (TypeError, ValueError):
                return math.nan

        parsed_times = (
            [timestamp(value) for value in timestamps]
            if isinstance(timestamps, Sequence)
            and not isinstance(timestamps, (str, bytes))
            else []
        )
        one_sweep_generation_ok = (
            isinstance(generations, Sequence)
            and not isinstance(generations, (str, bytes))
            and len(generations) == len(positions)
            and all(bool(str(value)) for value in generations)
            and len({str(value) for value in generations if str(value)}) == 1
            and len(parsed_times) == len(positions)
            and all(math.isfinite(value) for value in parsed_times)
            and all(right > left for left, right in zip(parsed_times, parsed_times[1:]))
        )
        geometry = identity_value.get("magnet_geometry")
        reference = identity_value.get("demag_reference")
        if not isinstance(geometry, Mapping) or not isinstance(reference, Mapping):
            demag_reference_ok = False
        else:
            geometry_revision = str(geometry.get("revision", ""))
            committed_at = timestamp(geometry.get("committed_at_utc"))
            generated_at = timestamp(reference.get("generated_at_utc"))
            demag_reference_ok = (
                bool(geometry_revision)
                and str(reference.get("geometry_revision", "")) == geometry_revision
                and math.isfinite(committed_at)
                and math.isfinite(generated_at)
                and generated_at >= committed_at
            )

    method_difference = _maximum_relative_difference(target, stress)
    independent_stress_difference = _maximum_relative_difference(stress, independent_stress)
    selection_differences = [
        abs(all_value - target_value)
        / max(abs(all_value), abs(target_value), 1.0e-300)
        for all_value, target_value in zip(all_body, target)
    ]
    all_body_ratios = [
        abs(all_value) / max(abs(target_value), 1.0e-300)
        for all_value, target_value in zip(all_body, target)
    ]
    target_integral = _trapezoid_integral(positions, target)
    stress_integral = _trapezoid_integral(positions, stress)
    work_difference = abs(target_integral - stress_integral) / max(
        abs(target_integral), abs(stress_integral), 1.0e-300
    )
    same_sign = all(a * b > 0.0 for a, b in zip(target, stress))
    same_trend = all(
        (right_a - left_a) * (right_b - left_b) >= 0.0
        for left_a, right_a, left_b, right_b in zip(
            target, target[1:], stress, stress[1:]
        )
    )

    replay = summary.get("replay")
    replay = replay if isinstance(replay, Mapping) else {}
    position_unit = str(summary.get("position_unit") or "")
    artifact_position_units = summary.get("artifact_position_units")
    artifact_position_units = (
        artifact_position_units if isinstance(artifact_position_units, Mapping) else None
    )
    artifact_units_match = artifact_position_units is None or (
        bool(artifact_position_units)
        and all(str(unit) == position_unit for unit in artifact_position_units.values())
    )
    parsed_replay = float(replay.get("parsed_max_abs", math.inf))
    checks = {
        "sample_count_sufficient": len(positions) >= minimum_sample_count,
        "positions_strictly_increase": increasing_positions,
        "dimension_and_force_unit_consistent": (
            str(summary.get("quantity_dimension") or ""),
            str(summary.get("force_unit") or ""),
        )
        in _DIMENSION_UNITS,
        "position_unit_recorded": position_unit in {"m", "mm"},
        "artifact_position_units_match_common_grid": artifact_units_match,
        "comparison_axis_recorded": summary.get("comparison_axis") in {"x", "y", "z"},
        "target_methods_share_sign": same_sign,
        "target_methods_share_stepwise_trend": same_trend,
        "target_method_closure_within_tolerance": method_difference
        <= tolerances["maximum_method_relative_difference"],
        "independent_stress_replay_within_tolerance": independent_stress_difference
        <= tolerances["maximum_independent_stress_relative_difference"],
        "selection_scope_is_materially_distinct": min(selection_differences)
        >= tolerances["minimum_selection_scope_relative_difference"],
        "all_body_control_is_not_target_body_force": max(all_body_ratios)
        <= tolerances["maximum_all_body_to_target_magnitude_ratio"],
        "force_position_integrals_close": work_difference
        <= tolerances["maximum_work_relative_difference"],
        "parsed_replay_is_exact_enough": math.isfinite(parsed_replay)
        and 0.0 <= parsed_replay
        <= tolerances["maximum_parsed_replay_absolute_difference"],
        "binary_nonlog_outputs_replay_exact": replay.get(
            "binary_nonlog_outputs_exact"
        )
        is True,
        "position_force_samples_share_one_sweep_generation": one_sweep_generation_ok,
        "demag_reference_matches_current_geometry_revision": demag_reference_ok,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "magnetic_force_method_profile_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "warnings": [] if identity_present else ["artifact_identity_not_recorded"],
        "metrics": {
            "sample_count": len(positions),
            "maximum_method_relative_difference": method_difference,
            "maximum_independent_stress_relative_difference": independent_stress_difference,
            "minimum_selection_scope_relative_difference": min(selection_differences),
            "maximum_all_body_to_target_magnitude_ratio": max(all_body_ratios),
            "target_force_position_integral": target_integral,
            "stress_force_position_integral": stress_integral,
            "force_position_integral_relative_difference": work_difference,
            "parsed_replay_maximum_absolute_difference": parsed_replay,
        },
        "tolerances": tolerances,
        "lesson": (
            "Force-method closure is meaningful only when the target body, closed "
            "stress surface, comparison axis, units, and dimensional convention are "
            "pinned. Bind every artifact's position grid to the same declared unit; equal "
            "numeric coordinates do not prove equal physical positions. Keep an all-body "
            "force as a negative control so a selection-scope "
            "error cannot masquerade as disagreement between force formulations."
        ),
    }
