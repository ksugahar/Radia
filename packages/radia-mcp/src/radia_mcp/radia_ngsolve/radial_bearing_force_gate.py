"""Solver-neutral excitation-symmetry gate for radial magnetic-bearing force."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


_ROLES = ("balanced", "positive_y", "negative_y")
_EXCITATION_ORDER = ("positive_x", "positive_y", "negative_x", "negative_y")


def _finite(value: object, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _vector(value: object, name: str, size: int) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    result = [_finite(item, name) for item in value]
    if len(result) != size:
        raise ValueError(f"{name} must contain exactly {size} values")
    return result


def radial_bearing_force_symmetry_gate(
    summary: Mapping[str, object],
    *,
    max_null_relative_force: float = 1.0e-3,
    max_transverse_relative_force: float = 1.0e-3,
    max_mirror_relative_error: float = 1.0e-3,
    max_replay_relative_error: float = 1.0e-9,
) -> dict[str, Any]:
    """Gate null excitation, mirrored force, transverse leakage, and replay.

    The three excitation cases must use the order ``+x,+y,-x,-y``.  The
    positive and negative y cases preserve the x-axis bias and exchange the
    two y-axis currents.  This makes force covariance a stronger check than a
    single magnitude regression and applies to any magnetic-body force route.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    limits = {
        "max_null_relative_force": max_null_relative_force,
        "max_transverse_relative_force": max_transverse_relative_force,
        "max_mirror_relative_error": max_mirror_relative_error,
        "max_replay_relative_error": max_replay_relative_error,
    }
    parsed_limits = {name: _finite(value, name) for name, value in limits.items()}
    if any(value < 0.0 for value in parsed_limits.values()):
        raise ValueError("relative tolerances must be nonnegative")

    order = tuple(str(item) for item in summary.get("excitation_order", []))
    cases = summary.get("cases")
    if not isinstance(cases, Sequence) or isinstance(cases, (str, bytes)):
        raise ValueError("cases must be an array")
    by_role: dict[str, dict[str, object]] = {}
    metadata_ok = True
    for index, case in enumerate(cases):
        if not isinstance(case, Mapping):
            raise ValueError(f"cases[{index}] must be an object")
        role = str(case.get("role", ""))
        if role in by_role:
            raise ValueError(f"duplicate role: {role}")
        force = _vector(case.get("force_n"), f"cases[{index}].force_n", 3)
        excitation = _vector(
            case.get("excitation_a"), f"cases[{index}].excitation_a", 4
        )
        metadata_ok = metadata_ok and case.get("force_unit") == "N"
        metadata_ok = metadata_ok and case.get("coordinate_frame") == "cartesian"
        by_role[role] = {"force": force, "excitation": excitation}
    if tuple(by_role) != _ROLES:
        raise ValueError(f"cases must be ordered exactly as {list(_ROLES)}")

    balanced = by_role["balanced"]
    positive = by_role["positive_y"]
    negative = by_role["negative_y"]
    balanced_excitation = balanced["excitation"]
    positive_excitation = positive["excitation"]
    negative_excitation = negative["excitation"]
    assert isinstance(balanced_excitation, list)
    assert isinstance(positive_excitation, list)
    assert isinstance(negative_excitation, list)
    excitation_scale = max(
        *(abs(value) for value in positive_excitation + negative_excitation), 1.0e-30
    )
    balanced_excitation_error = max(
        abs(value - balanced_excitation[0]) for value in balanced_excitation
    ) / excitation_scale
    x_bias_error = max(
        abs(positive_excitation[0] - positive_excitation[2]),
        abs(negative_excitation[0] - negative_excitation[2]),
        abs(positive_excitation[0] - negative_excitation[0]),
    ) / excitation_scale
    y_mirror_error = max(
        abs(positive_excitation[1] - negative_excitation[3]),
        abs(positive_excitation[3] - negative_excitation[1]),
    ) / excitation_scale
    y_bias_error = max(
        abs(
            positive_excitation[1]
            + positive_excitation[3]
            - 2.0 * positive_excitation[0]
        ),
        abs(
            negative_excitation[1]
            + negative_excitation[3]
            - 2.0 * negative_excitation[0]
        ),
    ) / excitation_scale

    balanced_force = balanced["force"]
    positive_force = positive["force"]
    negative_force = negative["force"]
    assert isinstance(balanced_force, list)
    assert isinstance(positive_force, list)
    assert isinstance(negative_force, list)
    dominant = max(abs(positive_force[1]), abs(negative_force[1]), 1.0e-30)
    null_relative = math.hypot(*balanced_force) / dominant
    transverse_relative = max(
        math.hypot(positive_force[0], positive_force[2]),
        math.hypot(negative_force[0], negative_force[2]),
    ) / dominant
    mirror_relative = abs(positive_force[1] + negative_force[1]) / dominant
    replay_error = _finite(
        summary.get("fresh_replay_relative_error"), "fresh_replay_relative_error"
    )
    method = str(summary.get("force_method", ""))

    checks = {
        "force_units_and_frame_recorded": metadata_ok,
        "excitation_order_recorded": order == _EXCITATION_ORDER,
        "balanced_excitation_is_equal": balanced_excitation_error <= 1.0e-12,
        "x_axis_bias_is_preserved": x_bias_error <= 1.0e-12,
        "y_axis_excitations_are_mirrored": y_mirror_error <= 1.0e-12,
        "y_axis_bias_is_preserved": y_bias_error <= 1.0e-12,
        "magnetic_body_force_method_recorded": method
        in {"weighted_stress_volume_integral", "coenergy_virtual_work"},
        "mirrored_axial_forces_reverse_sign": positive_force[1] > 0.0
        and negative_force[1] < 0.0,
        "balanced_force_is_near_zero": null_relative
        <= parsed_limits["max_null_relative_force"],
        "transverse_force_is_bounded": transverse_relative
        <= parsed_limits["max_transverse_relative_force"],
        "mirrored_axial_magnitudes_match": mirror_relative
        <= parsed_limits["max_mirror_relative_error"],
        "fresh_replay_is_stable": replay_error
        <= parsed_limits["max_replay_relative_error"],
    }
    return {
        "policy": "radial_bearing_force_symmetry_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "dominant_axial_force_n": dominant,
            "balanced_relative_force": null_relative,
            "maximum_transverse_relative_force": transverse_relative,
            "mirrored_axial_relative_error": mirror_relative,
            "balanced_excitation_relative_error": balanced_excitation_error,
            "x_bias_relative_error": x_bias_error,
            "y_mirror_relative_error": y_mirror_error,
            "y_bias_relative_error": y_bias_error,
            "fresh_replay_relative_error": replay_error,
        },
        "lesson": (
            "For magnetic-body force, pair the weighted-stress or coenergy result with "
            "an equal-current null control and a mirrored excitation. Gate sign reversal, "
            "magnitude covariance, and transverse leakage before accepting magnitude."
        ),
    }
