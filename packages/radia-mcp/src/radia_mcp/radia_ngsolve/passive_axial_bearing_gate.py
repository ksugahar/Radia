"""Solver-neutral force and signed-stiffness gate for a passive axial bearing."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any



def _numpy():
    import numpy as np

    return np


def _array(value: object, name: str):
    np = _numpy()
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    result = np.asarray(value, dtype=float).ravel()
    if result.size < 5 or result.size % 2 == 0 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain an odd number of at least five finite values")
    return result


def _finite(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed) or parsed < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return parsed


def passive_axial_bearing_stiffness_gate(
    summary: Mapping[str, object],
) -> dict[str, Any]:
    """Gate symmetry, action-reaction, local stability sign, and replay."""
    np = _numpy()
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    position = _array(summary.get("position_m"), "position_m")
    force0 = _array(summary.get("force_object_0_n"), "force_object_0_n")
    force1 = _array(summary.get("force_object_1_n"), "force_object_1_n")
    if not (position.size == force0.size == force1.size):
        raise ValueError("position and force arrays must have equal length")
    if not np.all(np.diff(position) > 0.0):
        raise ValueError("position_m must be strictly increasing")

    expected = str(summary.get("expected_axial_stability", "")).strip().lower()
    if expected not in {"stable", "unstable"}:
        raise ValueError("expected_axial_stability must be stable or unstable")
    replay_error = _finite(
        summary.get("fresh_replay_relative_error"), "fresh_replay_relative_error"
    )
    center = position.size // 2
    if abs(float(position[center])) > 1.0e-12 * max(float(np.max(np.abs(position))), 1.0):
        raise ValueError("position_m must contain zero at its center")

    scale = max(float(np.max(np.abs(force0))), float(np.max(np.abs(force1))), 1.0e-30)
    axis_scale = max(float(np.max(np.abs(position))), 1.0e-30)
    axis_symmetry = float(np.max(np.abs(position + position[::-1])) / axis_scale)
    force_odd_symmetry = float(np.max(np.abs(force0 + force0[::-1])) / scale)
    action_reaction = float(np.max(np.abs(force0 + force1)) / scale)
    center_force = abs(float(force0[center])) / scale
    signed_stiffness = float(
        (force0[center + 1] - force0[center - 1])
        / (position[center + 1] - position[center - 1])
    )
    observed = "unstable" if signed_stiffness > 0.0 else "stable"
    local_indices = [index for index in range(center - 2, center + 3) if index != center]
    local_products = np.asarray(
        [position[index] * force0[index] for index in local_indices], dtype=float
    )
    local_sign_matches = bool(
        np.all(local_products > 0.0)
        if observed == "unstable"
        else np.all(local_products < 0.0)
    )

    checks = {
        "position_axis_is_symmetric": axis_symmetry <= 1.0e-12,
        "force_is_nearly_odd": force_odd_symmetry <= 1.0e-2,
        "action_reaction_closes": action_reaction <= 2.0e-2,
        "center_force_is_near_zero": center_force <= 5.0e-3,
        "local_force_sign_matches_signed_stiffness": local_sign_matches,
        "observed_stability_matches_expectation": observed == expected,
        "fresh_force_curve_replay_is_stable": replay_error <= 1.0e-9,
    }
    checks = {name: bool(value) for name, value in checks.items()}
    return {
        "policy": "passive_axial_bearing_stiffness_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "point_count": int(position.size),
            "position_start_m": float(position[0]),
            "position_stop_m": float(position[-1]),
            "position_axis_symmetry_relative": axis_symmetry,
            "force_odd_symmetry_relative": force_odd_symmetry,
            "action_reaction_relative": action_reaction,
            "center_force_relative": center_force,
            "center_signed_stiffness_df_dx_n_per_m": signed_stiffness,
            "center_restoring_stiffness_minus_df_dx_n_per_m": -signed_stiffness,
            "observed_axial_stability": observed,
            "fresh_replay_relative_error": replay_error,
        },
        "lesson": (
            "Do not label a passive magnetic bearing restoring from force magnitude alone. "
            "Sweep a signed displacement, check odd symmetry and action-reaction, then "
            "classify local stability from dF/dx. Positive dF/dx means negative restoring "
            "stiffness and axial instability even when the centered force is nearly zero."
        ),
    }
