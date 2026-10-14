"""Solver-independent force reversal gate for a facing permanent-magnet pair."""
from __future__ import annotations

import json
import math


def permanent_magnet_force_pair_gate(
    summary_json: str,
    magnitude_relative_tolerance: float = 2.0e-2,
    off_axis_relative_tolerance: float = 1.0e-3,
) -> dict:
    """Check attraction/repulsion reversal without depending on a solver format."""

    if magnitude_relative_tolerance <= 0.0 or off_axis_relative_tolerance <= 0.0:
        raise ValueError("tolerances must be positive")
    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    if not isinstance(summary, dict):
        raise ValueError("summary_json must decode to an object")

    rows = summary.get("cases")
    if not isinstance(rows, list) or len(rows) != 2:
        raise ValueError("cases must contain exactly two records")
    indexed = {row.get("pole_relation"): row for row in rows if isinstance(row, dict)}
    like = indexed.get("like")
    opposite = indexed.get("opposite")
    if like is None or opposite is None:
        raise ValueError("cases must contain pole_relation 'like' and 'opposite'")

    def vector(row: dict, name: str) -> list[float]:
        value = row.get(name)
        if not isinstance(value, list) or len(value) != 3:
            raise ValueError(f"{name} must be a three-component list")
        out = [float(component) for component in value]
        if not all(math.isfinite(component) for component in out):
            raise ValueError(f"{name} components must be finite")
        return out

    axis = summary.get("interaction_axis")
    if axis not in {"x", "y", "z"}:
        raise ValueError("interaction_axis must be x, y, or z")
    reference_length_m = float(summary.get("reference_length_m", 0.0))
    if not math.isfinite(reference_length_m) or reference_length_m <= 0.0:
        raise ValueError("reference_length_m must be positive and finite")
    axis_index = {"x": 0, "y": 1, "z": 2}[axis]
    force_like = vector(like, "force_N")
    force_opposite = vector(opposite, "force_N")
    torque_like = vector(like, "torque_Nm")
    torque_opposite = vector(opposite, "torque_Nm")
    axial_like = force_like[axis_index]
    axial_opposite = force_opposite[axis_index]
    scale = max(abs(axial_like), abs(axial_opposite), 1.0e-300)
    magnitude_mismatch = abs(abs(axial_like) - abs(axial_opposite)) / scale
    off_axis = max(
        *(abs(value) for i, value in enumerate(force_like) if i != axis_index),
        *(abs(value) for i, value in enumerate(force_opposite) if i != axis_index),
    ) / scale
    torque_residual = max(*(abs(value) for value in torque_like), *(abs(value) for value in torque_opposite)) / (scale * reference_length_m)

    checks = {
        "force_unit_is_N": summary.get("force_unit") == "N",
        "torque_unit_is_Nm": summary.get("torque_unit") == "N*m",
        "component_frame_is_global_cartesian": summary.get("component_frame") == "global_cartesian",
        "sign_convention_explicit": summary.get("positive_axis_interaction") == "repulsion",
        "like_poles_repel": axial_like > 0.0,
        "opposite_poles_attract": axial_opposite < 0.0,
        "axial_magnitudes_match": magnitude_mismatch <= magnitude_relative_tolerance,
        "off_axis_force_small": off_axis <= off_axis_relative_tolerance,
        "torque_residual_small": torque_residual <= off_axis_relative_tolerance,
    }
    return {
        "policy": "permanent_magnet_force_pair_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "axial_magnitude_relative_mismatch": magnitude_mismatch,
            "off_axis_force_relative_max": off_axis,
            "torque_relative_max": torque_residual,
        },
        "notes": [
            "like/opposite pole cases must reuse geometry, selection, and force convention",
            "the gate checks force reversal and symmetry; it does not prescribe a force formulation",
        ],
    }
