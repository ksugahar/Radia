"""Solver-neutral force gate for background and material-contrast cases."""

from __future__ import annotations

import math
from typing import Any


def material_contrast_force_gate(
    cases: list[dict[str, Any]],
    *,
    interaction_axis: str = "x",
    max_background_relative_force: float = 0.01,
    max_transverse_relative_force: float = 1.0e-6,
    min_stronger_repulsion_ratio: float = 1.5,
) -> dict[str, Any]:
    """Gate force sign, null background, transverse residual, and contrast order."""

    if not isinstance(cases, list) or len(cases) != 4:
        raise ValueError("cases must contain exactly four role records")
    axis_index = {"x": 0, "y": 1, "z": 2}.get(interaction_axis)
    if axis_index is None:
        raise ValueError("interaction_axis must be x, y, or z")
    relative_tolerances = [
        max_background_relative_force,
        max_transverse_relative_force,
    ]
    if any(
        not math.isfinite(float(value)) or float(value) < 0.0
        for value in relative_tolerances
    ):
        raise ValueError("relative tolerances must be finite and nonnegative")
    if (
        not math.isfinite(float(min_stronger_repulsion_ratio))
        or float(min_stronger_repulsion_ratio) <= 1.0
    ):
        raise ValueError("min_stronger_repulsion_ratio must be finite and greater than one")

    expected_roles = {"background", "attractive", "repulsive_low", "repulsive_high"}
    parsed: dict[str, list[float]] = {}
    metadata_ok = True
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise ValueError(f"case {index} must be an object")
        role = str(case.get("role") or "").strip()
        try:
            force = [float(value) for value in case["force_n"]]
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"case {index}.force_n must be a numeric vector") from exc
        if len(force) != 3 or not all(math.isfinite(value) for value in force):
            raise ValueError(f"case {index}.force_n must contain three finite values")
        if role in parsed:
            raise ValueError(f"duplicate role: {role}")
        parsed[role] = force
        metadata_ok = metadata_ok and case.get("force_unit") == "N" and case.get(
            "coordinate_frame"
        ) == "cartesian"
    if set(parsed) != expected_roles:
        raise ValueError(f"roles must be exactly {sorted(expected_roles)}")

    axial = {role: vector[axis_index] for role, vector in parsed.items()}
    dominant = max(abs(value) for role, value in axial.items() if role != "background")
    background_relative = abs(axial["background"]) / max(dominant, 1.0e-30)
    transverse_relative = {}
    for role, vector in parsed.items():
        transverse = math.sqrt(
            sum(value * value for index, value in enumerate(vector) if index != axis_index)
        )
        transverse_relative[role] = transverse / max(abs(axial[role]), dominant * 1.0e-6)
    stronger_ratio = axial["repulsive_high"] / max(axial["repulsive_low"], 1.0e-30)

    checks = {
        "force_units_and_frame_recorded": metadata_ok,
        "nonbackground_force_scale_is_positive": dominant > 0.0,
        "background_force_is_small": background_relative <= max_background_relative_force,
        "material_contrast_reverses_axial_force": axial["attractive"] < 0.0
        and axial["repulsive_low"] > 0.0
        and axial["repulsive_high"] > 0.0,
        "stronger_repulsive_contrast_increases_force": stronger_ratio
        >= min_stronger_repulsion_ratio,
        "transverse_force_is_bounded": max(transverse_relative.values())
        <= max_transverse_relative_force,
    }
    return {
        "policy": "material_contrast_force_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "interaction_axis": interaction_axis,
        "metrics": {
            "axial_force_by_role_n": axial,
            "dominant_force_n": dominant,
            "background_relative_force": background_relative,
            "transverse_relative_force_by_role": transverse_relative,
            "stronger_repulsion_ratio": stronger_ratio,
        },
        "lesson": (
            "Validate magnetic-material force with a null/background control, an attractive "
            "contrast, and two increasing repulsive contrasts. Gate sign and contrast order "
            "before magnitude regression, and normalize transverse leakage by the dominant "
            "axial force so the near-zero background does not create a false failure."
        ),
    }
