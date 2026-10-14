"""Solver-independent force-vector rotation covariance gate."""
from __future__ import annotations

import json
import math


def evaluate_force_rotation_covariance(
    reference_force: dict,
    rotated_force: dict,
    rotation_deg: float,
    relative_tolerance: float = 1.0e-3,
) -> dict:
    """Compare a 2D force vector with the expected rigid rotation."""
    if relative_tolerance <= 0:
        raise ValueError("relative_tolerance must be positive")
    try:
        fx0 = float(reference_force["Fx"])
        fy0 = float(reference_force["Fy"])
        fx1 = float(rotated_force["Fx"])
        fy1 = float(rotated_force["Fy"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("force objects must contain numeric Fx and Fy") from exc
    angle = math.radians(float(rotation_deg))
    expected = {
        "Fx": math.cos(angle) * fx0 - math.sin(angle) * fy0,
        "Fy": math.sin(angle) * fx0 + math.cos(angle) * fy0,
    }
    reference_magnitude = math.hypot(fx0, fy0)
    rotated_magnitude = math.hypot(fx1, fy1)
    scale = max(reference_magnitude, rotated_magnitude, 1.0e-30)
    vector_relative_error = math.hypot(fx1 - expected["Fx"], fy1 - expected["Fy"]) / scale
    magnitude_relative_error = abs(rotated_magnitude - reference_magnitude) / scale
    checks = {
        "reference_force_nonzero": reference_magnitude > 0.0,
        "vector_rotates_covariantly": vector_relative_error <= relative_tolerance,
        "force_magnitude_preserved": magnitude_relative_error <= relative_tolerance,
    }
    return {
        "schema": "radia-motor-force-rotation-covariance/v1",
        "policy": "force_vector_must_follow_geometry_excitation_rotation",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "rotation_deg": float(rotation_deg),
        "reference_force": {"Fx": fx0, "Fy": fy0},
        "rotated_force": {"Fx": fx1, "Fy": fy1},
        "expected_rotated_force": expected,
        "reference_magnitude": reference_magnitude,
        "rotated_magnitude": rotated_magnitude,
        "vector_relative_error": vector_relative_error,
        "magnitude_relative_error": magnitude_relative_error,
        "relative_tolerance": relative_tolerance,
        "checks": checks,
    }


def force_rotation_covariance_gate(
    reference_force_json: str,
    rotated_force_json: str,
    rotation_deg: float,
    relative_tolerance: float = 1.0e-3,
) -> str:
    try:
        reference = json.loads(reference_force_json)
        rotated = json.loads(rotated_force_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"force inputs must be valid JSON: {exc.msg}") from exc
    if not isinstance(reference, dict) or not isinstance(rotated, dict):
        raise ValueError("force inputs must decode to objects")
    return json.dumps(
        evaluate_force_rotation_covariance(
            reference, rotated, rotation_deg, relative_tolerance
        ),
        indent=2,
        sort_keys=True,
    )
