"""Solver-neutral analytic gate for a finite solenoid surface-current sheet."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _number(value: object, name: str, *, positive: bool = False) -> float:
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
    parsed = [_number(item, f"{name}[{index}]") for index, item in enumerate(value)]
    if expected is not None and len(parsed) != expected:
        raise ValueError(f"{name} must contain {expected} values")
    return parsed


def finite_solenoid_surface_current_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Select a mesh by analytic axial-field error and gate signed linearity."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    geometry = summary.get("geometry")
    units = summary.get("units")
    tolerances = summary.get("gate_tolerances")
    linearity = summary.get("linearity")
    levels = summary.get("levels")
    if not all(isinstance(item, Mapping) for item in (geometry, units, tolerances, linearity)):
        raise ValueError("geometry, units, gate_tolerances, and linearity must be objects")
    if not isinstance(levels, list) or len(levels) < 2:
        raise ValueError("levels must contain at least two mesh records")

    radius = _number(geometry.get("radius"), "geometry.radius", positive=True)
    length = _number(geometry.get("length"), "geometry.length", positive=True)
    center_z = _number(geometry.get("center_z"), "geometry.center_z")
    current_density = _number(summary.get("current_density_a_per_m"), "current_density_a_per_m")
    if current_density == 0.0:
        raise ValueError("current_density_a_per_m must be nonzero")
    z = _array(summary.get("profile_z"), "profile_z")
    minimum_samples = int(_number(tolerances.get("minimum_samples"), "minimum_samples", positive=True))
    if len(z) < minimum_samples or len(z) % 2 != 1:
        raise ValueError("profile_z must contain an odd number of sufficiently many samples")
    if any(right <= left for left, right in zip(z, z[1:])):
        raise ValueError("profile_z must be strictly increasing")

    mu0 = 4.0 * math.pi * 1.0e-7
    relative_z = [value - center_z for value in z]
    analytic = [
        mu0 * current_density / 2.0 * (
            (value + length / 2.0) / math.sqrt(radius**2 + (value + length / 2.0) ** 2)
            - (value - length / 2.0) / math.sqrt(radius**2 + (value - length / 2.0) ** 2)
        )
        for value in relative_z
    ]
    analytic_center = mu0 * current_density * length / (
        2.0 * math.sqrt(radius**2 + (length / 2.0) ** 2)
    )
    scale = abs(analytic_center)
    center_index = min(range(len(z)), key=lambda index: abs(z[index] - center_z))

    parsed_levels: list[dict[str, Any]] = []
    labels: set[str] = set()
    for index, level in enumerate(levels):
        if not isinstance(level, Mapping):
            raise ValueError(f"levels[{index}] must be an object")
        label = str(level.get("label") or "").strip()
        if not label or label in labels:
            raise ValueError("mesh level labels must be nonempty and unique")
        labels.add(label)
        bx = _array(level.get("Bx_T"), f"levels[{index}].Bx_T", len(z))
        by = _array(level.get("By_T"), f"levels[{index}].By_T", len(z))
        bz = _array(level.get("Bz_T"), f"levels[{index}].Bz_T", len(z))
        tet_count = int(_number(level.get("tet_count"), f"levels[{index}].tet_count", positive=True))
        minimum_quality = _number(level.get("minimum_mesh_quality"), f"levels[{index}].minimum_mesh_quality", positive=True)
        errors = [(value - reference) / scale for value, reference in zip(bz, analytic)]
        parsed_levels.append({
            "label": label,
            "tet_count": tet_count,
            "minimum_mesh_quality": minimum_quality,
            "center_relative_error": abs(bz[center_index] - analytic_center) / scale,
            "maximum_error_normalized_by_center": max(abs(value) for value in errors),
            "rms_error_normalized_by_center": math.sqrt(sum(value * value for value in errors) / len(errors)),
            "mirror_symmetry_relative": max(abs(left - right) for left, right in zip(bz, reversed(bz))) / scale,
            "transverse_leakage_relative": max(math.hypot(x, y) for x, y in zip(bx, by)) / scale,
        })
    best = min(parsed_levels, key=lambda item: item["rms_error_normalized_by_center"])
    baseline = next((item for item in parsed_levels if item["label"] == "source_auto"), None)
    if baseline is None:
        raise ValueError("levels must include source_auto")
    improvement = baseline["maximum_error_normalized_by_center"] / max(best["maximum_error_normalized_by_center"], 1.0e-300)

    k_values = _array(linearity.get("current_density_a_per_m"), "linearity.current_density_a_per_m", 3)
    bz_rows_raw = linearity.get("Bz_T")
    if not isinstance(bz_rows_raw, Sequence) or len(bz_rows_raw) != 3:
        raise ValueError("linearity.Bz_T must contain three profiles")
    bz_rows = [_array(row, f"linearity.Bz_T[{index}]", len(z)) for index, row in enumerate(bz_rows_raw)]
    expected_k = [current_density, 2.0 * current_density, -current_density]
    k_contract = all(math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-15) for actual, expected in zip(k_values, expected_k))
    k2_error = max(abs(doubled - 2.0 * base) for base, doubled in zip(bz_rows[0], bz_rows[1])) / (2.0 * scale)
    kminus_error = max(abs(negative + base) for base, negative in zip(bz_rows[0], bz_rows[2])) / scale

    limits = {
        "center": _number(tolerances.get("maximum_center_relative_error"), "maximum_center_relative_error", positive=True),
        "maximum": _number(tolerances.get("maximum_profile_error_normalized_by_center"), "maximum_profile_error_normalized_by_center", positive=True),
        "rms": _number(tolerances.get("maximum_rms_error_normalized_by_center"), "maximum_rms_error_normalized_by_center", positive=True),
        "symmetry": _number(tolerances.get("maximum_mirror_symmetry_relative"), "maximum_mirror_symmetry_relative", positive=True),
        "transverse": _number(tolerances.get("maximum_transverse_leakage_relative"), "maximum_transverse_leakage_relative", positive=True),
        "linearity": _number(tolerances.get("maximum_linearity_relative_error"), "maximum_linearity_relative_error", positive=True),
        "improvement": _number(tolerances.get("minimum_baseline_to_best_improvement"), "minimum_baseline_to_best_improvement", positive=True),
    }
    axis_symmetry = max(abs((left - center_z) + (right - center_z)) for left, right in zip(z, reversed(z))) / length
    checks = {
        "units_explicit": units.get("length") in {"m", "cm", "mm"} and units.get("field") == "T" and units.get("surface_current_density") == "A/m",
        "axis_is_centered_and_symmetric": axis_symmetry <= 1.0e-12,
        "mesh_metadata_is_physical": all(row["tet_count"] > 0 and 0.0 < row["minimum_mesh_quality"] <= 1.0 for row in parsed_levels),
        "best_center_matches_analytic_field": best["center_relative_error"] <= limits["center"],
        "best_profile_max_error_is_bounded": best["maximum_error_normalized_by_center"] <= limits["maximum"],
        "best_profile_rms_error_is_bounded": best["rms_error_normalized_by_center"] <= limits["rms"],
        "best_profile_is_mirror_symmetric": best["mirror_symmetry_relative"] <= limits["symmetry"],
        "best_profile_transverse_leakage_is_bounded": best["transverse_leakage_relative"] <= limits["transverse"],
        "analytic_profile_improves_over_source_mesh": improvement >= limits["improvement"],
        "signed_current_density_contract_recorded": k_contract,
        "positive_and_negative_scaling_is_linear": max(k2_error, kminus_error) <= limits["linearity"],
        "center_field_sign_tracks_current": bz_rows[0][center_index] * current_density > 0.0 and bz_rows[1][center_index] * current_density > 0.0 and bz_rows[2][center_index] * current_density < 0.0,
    }
    issues = [name for name, ok in checks.items() if not ok]
    largest = max(parsed_levels, key=lambda item: item["tet_count"])
    return {
        "policy": "finite_solenoid_surface_current_profile_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "analytic_center_field_T": analytic_center,
            "best_level": best,
            "baseline_level": baseline,
            "baseline_to_best_max_error_improvement": improvement,
            "largest_tet_level_label": largest["label"],
            "best_is_largest_tet_mesh": best["label"] == largest["label"],
            "k2_linearity_relative_error": k2_error,
            "kminus1_linearity_relative_error": kminus_error,
        },
        "lesson": "Validate an azimuthal surface-current sheet against the finite-solenoid axis formula, select the mesh by field-profile error rather than element count, and verify positive scaling and sign reversal.",
    }
