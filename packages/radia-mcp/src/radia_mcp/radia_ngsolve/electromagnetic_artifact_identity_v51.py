"""Incremental-material and weighted-force artifact identity checks for v51."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


INCREMENTAL = "incremental_frozen_bias_harmonic_tangent_branch_owner_identity"
WEIGHTED_FORCE = "weighted_stress_mask_air_axisym_factor_force_frame_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_positive(value: object) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and float(value) > 0.0
    )


def _finite_vector(value: object, length: int) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) == length
        and all(
            isinstance(item, (int, float))
            and not isinstance(item, bool)
            and math.isfinite(float(item))
            for item in value
        )
    )


def _incremental_ok(row: Mapping[str, object]) -> bool:
    tangent = row.get("tangent_permeability_relative")
    frequency = row.get("harmonic_frequency_hz")
    tangent_ok = (
        isinstance(tangent, Mapping)
        and set(tangent) == {"radial", "tangential"}
        and all(_finite_positive(value) for value in tangent.values())
    )
    return (
        _generations(
            row,
            "bias_generation",
            "harmonic_generation",
            "tangent_generation",
            "branch_generation",
            "owner_generation",
            "result_generation",
        )
        and row.get("analysis_mode") == "incremental_permeability"
        and row.get("result_analysis_mode") == row.get("analysis_mode")
        and _digest(row.get("frozen_bias_solution_sha256"))
        and row.get("result_frozen_bias_solution_sha256") == row.get("frozen_bias_solution_sha256")
        and _finite_positive(frequency)
        and row.get("result_harmonic_frequency_hz") == frequency
        and tangent_ok
        and row.get("result_tangent_permeability_relative") == tangent
        and row.get("branch_state") in {"ascending_major_loop", "descending_major_loop", "recoil_branch"}
        and row.get("result_branch_state") == row.get("branch_state")
        and str(row.get("operating_point_owner") or "").startswith("operating-point:")
        and row.get("result_operating_point_owner") == row.get("operating_point_owner")
        and _result(row)
    )


def _weighted_force_ok(row: Mapping[str, object]) -> bool:
    air_elements = row.get("air_element_ids")
    radius = row.get("axisymmetric_radius_m")
    factor = row.get("axisymmetric_factor_m")
    force = row.get("force_n")
    return (
        _generations(
            row,
            "mask_generation",
            "air_generation",
            "axisym_generation",
            "force_generation",
            "frame_generation",
            "owner_generation",
            "result_generation",
        )
        and _digest(row.get("weighted_stress_mask_sha256"))
        and row.get("result_weighted_stress_mask_sha256") == row.get("weighted_stress_mask_sha256")
        and isinstance(air_elements, list)
        and bool(air_elements)
        and all(isinstance(item, int) and not isinstance(item, bool) and item > 0 for item in air_elements)
        and air_elements == sorted(set(air_elements))
        and row.get("result_air_element_ids") == air_elements
        and _finite_positive(radius)
        and row.get("result_axisymmetric_radius_m") == radius
        and _finite_positive(factor)
        and math.isclose(float(factor), 2.0 * math.pi * float(radius), rel_tol=1e-12, abs_tol=1e-15)
        and row.get("result_axisymmetric_factor_m") == factor
        and _finite_vector(force, 2)
        and row.get("result_force_n") == force
        and row.get("force_frame") == "global_rz"
        and row.get("result_force_frame") == row.get("force_frame")
        and str(row.get("force_owner") or "").startswith("force:")
        and row.get("result_force_owner") == row.get("force_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    incremental = identity.get(INCREMENTAL)
    weighted_force = identity.get(WEIGHTED_FORCE)
    if incremental is not None:
        checks["v51_incremental_bias_harmonic_tangent_branch_owner"] = (
            isinstance(incremental, Mapping) and _incremental_ok(incremental)
        )
    if weighted_force is not None:
        checks["v51_weighted_stress_air_axisym_factor_force_frame_owner"] = (
            isinstance(weighted_force, Mapping) and _weighted_force_ok(weighted_force)
        )
    return checks
