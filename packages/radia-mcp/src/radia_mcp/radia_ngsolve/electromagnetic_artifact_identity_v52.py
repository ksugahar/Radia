"""Magnetic-pressure and frozen-permeability identity checks for v52."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


MAGNETIC_PRESSURE = "magnetic_pressure_fieldjump_normal_traction_owner_identity"
FROZEN_INDUCTANCE = "frozen_permeability_incremental_inductance_bias_perturbation_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _finite_vector(value: object, length: int) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) == length
        and all(_finite(item) for item in value)
    )


def _magnetic_pressure_ok(row: Mapping[str, object]) -> bool:
    field_jump = row.get("field_jump")
    normal = row.get("boundary_normal")
    traction = row.get("traction_n_per_m2")
    field_ok = (
        isinstance(field_jump, Mapping)
        and set(field_jump) == {"normal_h_a_per_m", "tangential_b_t"}
        and all(_finite(value) for value in field_jump.values())
    )
    normal_ok = _finite_vector(normal, 3) and math.isclose(
        sum(float(value) ** 2 for value in normal), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12
    )
    return (
        _generations(
            row,
            "field_generation",
            "normal_generation",
            "traction_generation",
            "owner_generation",
            "result_generation",
        )
        and field_ok
        and row.get("result_field_jump") == field_jump
        and normal_ok
        and row.get("result_boundary_normal") == normal
        and _finite_vector(traction, 3)
        and row.get("result_traction_n_per_m2") == traction
        and row.get("integration_measure") == "surface_area_m2"
        and row.get("result_integration_measure") == row.get("integration_measure")
        and str(row.get("field_owner") or "").startswith("field:")
        and row.get("result_field_owner") == row.get("field_owner")
        and _result(row)
    )


def _frozen_inductance_ok(row: Mapping[str, object]) -> bool:
    bias = row.get("bias_point")
    perturbation = row.get("perturbation")
    bias_current = bias.get("current_a") if isinstance(bias, Mapping) else None
    delta_current = perturbation.get("delta_current_a") if isinstance(perturbation, Mapping) else None
    frequency = perturbation.get("frequency_hz") if isinstance(perturbation, Mapping) else None
    bias_ok = (
        isinstance(bias, Mapping)
        and set(bias) == {"current_a", "solution_sha256"}
        and _finite(bias_current)
        and _digest(bias.get("solution_sha256"))
    )
    perturbation_ok = (
        isinstance(perturbation, Mapping)
        and set(perturbation) == {"delta_current_a", "frequency_hz"}
        and _finite(delta_current)
        and 0.0 < abs(float(delta_current)) <= max(1.0e-6, 0.1 * abs(float(bias_current)))
        and _finite(frequency)
        and float(frequency) > 0.0
    ) if bias_ok else False
    inductance = row.get("incremental_inductance_h")
    return (
        _generations(
            row,
            "permeability_generation",
            "bias_generation",
            "perturbation_generation",
            "inductance_generation",
            "owner_generation",
            "result_generation",
        )
        and row.get("permeability_mode") == "frozen_at_bias"
        and row.get("result_permeability_mode") == row.get("permeability_mode")
        and bias_ok
        and row.get("result_bias_point") == bias
        and perturbation_ok
        and row.get("result_perturbation") == perturbation
        and _finite(inductance)
        and float(inductance) > 0.0
        and row.get("result_incremental_inductance_h") == inductance
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    magnetic_pressure = identity.get(MAGNETIC_PRESSURE)
    frozen_inductance = identity.get(FROZEN_INDUCTANCE)
    if magnetic_pressure is not None:
        checks["v52_magnetic_pressure_field_jump_normal_traction_owner"] = (
            isinstance(magnetic_pressure, Mapping) and _magnetic_pressure_ok(magnetic_pressure)
        )
    if frozen_inductance is not None:
        checks["v52_frozen_permeability_bias_perturbation_inductance_owner"] = (
            isinstance(frozen_inductance, Mapping) and _frozen_inductance_ok(frozen_inductance)
        )
    return checks
