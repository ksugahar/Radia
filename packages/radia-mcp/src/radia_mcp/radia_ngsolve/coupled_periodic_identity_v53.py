"""Electrothermal-contact and Floquet-periodic identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .conservation_identity_v54 import validate_public_v54_identity


ELECTROTHERMAL = "electrothermal_contact_resistance_power_heat_reciprocity_time_owner_identity"
FLOQUET = "floquet_phase_wavevector_boundarypair_orientation_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _number(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(float(value)):
        return False
    return (not positive or float(value) > 0.0) and (not nonnegative or float(value) >= 0.0)


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-12)


def _result_identity(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _electrothermal_ok(row: Mapping[str, object]) -> bool:
    resistance = row.get("contact_resistance_ohm")
    current = row.get("contact_current_a")
    electric_power = row.get("electric_power_w")
    deposited_heat = row.get("deposited_heat_w")
    time_s = row.get("time_s")
    return (
        _generation(row, "contact_generation", "electric_generation", "thermal_generation", "time_generation", "owner_generation", "result_generation")
        and _number(resistance, positive=True)
        and _number(current)
        and _number(electric_power, nonnegative=True)
        and _number(deposited_heat, nonnegative=True)
        and _number(time_s, nonnegative=True)
        and _close(electric_power, float(current) ** 2 * float(resistance))
        and _close(deposited_heat, electric_power)
        and all(row.get("result_" + name) == row.get(name) for name in ("contact_resistance_ohm", "contact_current_a", "electric_power_w", "deposited_heat_w", "time_s"))
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity(row)
    )


def _vector(value: object) -> list[float] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        return None
    if not all(_number(item) for item in value):
        return None
    return [float(item) for item in value]


def _floquet_ok(row: Mapping[str, object]) -> bool:
    wave_vector = _vector(row.get("wave_vector_per_m"))
    translation = _vector(row.get("translation_m"))
    phase = row.get("phase_rad")
    pair = row.get("boundary_pair")
    if wave_vector is None or translation is None or not _number(phase) or not isinstance(pair, Mapping):
        return False
    expected_phase = sum(left * right for left, right in zip(wave_vector, translation))
    wrapped_error = math.atan2(math.sin(float(phase) - expected_phase), math.cos(float(phase) - expected_phase))
    pair_ok = set(pair) == {"source", "destination"} and all(isinstance(pair[name], str) and pair[name] for name in pair) and pair["source"] != pair["destination"]
    return (
        _generation(row, "phase_generation", "wavevector_generation", "pair_generation", "orientation_generation", "owner_generation", "result_generation")
        and math.isclose(wrapped_error, 0.0, abs_tol=1.0e-12)
        and pair_ok
        and row.get("orientation") == "source_to_destination"
        and all(row.get("result_" + name) == row.get(name) for name in ("phase_rad", "wave_vector_per_m", "translation_m", "boundary_pair", "orientation"))
        and str(row.get("field_owner") or "").startswith("field:")
        and row.get("result_field_owner") == row.get("field_owner")
        and _result_identity(row)
    )


def validate_public_v53_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    v54 = validate_public_v54_identity(payload)
    if v54:
        checks.update(v54["checks"])
    electrothermal = payload.get(ELECTROTHERMAL)
    floquet = payload.get(FLOQUET)
    if electrothermal is not None:
        checks["v53_electrothermal_contact_power_heat_time_owner"] = isinstance(electrothermal, Mapping) and _electrothermal_ok(electrothermal)
    if floquet is not None:
        checks["v53_floquet_phase_wavevector_pair_orientation_owner"] = isinstance(floquet, Mapping) and _floquet_ok(floquet)
    if not checks:
        return {}
    return {"policy": "coupled_periodic_identity_v53", "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, accepted in checks.items() if not accepted]}
