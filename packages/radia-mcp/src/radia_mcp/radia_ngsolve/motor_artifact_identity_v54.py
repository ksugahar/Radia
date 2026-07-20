"""Torque-angle and irreversible-demagnetization artifact checks for v54."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .motor_artifact_identity_v55 import validate_public_identity as validate_public_v55_identity


TORQUE = "torqueripple_harmonic_mechanical_electrical_angle_polepair_owner_identity"
DEMAG = "demag_irreversible_knee_temperature_currentvector_recovery_owner_identity"


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


def _numeric_sequence(value: object, length: int | None = None) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and bool(value) and (length is None or len(value) == length) and all(_finite(item) for item in value)


def _torque_ok(row: Mapping[str, object]) -> bool:
    pole_pairs = row.get("pole_pairs")
    mechanical = row.get("mechanical_angles_deg")
    electrical = row.get("electrical_angles_deg")
    angles_ok = (
        isinstance(pole_pairs, int)
        and not isinstance(pole_pairs, bool)
        and pole_pairs > 0
        and _numeric_sequence(mechanical)
        and _numeric_sequence(electrical, len(mechanical))
        and all(math.isclose(float(electric) % 360.0, (float(mechanic) * pole_pairs) % 360.0, rel_tol=0.0, abs_tol=1.0e-10) for mechanic, electric in zip(mechanical, electrical))
    )
    harmonics = row.get("torque_harmonics")
    harmonics_ok = isinstance(harmonics, Sequence) and not isinstance(harmonics, (str, bytes)) and bool(harmonics)
    orders: set[int] = set()
    if harmonics_ok:
        for harmonic in harmonics:
            if not isinstance(harmonic, Mapping) or set(harmonic) != {"order", "amplitude_nm", "phase_electrical_deg"}:
                harmonics_ok = False
                break
            order = harmonic["order"]
            if not (
                isinstance(order, int) and not isinstance(order, bool) and order > 0 and order not in orders
                and _finite(harmonic["amplitude_nm"]) and float(harmonic["amplitude_nm"]) >= 0.0
                and _finite(harmonic["phase_electrical_deg"]) and -180.0 <= float(harmonic["phase_electrical_deg"]) <= 180.0
            ):
                harmonics_ok = False
                break
            orders.add(order)
    return (
        _generations(row, "harmonic_generation", "mechanical_generation", "electrical_generation", "polepair_generation", "owner_generation", "result_generation")
        and angles_ok
        and row.get("result_pole_pairs") == pole_pairs
        and row.get("result_mechanical_angles_deg") == mechanical
        and row.get("result_electrical_angles_deg") == electrical
        and harmonics_ok
        and row.get("result_torque_harmonics") == harmonics
        and str(row.get("result_owner") or "").startswith("result:")
        and row.get("accepted_result_owner") == row.get("result_owner")
        and _result(row)
    )


def _demag_ok(row: Mapping[str, object]) -> bool:
    knee = row.get("knee_criterion")
    knee_ok = (
        isinstance(knee, Mapping)
        and set(knee) == {"b_t", "h_a_per_m", "criterion"}
        and _finite(knee["b_t"]) and float(knee["b_t"]) >= 0.0
        and _finite(knee["h_a_per_m"]) and float(knee["h_a_per_m"]) < 0.0
        and knee["criterion"] == "operating_point_below_knee"
    )
    current = row.get("current_vector_abc_a")
    current_ok = _numeric_sequence(current, 3) and math.isclose(sum(float(value) for value in current), 0.0, rel_tol=0.0, abs_tol=1.0e-10)
    fraction = row.get("irreversible_demag_fraction")
    recovery = row.get("post_recovery_remanence_fraction")
    fraction_ok = _finite(fraction) and 0.0 <= float(fraction) <= 1.0 and _finite(recovery) and math.isclose(float(recovery), 1.0 - float(fraction), rel_tol=0.0, abs_tol=1.0e-12)
    expected_state = "reversible" if fraction_ok and float(fraction) == 0.0 else "fully_demagnetized" if fraction_ok and float(fraction) == 1.0 else "partially_demagnetized" if fraction_ok else None
    return (
        _generations(row, "knee_generation", "temperature_generation", "current_generation", "recovery_generation", "owner_generation", "result_generation")
        and knee_ok
        and row.get("result_knee_criterion") == knee
        and _finite(row.get("temperature_c")) and float(row["temperature_c"]) > -273.15
        and row.get("result_temperature_c") == row.get("temperature_c")
        and current_ok
        and row.get("result_current_vector_abc_a") == current
        and fraction_ok
        and row.get("result_irreversible_demag_fraction") == fraction
        and row.get("result_post_recovery_remanence_fraction") == recovery
        and row.get("recovery_state") == expected_state
        and row.get("result_recovery_state") == expected_state
        and str(row.get("magnet_owner") or "").startswith("magnet:")
        and row.get("result_magnet_owner") == row.get("magnet_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks = validate_public_v55_identity(identity)
    torque = identity.get(TORQUE)
    demag = identity.get(DEMAG)
    if torque is not None:
        checks["motor_v54_torque_harmonic_mechanical_electrical_polepair_owner"] = isinstance(torque, Mapping) and _torque_ok(torque)
    if demag is not None:
        checks["motor_v54_demag_knee_temperature_current_recovery_owner"] = isinstance(demag, Mapping) and _demag_ok(demag)
    return checks
