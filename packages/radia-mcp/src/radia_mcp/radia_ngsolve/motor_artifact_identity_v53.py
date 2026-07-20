"""Skewed-machine torque and magnet demagnetization identity checks for v53."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .motor_artifact_identity_v54 import validate_public_identity as validate_public_v54_identity


SKEW = "skew_slice_weight_angle_harmonic_torque_rotor_owner_identity"
DEMAG = "magnet_demag_operatingpoint_temperature_recoil_irreversible_owner_identity"


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


def _numeric_list(value: object) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and bool(value) and all(_finite(item) for item in value)


def _skew_ok(row: Mapping[str, object]) -> bool:
    weights = row.get("slice_weights")
    angles = row.get("skew_angles_mechanical_deg")
    harmonics = row.get("harmonic_torque")
    slices_ok = (
        _numeric_list(weights)
        and _numeric_list(angles)
        and len(weights) == len(angles)
        and len(weights) >= 2
        and all(float(weight) > 0.0 for weight in weights)
        and math.isclose(sum(float(weight) for weight in weights), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(float(left) < float(right) for left, right in zip(angles, angles[1:]))
        and math.isclose(sum(float(weight) * float(angle) for weight, angle in zip(weights, angles)), 0.0, rel_tol=0.0, abs_tol=1.0e-12)
    )
    harmonics_ok = isinstance(harmonics, Sequence) and not isinstance(harmonics, (str, bytes)) and bool(harmonics)
    if harmonics_ok:
        orders: set[int] = set()
        for harmonic in harmonics:
            if not isinstance(harmonic, Mapping) or set(harmonic) != {"order", "amplitude_nm", "phase_deg"}:
                harmonics_ok = False
                break
            order = harmonic["order"]
            if not isinstance(order, int) or isinstance(order, bool) or order <= 0 or order in orders or not _finite(harmonic["amplitude_nm"]) or float(harmonic["amplitude_nm"]) < 0.0 or not _finite(harmonic["phase_deg"]) or not -180.0 <= float(harmonic["phase_deg"]) <= 180.0:
                harmonics_ok = False
                break
            orders.add(order)
    return (
        _generations(row, "slice_generation", "angle_generation", "harmonic_generation", "owner_generation", "result_generation")
        and slices_ok
        and row.get("result_slice_weights") == weights
        and row.get("result_skew_angles_mechanical_deg") == angles
        and harmonics_ok
        and row.get("result_harmonic_torque") == harmonics
        and str(row.get("rotor_owner") or "").startswith("rotor:")
        and row.get("result_rotor_owner") == row.get("rotor_owner")
        and _result(row)
    )


def _demag_ok(row: Mapping[str, object]) -> bool:
    operating = row.get("operating_point")
    recoil = row.get("recoil_line")
    fraction = row.get("irreversible_demag_fraction")
    state = row.get("irreversible_state")
    operating_ok = isinstance(operating, Mapping) and set(operating) == {"b_t", "h_a_per_m"} and _finite(operating["b_t"]) and float(operating["b_t"]) >= 0.0 and _finite(operating["h_a_per_m"]) and float(operating["h_a_per_m"]) < 0.0
    recoil_ok = isinstance(recoil, Mapping) and set(recoil) == {"relative_permeability", "coercivity_a_per_m"} and _finite(recoil["relative_permeability"]) and float(recoil["relative_permeability"]) >= 1.0 and _finite(recoil["coercivity_a_per_m"]) and float(recoil["coercivity_a_per_m"]) > 0.0
    fraction_ok = _finite(fraction) and 0.0 <= float(fraction) <= 1.0
    state_ok = fraction_ok and ((float(fraction) == 0.0 and state == "reversible") or (0.0 < float(fraction) < 1.0 and state == "partially_demagnetized") or (float(fraction) == 1.0 and state == "fully_demagnetized"))
    return (
        _generations(row, "operating_generation", "temperature_generation", "recoil_generation", "irreversible_generation", "owner_generation", "result_generation")
        and operating_ok
        and row.get("result_operating_point") == operating
        and _finite(row.get("temperature_c"))
        and float(row["temperature_c"]) > -273.15
        and row.get("result_temperature_c") == row.get("temperature_c")
        and recoil_ok
        and row.get("result_recoil_line") == recoil
        and state_ok
        and row.get("result_irreversible_demag_fraction") == fraction
        and row.get("result_irreversible_state") == state
        and str(row.get("magnet_owner") or "").startswith("magnet:")
        and row.get("result_magnet_owner") == row.get("magnet_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks = validate_public_v54_identity(identity)
    skew = identity.get(SKEW)
    demag = identity.get(DEMAG)
    if skew is not None:
        checks["motor_v53_skew_weight_angle_harmonic_torque_rotor_owner"] = isinstance(skew, Mapping) and _skew_ok(skew)
    if demag is not None:
        checks["motor_v53_demag_operating_temperature_recoil_irreversible_owner"] = isinstance(demag, Mapping) and _demag_ok(demag)
    return checks
