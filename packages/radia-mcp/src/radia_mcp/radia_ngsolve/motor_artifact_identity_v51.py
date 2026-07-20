"""Torque-ripple and winding-loss artifact identity checks for v51."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .motor_artifact_identity_v52 import validate_public_identity as validate_public_v52_identity


TORQUE = "torque_ripple_rotor_angle_electrical_mechanical_period_fft_window_owner_identity"
WINDING = "winding_temperature_resistance_endturn_length_fillfactor_copperloss_owner_identity"


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


def _finite_vector(value: object) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and bool(value) and all(_finite(item) for item in value)


def _torque_ok(row: Mapping[str, object]) -> bool:
    pole_pairs = row.get("pole_pairs")
    mechanical = row.get("rotor_angles_mechanical_deg")
    electrical = row.get("rotor_angles_electrical_deg")
    mechanical_period = row.get("mechanical_period_deg")
    electrical_period = row.get("electrical_period_deg")
    window = row.get("sample_window_mechanical_deg")
    torque = row.get("torque_samples_nm")
    harmonics = row.get("fft_harmonics")
    sample_count = len(torque) if isinstance(torque, Sequence) and not isinstance(torque, (str, bytes)) else 0
    angle_ok = (
        isinstance(pole_pairs, int)
        and not isinstance(pole_pairs, bool)
        and pole_pairs > 0
        and row.get("result_pole_pairs") == pole_pairs
        and _finite_vector(mechanical)
        and len(mechanical) >= 5
        and all(float(left) < float(right) for left, right in zip(mechanical, mechanical[1:]))
        and row.get("result_rotor_angles_mechanical_deg") == mechanical
        and _finite_vector(electrical)
        and len(electrical) == len(mechanical)
        and all(math.isclose(float(e), float(m) * pole_pairs, rel_tol=1e-12, abs_tol=1e-12) for m, e in zip(mechanical, electrical))
        and row.get("result_rotor_angles_electrical_deg") == electrical
        and _finite(mechanical_period)
        and float(mechanical_period) > 0.0
        and row.get("result_mechanical_period_deg") == mechanical_period
        and _finite(electrical_period)
        and math.isclose(float(electrical_period), float(mechanical_period) * pole_pairs, rel_tol=1e-12, abs_tol=1e-12)
        and row.get("result_electrical_period_deg") == electrical_period
        and math.isclose(float(mechanical[0]), 0.0, abs_tol=1e-12)
        and math.isclose(float(mechanical[-1]), float(mechanical_period), rel_tol=1e-12, abs_tol=1e-12)
        and isinstance(window, Sequence)
        and not isinstance(window, (str, bytes))
        and len(window) == 2
        and list(window) == [mechanical[0], mechanical[-1]]
        and row.get("result_sample_window_mechanical_deg") == window
    )
    torque_ok = (
        _finite_vector(torque)
        and isinstance(mechanical, Sequence)
        and len(torque) == len(mechanical)
        and math.isclose(float(torque[0]), float(torque[-1]), rel_tol=1e-12, abs_tol=1e-12)
        and row.get("result_torque_samples_nm") == torque
    )
    harmonic_ok = (
        isinstance(harmonics, list)
        and bool(harmonics)
        and all(
            isinstance(item, Mapping)
            and isinstance(item.get("order"), int)
            and not isinstance(item.get("order"), bool)
            and sample_count >= 3
            and 0 <= item["order"] <= (sample_count - 1) // 2
            and _finite(item.get("amplitude_nm"))
            and float(item["amplitude_nm"]) >= 0.0
            for item in harmonics
        )
        and len({item["order"] for item in harmonics}) == len(harmonics)
        and any(item["order"] == 0 for item in harmonics)
        and row.get("result_fft_harmonics") == harmonics
    )
    return (
        _generations(row, "rotor_generation", "angle_generation", "period_generation", "window_generation", "fft_generation", "owner_generation", "result_generation")
        and angle_ok
        and torque_ok
        and harmonic_ok
        and str(row.get("torque_owner") or "").startswith("torque:")
        and row.get("result_torque_owner") == row.get("torque_owner")
        and _result(row)
    )


def _winding_ok(row: Mapping[str, object]) -> bool:
    reference_temperature = row.get("reference_temperature_c")
    temperature = row.get("winding_temperature_c")
    coefficient = row.get("copper_temperature_coefficient_per_k")
    reference_resistance = row.get("resistance_reference_ohm")
    resistance = row.get("resistance_at_temperature_ohm")
    active_length = row.get("active_length_m")
    end_turn_length = row.get("end_turn_length_m")
    fill_factor = row.get("slot_fill_factor")
    current = row.get("current_rms_a")
    loss = row.get("copper_loss_w")
    numeric = all(_finite(value) for value in (reference_temperature, temperature, coefficient, reference_resistance, resistance, active_length, end_turn_length, fill_factor, current, loss))
    if not numeric:
        return False
    expected_resistance = float(reference_resistance) * (1.0 + float(coefficient) * (float(temperature) - float(reference_temperature)))
    expected_loss = float(current) ** 2 * float(resistance)
    return (
        _generations(row, "temperature_generation", "resistance_generation", "length_generation", "fill_generation", "loss_generation", "owner_generation", "result_generation")
        and float(reference_temperature) > -273.15
        and float(temperature) > -273.15
        and row.get("result_winding_temperature_c") == temperature
        and float(coefficient) > 0.0
        and float(reference_resistance) > 0.0
        and float(resistance) > 0.0
        and math.isclose(float(resistance), expected_resistance, rel_tol=1e-12, abs_tol=1e-15)
        and row.get("result_resistance_at_temperature_ohm") == resistance
        and float(active_length) > 0.0
        and float(end_turn_length) > 0.0
        and row.get("result_end_turn_length_m") == end_turn_length
        and 0.0 < float(fill_factor) <= 1.0
        and row.get("result_slot_fill_factor") == fill_factor
        and float(current) >= 0.0
        and float(loss) >= 0.0
        and math.isclose(float(loss), expected_loss, rel_tol=1e-12, abs_tol=1e-12)
        and row.get("result_copper_loss_w") == loss
        and str(row.get("winding_owner") or "").startswith("winding:")
        and row.get("result_winding_owner") == row.get("winding_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks = validate_public_v52_identity(identity)
    torque = identity.get(TORQUE)
    winding = identity.get(WINDING)
    if torque is not None:
        checks["motor_v51_torque_ripple_angle_period_fft_window_owner"] = isinstance(torque, Mapping) and _torque_ok(torque)
    if winding is not None:
        checks["motor_v51_winding_temperature_resistance_length_fill_loss_owner"] = isinstance(winding, Mapping) and _winding_ok(winding)
    return checks
