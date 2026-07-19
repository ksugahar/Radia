"""Motor demagnetization and iron-loss artifact identity checks for v49."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


DEMAG = "demag_temperature_current_angle_irreversible_magnet_state_owner_identity"
IRON = "iron_loss_harmonic_time_frequency_hysteresis_eddy_excess_coefficient_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_vector(value: object, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _magnet_state(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(str(owner).startswith("magnet:") for owner in value)
        and all(isinstance(state, (int, float)) and math.isfinite(float(state)) and 0.0 <= float(state) <= 1.0 for state in value.values())
    )


def _demag_ok(row: Mapping[str, object]) -> bool:
    temperature = row.get("temperature_c")
    currents = row.get("phase_current_a")
    angle = row.get("rotor_angle_electrical_deg")
    state = row.get("irreversible_magnet_state")
    owner = str(row.get("operating_point_owner") or "")
    return (
        _generations(row, "temperature_generation", "current_generation", "angle_generation", "state_generation", "result_generation")
        and isinstance(temperature, (int, float))
        and math.isfinite(float(temperature))
        and -273.15 < float(temperature) <= 500.0
        and row.get("result_temperature_c") == temperature
        and _finite_vector(currents, 3)
        and math.isclose(sum(float(current) for current in currents), 0.0, rel_tol=0.0, abs_tol=1.0e-9)
        and row.get("result_phase_current_a") == currents
        and isinstance(angle, (int, float))
        and math.isfinite(float(angle))
        and row.get("result_rotor_angle_electrical_deg") == angle
        and _magnet_state(state)
        and row.get("result_irreversible_magnet_state") == state
        and owner.startswith("operating-point:")
        and row.get("result_operating_point_owner") == owner
        and _result(row)
    )


def _harmonics(value: object) -> bool:
    if not isinstance(value, list) or not value or not all(isinstance(row, Mapping) for row in value):
        return False
    orders: list[int] = []
    fundamentals: list[float] = []
    for row in value:
        order = row.get("order")
        frequency = row.get("frequency_hz")
        if not isinstance(order, int) or order <= 0 or not isinstance(frequency, (int, float)) or not math.isfinite(float(frequency)) or float(frequency) <= 0.0:
            return False
        losses = [row.get("hysteresis_w"), row.get("eddy_w"), row.get("excess_w")]
        if not all(isinstance(loss, (int, float)) and math.isfinite(float(loss)) and float(loss) >= 0.0 for loss in losses):
            return False
        orders.append(order)
        fundamentals.append(float(frequency) / order)
    return len(set(orders)) == len(orders) and all(
        math.isclose(value, fundamentals[0], rel_tol=1.0e-12, abs_tol=1.0e-12) for value in fundamentals
    )


def _coefficients(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and set(value) == {"hysteresis", "eddy", "excess"}
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) and float(item) >= 0.0 for item in value.values())
    )


def _iron_ok(row: Mapping[str, object]) -> bool:
    window = row.get("time_window_s")
    harmonics = row.get("harmonic_rows")
    coefficients = row.get("loss_coefficients")
    owner = str(row.get("loss_owner") or "")
    return (
        _generations(row, "time_generation", "frequency_generation", "harmonic_generation", "coefficient_generation", "result_generation")
        and _finite_vector(window, 2)
        and float(window[1]) > float(window[0])
        and row.get("result_time_window_s") == window
        and _harmonics(harmonics)
        and row.get("result_harmonic_rows") == harmonics
        and _coefficients(coefficients)
        and row.get("result_loss_coefficients") == coefficients
        and owner.startswith("loss-table:")
        and row.get("result_loss_owner") == owner
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    """Validate optional v49 motor operating-point and iron-loss records."""
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    demag = identity.get(DEMAG)
    iron = identity.get(IRON)
    if demag is not None:
        checks["motor_v49_demag_temperature_current_angle_state_owner"] = isinstance(demag, Mapping) and _demag_ok(demag)
    if iron is not None:
        checks["motor_v49_iron_loss_harmonic_time_frequency_coefficient_owner"] = isinstance(iron, Mapping) and _iron_ok(iron)
    return checks
