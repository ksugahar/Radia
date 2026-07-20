"""Motor-map power balance and induction-machine slip identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping

MAP = "motormap_speed_torque_input_output_loss_efficiency_owner_identity"
INDUCTION = "inductionmotor_slip_synchronousspeed_rotorfrequency_torque_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _number(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    number = float(value)
    return math.isfinite(number) and (not positive or number > 0.0) and (not nonnegative or number >= 0.0)


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(float(left), float(right), rel_tol=1.0e-9, abs_tol=1.0e-9)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _map_ok(row: Mapping[str, object]) -> bool:
    losses = row.get("loss_components_w")
    losses_ok = isinstance(losses, Mapping) and bool(losses) and all(isinstance(name, str) and name.endswith("_w") and _number(value, nonnegative=True) for name, value in losses.items())
    speed = row.get("speed_rpm"); torque = row.get("torque_nm")
    expected_output = float(torque) * float(speed) * 2.0 * math.pi / 60.0 if _number(speed) and _number(torque) else math.nan
    expected_input = expected_output + sum(float(value) for value in losses.values()) if losses_ok else math.nan
    expected_efficiency = expected_output / expected_input if expected_input > 0.0 else math.nan
    names = ("speed_rpm", "torque_nm", "input_power_w", "output_power_w", "loss_components_w", "efficiency")
    return (
        _generations(row, "speed_generation", "torque_generation", "input_generation", "output_generation", "loss_generation", "efficiency_generation", "owner_generation", "result_generation")
        and losses_ok and _number(row.get("input_power_w"), positive=True) and _number(row.get("output_power_w"), nonnegative=True)
        and _close(row.get("output_power_w"), expected_output) and _close(row.get("input_power_w"), expected_input)
        and _number(row.get("efficiency"), nonnegative=True) and float(row["efficiency"]) <= 1.0 and _close(row.get("efficiency"), expected_efficiency)
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and str(row.get("result_owner") or "").startswith("result:") and row.get("accepted_result_owner") == row.get("result_owner") and _result(row)
    )


def _induction_ok(row: Mapping[str, object]) -> bool:
    frequency = row.get("supply_frequency_hz"); poles = row.get("pole_count"); rotor = row.get("rotor_speed_rpm"); torque = row.get("torque_nm")
    base_ok = _number(frequency, positive=True) and isinstance(poles, int) and not isinstance(poles, bool) and poles >= 2 and poles % 2 == 0 and _number(rotor, nonnegative=True)
    expected_sync = 120.0 * float(frequency) / poles if base_ok else math.nan
    expected_slip = (expected_sync - float(rotor)) / expected_sync if base_ok else math.nan
    expected_state = "motoring" if 0.0 < expected_slip < 1.0 and _number(torque, positive=True) else "invalid"
    names = ("supply_frequency_hz", "pole_count", "synchronous_speed_rpm", "rotor_speed_rpm", "slip", "rotor_frequency_hz", "torque_nm", "torque_state")
    return (
        _generations(row, "frequency_generation", "pole_generation", "speed_generation", "slip_generation", "rotorfrequency_generation", "torque_generation", "owner_generation", "result_generation")
        and base_ok and _close(row.get("synchronous_speed_rpm"), expected_sync) and _close(row.get("slip"), expected_slip)
        and _close(row.get("rotor_frequency_hz"), expected_slip * float(frequency)) and expected_state == "motoring" and row.get("torque_state") == expected_state
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and str(row.get("result_owner") or "").startswith("result:") and row.get("accepted_result_owner") == row.get("result_owner") and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if identity.get(MAP) is not None:
        checks["motor_v56_map_power_loss_efficiency_owner"] = isinstance(identity[MAP], Mapping) and _map_ok(identity[MAP])
    if identity.get(INDUCTION) is not None:
        checks["motor_v56_induction_slip_rotorfrequency_torque_owner"] = isinstance(identity[INDUCTION], Mapping) and _induction_ok(identity[INDUCTION])
    return checks
