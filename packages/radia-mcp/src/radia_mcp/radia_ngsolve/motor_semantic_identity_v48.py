"""Solver-neutral skew aggregation and PWM timeline identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


SKEW = "skew_slice_weight_rotor_angle_torque_harmonic_phase_aggregation_identity"
PWM = "pwm_carrier_control_sample_switch_state_current_voltage_loss_owner_identity"


def _sha(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _digest(row: Mapping[str, object]) -> bool:
    return _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_vector(value: object, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _skew_ok(row: Mapping[str, object]) -> bool:
    slices = row.get("slice_ids")
    weights = row.get("slice_weights")
    angles = row.get("rotor_angles_deg")
    harmonics = row.get("torque_harmonic_phasors_nm")
    phases = row.get("phase_origins_deg")
    count = len(slices) if isinstance(slices, list) else 0
    return (
        _generation(row, ("slice_generation", "angle_generation", "harmonic_generation", "phase_generation", "result_generation"))
        and count > 0
        and len(set(slices)) == count
        and _finite_vector(weights, count)
        and all(float(weight) >= 0.0 for weight in weights)
        and math.isclose(sum(float(weight) for weight in weights), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and _finite_vector(angles, count)
        and isinstance(harmonics, list)
        and len(harmonics) == count
        and all(_finite_vector(phasor, 2) for phasor in harmonics)
        and _finite_vector(phases, count)
        and row.get("result_slice_ids") == slices
        and row.get("result_slice_weights") == weights
        and row.get("result_rotor_angles_deg") == angles
        and row.get("result_torque_harmonic_phasors_nm") == harmonics
        and row.get("result_phase_origins_deg") == phases
        and str(row.get("machine_state_owner") or "").startswith("machine-state:")
        and row.get("result_machine_state_owner") == row.get("machine_state_owner")
        and _digest(row)
    )


def _matrix_rows(value: object, count: int, width: int) -> bool:
    return isinstance(value, list) and len(value) == count and all(_finite_vector(row, width) for row in value)


def _pwm_ok(row: Mapping[str, object]) -> bool:
    times = row.get("sample_times_s")
    states = row.get("switch_states")
    currents = row.get("phase_current_a")
    voltages = row.get("phase_voltage_v")
    losses = row.get("loss_w")
    count = len(times) if isinstance(times, list) else 0
    return (
        _generation(row, ("carrier_generation", "control_generation", "switch_generation", "electrical_generation", "loss_generation", "result_generation"))
        and _finite_vector(times)
        and all(float(times[index]) < float(times[index + 1]) for index in range(count - 1))
        and isinstance(row.get("carrier_frequency_hz"), (int, float))
        and float(row["carrier_frequency_hz"]) > 0.0
        and isinstance(row.get("control_sample_divider"), int)
        and int(row["control_sample_divider"]) > 0
        and isinstance(states, list)
        and len(states) == count
        and all(isinstance(state, list) and len(state) == 3 and set(state) <= {0, 1} for state in states)
        and _matrix_rows(currents, count, 3)
        and _matrix_rows(voltages, count, 3)
        and _finite_vector(losses, count)
        and all(float(loss) >= 0.0 for loss in losses)
        and row.get("result_sample_times_s") == times
        and row.get("result_carrier_frequency_hz") == row.get("carrier_frequency_hz")
        and row.get("result_control_sample_divider") == row.get("control_sample_divider")
        and row.get("result_switch_states") == states
        and row.get("result_phase_current_a") == currents
        and row.get("result_phase_voltage_v") == voltages
        and row.get("result_loss_w") == losses
        and str(row.get("timeline_owner") or "").startswith("timeline:")
        and row.get("result_timeline_owner") == row.get("timeline_owner")
        and _digest(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    skew = identity.get(SKEW)
    pwm = identity.get(PWM)
    if skew is not None:
        checks["motor_v48_skew_slice_angle_harmonic_phase_owner"] = isinstance(skew, Mapping) and _skew_ok(skew)
    if pwm is not None:
        checks["motor_v48_pwm_timeline_switch_electrical_loss_owner"] = isinstance(pwm, Mapping) and _pwm_ok(pwm)
    return checks
