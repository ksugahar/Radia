"""Iron-loss waveform and cogging-period identity checks for v52."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


IRON_LOSS = "ironloss_fft_window_harmonic_rotation_frequency_coefficient_waveform_owner_identity"
COGGING = "cogging_slot_pole_period_sampling_phase_torque_owner_identity"


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
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and bool(value)
        and all(_finite(item) for item in value)
    )


def _iron_loss_ok(row: Mapping[str, object]) -> bool:
    sample_count = row.get("sample_count")
    harmonics = row.get("harmonics")
    rotation_frequency = row.get("rotation_frequency_hz")
    pole_pairs = row.get("pole_pairs")
    electrical_frequency = row.get("electrical_frequency_hz")
    coefficients = row.get("loss_coefficients")
    sample_ok = (
        isinstance(sample_count, int)
        and not isinstance(sample_count, bool)
        and sample_count >= 64
        and sample_count & (sample_count - 1) == 0
    )
    harmonics_ok = (
        isinstance(harmonics, list)
        and bool(harmonics)
        and all(
            isinstance(item, Mapping)
            and isinstance(item.get("order"), int)
            and not isinstance(item.get("order"), bool)
            and 1 <= item["order"] <= sample_count // 2
            and _finite(item.get("amplitude_t"))
            and float(item["amplitude_t"]) >= 0.0
            for item in harmonics
        ) if sample_ok else False
    )
    coefficients_ok = (
        isinstance(coefficients, Mapping)
        and set(coefficients) == {"hysteresis", "classical_eddy", "excess"}
        and all(_finite(value) and float(value) >= 0.0 for value in coefficients.values())
    )
    return (
        _generations(row, "window_generation", "harmonic_generation", "frequency_generation", "coefficient_generation", "owner_generation", "result_generation")
        and row.get("fft_window") == "hann_periodic"
        and row.get("result_fft_window") == row.get("fft_window")
        and sample_ok
        and row.get("result_sample_count") == sample_count
        and harmonics_ok
        and len({item["order"] for item in harmonics}) == len(harmonics)
        and row.get("result_harmonics") == harmonics
        and _finite(rotation_frequency)
        and float(rotation_frequency) > 0.0
        and row.get("result_rotation_frequency_hz") == rotation_frequency
        and isinstance(pole_pairs, int)
        and not isinstance(pole_pairs, bool)
        and pole_pairs > 0
        and row.get("result_pole_pairs") == pole_pairs
        and _finite(electrical_frequency)
        and math.isclose(float(electrical_frequency), float(rotation_frequency) * pole_pairs, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and row.get("result_electrical_frequency_hz") == electrical_frequency
        and coefficients_ok
        and row.get("result_loss_coefficients") == coefficients
        and str(row.get("waveform_owner") or "").startswith("waveform:")
        and row.get("result_waveform_owner") == row.get("waveform_owner")
        and _result(row)
    )


def _cogging_ok(row: Mapping[str, object]) -> bool:
    slots = row.get("slot_count")
    poles = row.get("pole_count")
    period = row.get("cogging_period_mechanical_deg")
    angles = row.get("sample_angles_mechanical_deg")
    torque = row.get("torque_samples_nm")
    counts_ok = all(isinstance(value, int) and not isinstance(value, bool) and value > 0 for value in (slots, poles))
    expected_period = 360.0 / math.lcm(slots, poles) if counts_ok else None
    angle_ok = (
        _finite_vector(angles)
        and len(angles) >= 5
        and math.isclose(float(angles[0]), 0.0, abs_tol=1.0e-12)
        and math.isclose(float(angles[-1]), float(period), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(float(left) < float(right) for left, right in zip(angles, angles[1:]))
        and all(
            math.isclose(float(right) - float(left), float(angles[1]) - float(angles[0]), rel_tol=1.0e-12, abs_tol=1.0e-12)
            for left, right in zip(angles, angles[1:])
        )
    ) if _finite(period) else False
    torque_ok = (
        _finite_vector(torque)
        and isinstance(angles, Sequence)
        and len(torque) == len(angles)
        and math.isclose(float(torque[0]), float(torque[-1]), rel_tol=1.0e-12, abs_tol=1.0e-12)
    )
    return (
        _generations(row, "period_generation", "sampling_generation", "phase_generation", "owner_generation", "result_generation")
        and counts_ok
        and row.get("result_slot_count") == slots
        and row.get("result_pole_count") == poles
        and _finite(period)
        and math.isclose(float(period), expected_period, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and row.get("result_cogging_period_mechanical_deg") == period
        and angle_ok
        and row.get("result_sample_angles_mechanical_deg") == angles
        and torque_ok
        and row.get("result_torque_samples_nm") == torque
        and row.get("phase_alignment") == "slot_center_to_pole_center"
        and row.get("result_phase_alignment") == row.get("phase_alignment")
        and str(row.get("torque_owner") or "").startswith("torque:")
        and row.get("result_torque_owner") == row.get("torque_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    iron_loss = identity.get(IRON_LOSS)
    cogging = identity.get(COGGING)
    if iron_loss is not None:
        checks["motor_v52_ironloss_fft_harmonics_frequency_coefficients_owner"] = (
            isinstance(iron_loss, Mapping) and _iron_loss_ok(iron_loss)
        )
    if cogging is not None:
        checks["motor_v52_cogging_slot_pole_period_sampling_phase_owner"] = (
            isinstance(cogging, Mapping) and _cogging_ok(cogging)
        )
    return checks
