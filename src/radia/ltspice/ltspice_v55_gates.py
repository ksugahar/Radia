"""AC-noise and switching-power artifact identity gates."""

from __future__ import annotations

import math
from collections.abc import Mapping


NOISE = "acnoise_density_transfer_bandwidth_owner_identity"
POWER = "switch_average_power_efficiency_window_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _number(value: object, *, nonnegative: bool = False, positive: bool = False) -> bool:
    valid = isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))
    return valid and (not nonnegative or float(value) >= 0.0) and (not positive or float(value) > 0.0)


def _vector(value: object, *, positive: bool = False, nonnegative: bool = False) -> list[float]:
    if not isinstance(value, list) or not value or not all(_number(item, positive=positive, nonnegative=nonnegative) for item in value):
        return []
    return [float(item) for item in value]


def _generation(contract: Mapping[str, object], *fields: str) -> bool:
    generation = str(contract.get("generation_id") or "")
    return bool(generation) and all(contract.get(field) == generation for field in fields)


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1.0e-10, abs_tol=1.0e-15)


def _noise_ok(contract: Mapping[str, object]) -> bool:
    frequency = _vector(contract.get("frequency_hz"), positive=True)
    input_density = _vector(contract.get("input_noise_density_v_per_sqrt_hz"), nonnegative=True)
    gain = _vector(contract.get("transfer_gain_v_per_v"))
    output_density = _vector(contract.get("output_noise_density_v_per_sqrt_hz"), nonnegative=True)
    vectors_ok = len(frequency) >= 2 and len(frequency) == len(input_density) == len(gain) == len(output_density)
    if vectors_ok:
        vectors_ok = all(left < right for left, right in zip(frequency, frequency[1:])) and all(
            _close(output, source * abs(transfer))
            for source, transfer, output in zip(input_density, gain, output_density)
        )
    bandwidth = contract.get("integration_bandwidth_hz")
    integrated = contract.get("integrated_output_noise_v_rms")
    if vectors_ok:
        expected_bandwidth = frequency[-1] - frequency[0]
        variance = sum(
            0.5 * (output_density[index] ** 2 + output_density[index + 1] ** 2)
            * (frequency[index + 1] - frequency[index])
            for index in range(len(frequency) - 1)
        )
        expected_integrated = math.sqrt(variance)
    else:
        expected_bandwidth = math.nan
        expected_integrated = math.nan
    return (
        _generation(contract, "density_generation_id", "gain_generation_id", "bandwidth_generation_id", "grid_generation_id", "owner_generation_id", "result_generation_id")
        and vectors_ok
        and _number(bandwidth, positive=True) and _close(float(bandwidth), expected_bandwidth)
        and _number(integrated, nonnegative=True) and _close(float(integrated), expected_integrated)
        and all(contract.get("result_" + field) == contract.get(field) for field in ("frequency_hz", "input_noise_density_v_per_sqrt_hz", "transfer_gain_v_per_v", "output_noise_density_v_per_sqrt_hz", "integration_bandwidth_hz", "integrated_output_noise_v_rms"))
        and str(contract.get("trace_owner") or "").startswith("trace:") and contract.get("result_trace_owner") == contract.get("trace_owner")
        and _digest(contract.get("result_sha256")) and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _power_ok(contract: Mapping[str, object]) -> bool:
    input_power = contract.get("average_input_power_w")
    output_power = contract.get("average_output_power_w")
    loss_power = contract.get("average_loss_power_w")
    efficiency = contract.get("efficiency")
    window = contract.get("measurement_window")
    window_ok = (
        isinstance(window, Mapping) and set(window) == {"start_s", "stop_s"}
        and _number(window.get("start_s"), nonnegative=True) and _number(window.get("stop_s"), positive=True)
        and float(window["start_s"]) < float(window["stop_s"])
    )
    balance_ok = (
        _number(input_power, positive=True) and _number(output_power, nonnegative=True)
        and _number(loss_power, nonnegative=True) and _number(efficiency, nonnegative=True)
        and float(output_power) <= float(input_power) and 0.0 <= float(efficiency) <= 1.0
        and _close(float(loss_power), float(input_power) - float(output_power))
        and _close(float(efficiency), float(output_power) / float(input_power))
    )
    return (
        _generation(contract, "input_generation_id", "output_generation_id", "loss_generation_id", "efficiency_generation_id", "window_generation_id", "owner_generation_id", "result_generation_id")
        and balance_ok and window_ok
        and all(contract.get("result_" + field) == contract.get(field) for field in ("average_input_power_w", "average_output_power_w", "average_loss_power_w", "efficiency", "measurement_window"))
        and str(contract.get("waveform_owner") or "").startswith("waveform:") and contract.get("result_waveform_owner") == contract.get("waveform_owner")
        and _digest(contract.get("result_sha256")) and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v55_identity(positive: Mapping[str, object]) -> bool:
    """Validate optional v55 identities, requiring the pair when either is present."""
    if not isinstance(positive, Mapping):
        return False
    noise = positive.get(NOISE)
    power = positive.get(POWER)
    if noise is None and power is None:
        return True
    return isinstance(noise, Mapping) and isinstance(power, Mapping) and _noise_ok(noise) and _power_ok(power)
