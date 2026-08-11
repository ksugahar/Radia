"""Noise and sampled steady-state replay identity gates."""

from __future__ import annotations

import math
from collections.abc import Mapping

from .ltspice_v54_gates import validate_ltspice_v54_identity


NOISE = "noise_spectraldensity_integrated_bandwidth_source_owner_identity"
STEADY = "switchmode_sampled_steadystate_cycle_phase_ripple_waveform_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(contract: Mapping[str, object], *names: str) -> bool:
    generation = str(contract.get("generation_id") or "")
    return bool(generation) and all(contract.get(name) == generation for name in names)


def _number(value: object, *, nonnegative: bool = False, positive: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(float(value)):
        return False
    return (not nonnegative or float(value) >= 0.0) and (not positive or float(value) > 0.0)


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-15)


def _numeric_vector(value: object, *, positive: bool = False, nonnegative: bool = False) -> list[float] | None:
    if not isinstance(value, list) or not value or not all(_number(item, positive=positive, nonnegative=nonnegative) for item in value):
        return None
    return [float(item) for item in value]


def _noise_ok(contract: Mapping[str, object]) -> bool:
    frequency = _numeric_vector(contract.get("frequency_hz"), positive=True)
    density = _numeric_vector(contract.get("spectral_density_v2_per_hz"), nonnegative=True)
    bandwidth = _numeric_vector(contract.get("integration_bandwidth_hz"), nonnegative=True)
    sources = contract.get("contributing_source_fraction")
    if not (
        frequency and density and len(frequency) == len(density) >= 2
        and all(right > left for left, right in zip(frequency, frequency[1:]))
        and bandwidth and len(bandwidth) == 2 and bandwidth[0] == frequency[0] and bandwidth[1] == frequency[-1]
        and isinstance(sources, Mapping) and bool(sources)
        and all(isinstance(name, str) and name.startswith("source:") and _number(value, nonnegative=True) for name, value in sources.items())
        and _close(sum(float(value) for value in sources.values()), 1.0)
    ):
        return False
    variance = sum(
        0.5 * (right_density + left_density) * (right_frequency - left_frequency)
        for left_frequency, right_frequency, left_density, right_density in zip(frequency, frequency[1:], density, density[1:])
    )
    integrated = contract.get("integrated_noise_v_rms")
    return (
        _generation(contract, "spectral_generation_id", "bandwidth_generation_id", "source_generation_id", "owner_generation_id", "result_generation_id")
        and _number(integrated, nonnegative=True) and _close(integrated, math.sqrt(variance))
        and all(contract.get("result_" + name) == contract.get(name) for name in ("frequency_hz", "spectral_density_v2_per_hz", "integration_bandwidth_hz", "integrated_noise_v_rms", "contributing_source_fraction"))
        and str(contract.get("trace_owner") or "").startswith("trace:")
        and contract.get("result_trace_owner") == contract.get("trace_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _steady_ok(contract: Mapping[str, object]) -> bool:
    cycles = contract.get("cycle_index")
    phase = contract.get("sample_phase_fraction")
    waveforms = contract.get("cycle_waveform_v")
    sampled = _numeric_vector(contract.get("sampled_output_v"))
    ripple = _numeric_vector(contract.get("ripple_pp_v"), nonnegative=True)
    if not (
        isinstance(cycles, list) and len(cycles) >= 3
        and all(isinstance(item, int) and not isinstance(item, bool) for item in cycles)
        and all(right == left + 1 for left, right in zip(cycles, cycles[1:]))
        and _number(phase) and 0.0 <= float(phase) <= 1.0
        and isinstance(waveforms, list) and len(waveforms) == len(cycles)
        and sampled and ripple and len(sampled) == len(ripple) == len(cycles)
    ):
        return False
    rows = [_numeric_vector(row) for row in waveforms]
    if not rows or not all(row and len(row) == len(rows[0]) >= 2 for row in rows):
        return False
    position = float(phase) * (len(rows[0]) - 1)
    sample_index = round(position)
    if not math.isclose(position, sample_index, abs_tol=1.0e-12):
        return False
    expected_sampled = [row[sample_index] for row in rows]
    expected_ripple = [max(row) - min(row) for row in rows]
    steady = max(sampled) - min(sampled) <= 0.1 * max(ripple)
    return (
        _generation(contract, "cycle_generation_id", "phase_generation_id", "ripple_generation_id", "waveform_generation_id", "owner_generation_id", "result_generation_id")
        and all(_close(actual, expected) for actual, expected in zip(sampled, expected_sampled))
        and all(_close(actual, expected) for actual, expected in zip(ripple, expected_ripple))
        and steady
        and all(contract.get("result_" + name) == contract.get(name) for name in ("cycle_index", "sample_phase_fraction", "cycle_waveform_v", "sampled_output_v", "ripple_pp_v"))
        and str(contract.get("waveform_owner") or "").startswith("waveform:")
        and contract.get("result_waveform_owner") == contract.get("waveform_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v53_identity(positive: Mapping[str, object]) -> bool:
    """Validate optional v53 identities, requiring the pair when either is present."""
    if not isinstance(positive, Mapping):
        return False
    noise = positive.get(NOISE)
    steady = positive.get(STEADY)
    if noise is None and steady is None:
        return validate_ltspice_v54_identity(positive)
    return isinstance(noise, Mapping) and isinstance(steady, Mapping) and _noise_ok(noise) and _steady_ok(steady) and validate_ltspice_v54_identity(positive)
