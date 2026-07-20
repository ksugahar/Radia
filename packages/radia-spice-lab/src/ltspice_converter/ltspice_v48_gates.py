"""Behavioral integration and Fourier artifact identity gates."""

from __future__ import annotations

import math
from collections.abc import Mapping


BEHAVIOR = "behavioral_source_discontinuity_event_timestep_reltol_waveform_owner_identity"
FOURIER = "fourier_harmonic_window_fundamental_phase_reference_trace_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(contract: Mapping[str, object], *names: str) -> bool:
    value = str(contract.get("generation_id") or "")
    return bool(value) and all(contract.get(name) == value for name in names)


def _finite(values: object, *, minimum: int = 1) -> bool:
    return isinstance(values, list) and len(values) >= minimum and all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) for value in values)


def _behavior_ok(contract: Mapping[str, object]) -> bool:
    events = contract.get("discontinuity_events_s")
    times = contract.get("accepted_timesteps_s")
    waveform = contract.get("waveform_values")
    reltol = contract.get("reltol")
    return (
        _generation(contract, "event_generation_id", "timestep_generation_id", "tolerance_generation_id", "waveform_generation_id", "result_generation_id")
        and _finite(events)
        and events == sorted(set(events))
        and contract.get("result_discontinuity_events_s") == events
        and _finite(times, minimum=2)
        and all(float(left) < float(right) for left, right in zip(times, times[1:]))
        and all(event in times for event in events)
        and contract.get("result_accepted_timesteps_s") == times
        and isinstance(reltol, (int, float))
        and 0.0 < float(reltol) <= 0.1
        and contract.get("result_reltol") == reltol
        and _finite(waveform, minimum=len(times))
        and len(waveform) == len(times)
        and contract.get("result_waveform_values") == waveform
        and str(contract.get("waveform_owner") or "").startswith("waveform:")
        and contract.get("result_waveform_owner") == contract.get("waveform_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _fourier_ok(contract: Mapping[str, object]) -> bool:
    window = contract.get("analysis_window_s")
    fundamental = contract.get("fundamental_hz")
    orders = contract.get("harmonic_orders")
    amplitudes = contract.get("harmonic_amplitudes")
    phases = contract.get("harmonic_phase_deg")
    cycles = (float(window[1]) - float(window[0])) * float(fundamental) if _finite(window, minimum=2) and len(window) == 2 and isinstance(fundamental, (int, float)) else 0.0
    return (
        _generation(contract, "window_generation_id", "fundamental_generation_id", "harmonic_generation_id", "phase_generation_id", "trace_generation_id", "result_generation_id")
        and _finite(window, minimum=2)
        and len(window) == 2
        and float(window[0]) < float(window[1])
        and contract.get("result_analysis_window_s") == window
        and isinstance(fundamental, (int, float))
        and math.isfinite(float(fundamental))
        and float(fundamental) > 0.0
        and contract.get("result_fundamental_hz") == fundamental
        and cycles >= 1.0
        and abs(cycles - round(cycles)) <= 1.0e-9
        and isinstance(orders, list)
        and bool(orders)
        and all(isinstance(order, int) and not isinstance(order, bool) and order > 0 for order in orders)
        and orders == sorted(set(orders))
        and contract.get("result_harmonic_orders") == orders
        and _finite(amplitudes, minimum=len(orders))
        and len(amplitudes) == len(orders)
        and all(float(value) >= 0.0 for value in amplitudes)
        and contract.get("result_harmonic_amplitudes") == amplitudes
        and _finite(phases, minimum=len(orders))
        and len(phases) == len(orders)
        and contract.get("result_harmonic_phase_deg") == phases
        and contract.get("phase_reference") == contract.get("result_phase_reference") == "cosine_at_window_start"
        and str(contract.get("trace_owner") or "").startswith("trace:")
        and contract.get("result_trace_owner") == contract.get("trace_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v48_identity(positive: Mapping[str, object]) -> bool:
    if not isinstance(positive, Mapping):
        return False
    behavior = positive.get(BEHAVIOR)
    fourier = positive.get(FOURIER)
    if behavior is None and fourier is None:
        return True
    return isinstance(behavior, Mapping) and isinstance(fourier, Mapping) and _behavior_ok(behavior) and _fourier_ok(fourier)
