"""AC-transfer and transient-event artifact identity gates."""

from __future__ import annotations

import math
from collections.abc import Mapping

from .ltspice_v52_gates import validate_ltspice_v52_identity


AC_TRANSFER = "ac_transfer_source_phase_db_groupdelay_owner_identity"
TRANSIENT_EVENT = "transient_event_interpolation_compression_edge_window_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(contract: Mapping[str, object], *names: str) -> bool:
    generation = str(contract.get("generation_id") or "")
    return bool(generation) and all(contract.get(name) == generation for name in names)


def _finite(values: object, *, length: int | None = None, positive: bool = False) -> bool:
    if not isinstance(values, list) or (length is not None and len(values) != length):
        return False
    return all(
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and (not positive or float(value) > 0.0)
        for value in values
    )


def _same(left: float, right: float) -> bool:
    return math.isclose(float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-12)


def _ac_transfer_ok(contract: Mapping[str, object]) -> bool:
    frequencies = contract.get("frequency_hz")
    magnitude = contract.get("transfer_magnitude")
    magnitude_db = contract.get("transfer_magnitude_db")
    phase = contract.get("unwrapped_phase_rad")
    group_delay = contract.get("group_delay_s")
    count = len(frequencies) if isinstance(frequencies, list) else 0
    numeric_ok = (
        _finite(frequencies, positive=True)
        and count >= 2
        and all(float(left) < float(right) for left, right in zip(frequencies, frequencies[1:]))
        and _finite(magnitude, length=count, positive=True)
        and _finite(magnitude_db, length=count)
        and _finite(phase, length=count)
        and _finite(group_delay, length=count - 1)
    )
    if not numeric_ok:
        return False
    expected_db = [20.0 * math.log10(float(value)) for value in magnitude]
    expected_delay = [
        -(float(right) - float(left)) / (2.0 * math.pi * (float(f_right) - float(f_left)))
        for left, right, f_left, f_right in zip(phase, phase[1:], frequencies, frequencies[1:])
    ]
    source_amplitude = contract.get("source_amplitude_v")
    source_phase = contract.get("source_phase_deg")
    return (
        _generation(
            contract,
            "source_generation_id", "db_generation_id", "phase_generation_id",
            "group_delay_generation_id", "trace_generation_id", "owner_generation_id",
            "result_generation_id",
        )
        and isinstance(source_amplitude, (int, float))
        and not isinstance(source_amplitude, bool)
        and math.isfinite(float(source_amplitude))
        and float(source_amplitude) > 0.0
        and isinstance(source_phase, (int, float))
        and not isinstance(source_phase, bool)
        and math.isfinite(float(source_phase))
        and contract.get("result_source_amplitude_v") == source_amplitude
        and contract.get("result_source_phase_deg") == source_phase
        and contract.get("db_convention") == "20log10_voltage_ratio"
        and contract.get("result_db_convention") == contract.get("db_convention")
        and contract.get("phase_unwrap") == "continuous_radians"
        and contract.get("result_phase_unwrap") == contract.get("phase_unwrap")
        and contract.get("result_frequency_hz") == frequencies
        and contract.get("result_transfer_magnitude") == magnitude
        and contract.get("result_transfer_magnitude_db") == magnitude_db
        and contract.get("result_unwrapped_phase_rad") == phase
        and contract.get("result_group_delay_s") == group_delay
        and all(_same(actual, expected) for actual, expected in zip(magnitude_db, expected_db))
        and all(_same(actual, expected) for actual, expected in zip(group_delay, expected_delay))
        and str(contract.get("trace_owner") or "").startswith("trace:")
        and contract.get("result_trace_owner") == contract.get("trace_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _transient_event_ok(contract: Mapping[str, object]) -> bool:
    time_s = contract.get("time_s")
    waveform = contract.get("waveform")
    edge = contract.get("edge_identity")
    window = contract.get("measurement_window_s")
    count = len(time_s) if isinstance(time_s, list) else 0
    if not (
        _finite(time_s)
        and count >= 2
        and all(float(left) < float(right) for left, right in zip(time_s, time_s[1:]))
        and _finite(waveform, length=count)
        and isinstance(edge, Mapping)
        and edge.get("direction") in {"rising", "falling"}
        and isinstance(edge.get("threshold"), (int, float))
        and isinstance(edge.get("occurrence"), int)
        and int(edge["occurrence"]) >= 1
        and _finite(window, length=2)
        and float(window[0]) < float(window[1])
    ):
        return False
    threshold = float(edge["threshold"])
    direction = str(edge["direction"])
    crossings: list[float] = []
    for index in range(count - 1):
        left, right = float(waveform[index]), float(waveform[index + 1])
        selected = left < threshold <= right if direction == "rising" else left > threshold >= right
        if selected and right != left:
            fraction = (threshold - left) / (right - left)
            crossings.append(float(time_s[index]) + fraction * (float(time_s[index + 1]) - float(time_s[index])))
    occurrence = int(edge["occurrence"])
    if len(crossings) < occurrence:
        return False
    event_time = contract.get("event_time_s")
    return (
        _generation(
            contract,
            "event_generation_id", "interpolation_generation_id", "compression_generation_id",
            "edge_generation_id", "window_generation_id", "waveform_generation_id",
            "owner_generation_id", "result_generation_id",
        )
        and contract.get("event_interpolation") == "linear_crossing"
        and contract.get("result_event_interpolation") == contract.get("event_interpolation")
        and contract.get("compression_mode") == "disabled"
        and contract.get("result_compression_mode") == contract.get("compression_mode")
        and contract.get("result_time_s") == time_s
        and contract.get("result_waveform") == waveform
        and contract.get("result_edge_identity") == edge
        and contract.get("result_measurement_window_s") == window
        and isinstance(event_time, (int, float))
        and _same(float(event_time), crossings[occurrence - 1])
        and float(window[0]) <= float(event_time) <= float(window[1])
        and contract.get("result_event_time_s") == event_time
        and str(contract.get("waveform_owner") or "").startswith("waveform:")
        and contract.get("result_waveform_owner") == contract.get("waveform_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v51_identity(positive: Mapping[str, object]) -> bool:
    if not isinstance(positive, Mapping):
        return False
    ac_transfer = positive.get(AC_TRANSFER)
    transient_event = positive.get(TRANSIENT_EVENT)
    if ac_transfer is None and transient_event is None:
        return validate_ltspice_v52_identity(positive)
    return (
        isinstance(ac_transfer, Mapping)
        and isinstance(transient_event, Mapping)
        and _ac_transfer_ok(ac_transfer)
        and _transient_event_ok(transient_event)
        and validate_ltspice_v52_identity(positive)
    )
