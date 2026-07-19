from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


_PORT = "v46_public_time_domain_port_wave_impedance_unit_scale_partial_trace_mismatch"
_MONITOR = "v46_public_field_monitor_coordinate_frame_sampling_window_nan_inf_mismatch"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _finite_sequence(value: object, *, minimum: int = 1) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) >= minimum
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _linked(row: Mapping[str, object], fields: tuple[str, ...], generation_key: str = "generation") -> bool:
    generation = str(row.get(generation_key, "")).strip()
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _port_ok(row: Mapping[str, object]) -> bool:
    times = row.get("time_s")
    return (
        _linked(row, ("time_generation", "trace_generation", "impedance_generation", "unit_scale_generation", "partial_generation", "result_generation"))
        and _finite_sequence(times, minimum=2)
        and times == row.get("result_time_s")
        and all(float(left) <= float(right) for left, right in zip(times, times[1:]))
        and _finite_sequence(row.get("port_wave_trace_v"), minimum=2)
        and row.get("port_wave_trace_v") == row.get("result_port_wave_trace_v")
        and row.get("wave_impedance_unit") == row.get("result_wave_impedance_unit") == "ohm"
        and row.get("unit_scale_to_si") == row.get("result_unit_scale_to_si") == 1.0
        and row.get("partial_trace") == row.get("result_partial_trace") is False
        and str(row.get("port_owner", "")).startswith("port:")
        and row.get("result_port_owner") == row.get("port_owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _monitor_ok(row: Mapping[str, object]) -> bool:
    window = row.get("sampling_window_s")
    return (
        _linked(row, ("frame_generation", "sampling_generation", "finite_generation", "monitor_generation", "result_generation"))
        and row.get("coordinate_frame") == row.get("result_coordinate_frame") == "global_cartesian"
        and _finite_sequence(window, minimum=2)
        and window == row.get("result_sampling_window_s")
        and float(window[0]) <= float(window[-1])
        and _finite_sequence(row.get("field_samples"), minimum=1)
        and row.get("field_samples") == row.get("result_field_samples")
        and row.get("nonfinite_sample_count") == row.get("result_nonfinite_sample_count") == 0
        and str(row.get("monitor_owner", "")).startswith("monitor:")
        and row.get("result_monitor_owner") == row.get("monitor_owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def validate_public_v46_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping):
        return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]
    checks: dict[str, bool] = {}
    ports = [row[_PORT] for row in rows if _PORT in row]
    monitors = [row[_MONITOR] for row in rows if _MONITOR in row]
    if ports:
        checks["cst_v46_time_domain_port_identity"] = len(ports) == len(rows) and all(isinstance(row, Mapping) and _port_ok(row) for row in ports)
    if monitors:
        checks["cst_v46_field_monitor_identity"] = len(monitors) == len(rows) and all(isinstance(row, Mapping) and _monitor_ok(row) for row in monitors)
    return checks
