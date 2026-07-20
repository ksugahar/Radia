"""Network renormalization and group-delay artifact checks for v51."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .wave_energy_identity_v52 import validate_public_v52_identity


S_PARAMETER = "sparameter_renormalization_complex_zref_modal_impedance_wavebasis_port_owner_identity"
GROUP_DELAY = "group_delay_unwrap_frequency_derivative_smoothing_window_trace_owner_identity"


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _complex_pair(value: object) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == 2 and all(_finite(item) for item in value)


def _s_parameter_ok(row: Mapping[str, object]) -> bool:
    zref = row.get("complex_reference_impedance_ohm")
    modal = row.get("modal_impedances_ohm")
    return (
        _generation(row, ("renormalization_generation", "zref_generation", "modal_generation", "wave_generation", "owner_generation", "result_generation"))
        and row.get("renormalized") is True
        and row.get("result_renormalized") is True
        and _complex_pair(zref)
        and float(zref[0]) > 0.0
        and row.get("result_complex_reference_impedance_ohm") == zref
        and isinstance(modal, Mapping)
        and bool(modal)
        and all(str(name).startswith("port") and _complex_pair(value) and float(value[0]) > 0.0 for name, value in modal.items())
        and row.get("result_modal_impedances_ohm") == modal
        and row.get("wave_basis") == "power_waves"
        and row.get("result_wave_basis") == row.get("wave_basis")
        and str(row.get("port_owner") or "").startswith("port:")
        and row.get("result_port_owner") == row.get("port_owner")
        and _result(row)
    )


def _group_delay_ok(row: Mapping[str, object]) -> bool:
    frequency = row.get("frequency_hz")
    phase = row.get("unwrapped_phase_deg")
    delay = row.get("group_delay_s")
    window = row.get("smoothing_window_points")
    vectors_ok = (
        isinstance(frequency, list)
        and len(frequency) >= 3
        and all(_finite(value) and float(value) > 0.0 for value in frequency)
        and all(float(left) < float(right) for left, right in zip(frequency, frequency[1:]))
        and isinstance(phase, list)
        and len(phase) == len(frequency)
        and all(_finite(value) for value in phase)
        and isinstance(delay, list)
        and len(delay) == len(frequency)
        and all(_finite(value) for value in delay)
    )
    derivative_ok = False
    if vectors_ok:
        expected = [
            -math.radians(float(right_phase) - float(left_phase)) / (2.0 * math.pi * (float(right_frequency) - float(left_frequency)))
            for left_phase, right_phase, left_frequency, right_frequency in zip(phase, phase[1:], frequency, frequency[1:])
        ]
        derivative_ok = all(math.isclose(float(value), expected[min(index, len(expected) - 1)], rel_tol=1e-9, abs_tol=1e-15) for index, value in enumerate(delay))
    return (
        _generation(row, ("unwrap_generation", "frequency_generation", "derivative_generation", "smoothing_generation", "owner_generation", "result_generation"))
        and vectors_ok
        and row.get("result_frequency_hz") == frequency
        and row.get("result_unwrapped_phase_deg") == phase
        and derivative_ok
        and row.get("result_group_delay_s") == delay
        and row.get("derivative_definition") == "minus_dphi_rad_domega"
        and row.get("result_derivative_definition") == row.get("derivative_definition")
        and isinstance(window, int)
        and not isinstance(window, bool)
        and 3 <= window <= len(frequency)
        and window % 2 == 1
        and row.get("result_smoothing_window_points") == window
        and str(row.get("trace_owner") or "").startswith("trace:")
        and row.get("result_trace_owner") == row.get("trace_owner")
        and _result(row)
    )


def validate_public_v51_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping):
        return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]
    checks = validate_public_v52_identity(payload)
    s_parameters = [row[S_PARAMETER] for row in rows if S_PARAMETER in row]
    delays = [row[GROUP_DELAY] for row in rows if GROUP_DELAY in row]
    if s_parameters:
        checks["network_v51_sparameter_zref_modal_wave_port_owner"] = len(s_parameters) == len(rows) and all(isinstance(row, Mapping) and _s_parameter_ok(row) for row in s_parameters)
    if delays:
        checks["network_v51_group_delay_unwrap_derivative_window_trace_owner"] = len(delays) == len(rows) and all(isinstance(row, Mapping) and _group_delay_ok(row) for row in delays)
    return checks
