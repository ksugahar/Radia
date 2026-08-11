from copy import deepcopy
import math

from radia.ltspice.ltspice_v51_gates import (
    AC_TRANSFER,
    TRANSIENT_EVENT,
    validate_ltspice_v51_identity,
)


CASE_IDS = {
    "v51_public_ac_transfer_source_amplitude_phase_db_convention_groupdelay_owner_mismatch",
    "v51_public_transient_event_interpolation_compression_edge_measure_window_owner_mismatch",
}


def _generation(prefix: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation_id": prefix, **{name: prefix for name in names}}


def _positive():
    frequencies = [10.0, 100.0, 1000.0]
    magnitude = [0.999, 0.9, 0.5]
    phase = [-0.01, -0.1, -0.8]
    delays = [
        -(right - left) / (2.0 * math.pi * (f_right - f_left))
        for left, right, f_left, f_right in zip(phase, phase[1:], frequencies, frequencies[1:])
    ]
    ac = {
        **_generation("ac-v51-test", (
            "source_generation_id", "db_generation_id", "phase_generation_id",
            "group_delay_generation_id", "trace_generation_id", "owner_generation_id",
            "result_generation_id",
        )),
        "source_amplitude_v": 1.0, "result_source_amplitude_v": 1.0,
        "source_phase_deg": 0.0, "result_source_phase_deg": 0.0,
        "db_convention": "20log10_voltage_ratio", "result_db_convention": "20log10_voltage_ratio",
        "phase_unwrap": "continuous_radians", "result_phase_unwrap": "continuous_radians",
        "frequency_hz": frequencies, "result_frequency_hz": frequencies,
        "transfer_magnitude": magnitude, "result_transfer_magnitude": magnitude,
        "transfer_magnitude_db": [20.0 * math.log10(value) for value in magnitude],
        "result_transfer_magnitude_db": [20.0 * math.log10(value) for value in magnitude],
        "unwrapped_phase_rad": phase, "result_unwrapped_phase_rad": phase,
        "group_delay_s": delays, "result_group_delay_s": delays,
        "trace_owner": "trace:ac-v51-test", "result_trace_owner": "trace:ac-v51-test",
        "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
    }
    edge = {"signal": "V(out)", "threshold": 0.5, "direction": "rising", "occurrence": 1}
    transient = {
        **_generation("tran-v51-test", (
            "event_generation_id", "interpolation_generation_id", "compression_generation_id",
            "edge_generation_id", "window_generation_id", "waveform_generation_id",
            "owner_generation_id", "result_generation_id",
        )),
        "time_s": [0.0, 1.0, 2.0], "result_time_s": [0.0, 1.0, 2.0],
        "waveform": [0.0, 0.25, 0.75], "result_waveform": [0.0, 0.25, 0.75],
        "event_interpolation": "linear_crossing", "result_event_interpolation": "linear_crossing",
        "compression_mode": "disabled", "result_compression_mode": "disabled",
        "edge_identity": edge, "result_edge_identity": edge,
        "measurement_window_s": [0.0, 2.0], "result_measurement_window_s": [0.0, 2.0],
        "event_time_s": 1.5, "result_event_time_s": 1.5,
        "waveform_owner": "waveform:tran-v51-test", "result_waveform_owner": "waveform:tran-v51-test",
        "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
    }
    return {AC_TRANSFER: ac, TRANSIENT_EVENT: transient}


def test_v51_public_replay_identity():
    assert validate_ltspice_v51_identity(_positive()) is True


def test_v51_public_rejects_ac_transfer_mutation():
    value = deepcopy(_positive())
    value[AC_TRANSFER]["result_db_convention"] = "10log10_power_ratio"
    assert validate_ltspice_v51_identity(value) is False


def test_v51_public_rejects_transient_event_mutation():
    value = deepcopy(_positive())
    value[TRANSIENT_EVENT]["event_time_s"] = 1.0
    value[TRANSIENT_EVENT]["result_event_time_s"] = 1.0
    assert validate_ltspice_v51_identity(value) is False
