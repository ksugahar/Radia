from copy import deepcopy

from ltspice_converter.ltspice_v48_gates import BEHAVIOR, FOURIER, validate_ltspice_v48_identity


PROMOTED_CASE_IDS = (
    "v48_public_behavioral_source_discontinuity_event_timestep_reltol_waveform_owner_mismatch",
    "v48_public_fourier_harmonic_window_fundamental_phase_reference_trace_owner_mismatch",
)


def _positive():
    behavior = "behavior-v48"
    fourier = "fourier-v48"
    events = [1.0e-6, 2.0e-6]
    times = [0.0, 0.5e-6, 1.0e-6, 1.5e-6, 2.0e-6, 2.5e-6]
    waveform = [0.0, 0.5, 1.0, 0.5, 0.0, -0.5]
    return {
        BEHAVIOR: {
            "generation_id": behavior, **{key: behavior for key in ("event_generation_id", "timestep_generation_id", "tolerance_generation_id", "waveform_generation_id", "result_generation_id")},
            "discontinuity_events_s": events, "result_discontinuity_events_s": events,
            "accepted_timesteps_s": times, "result_accepted_timesteps_s": times,
            "reltol": 1.0e-3, "result_reltol": 1.0e-3,
            "waveform_values": waveform, "result_waveform_values": waveform,
            "waveform_owner": "waveform:test", "result_waveform_owner": "waveform:test",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        FOURIER: {
            "generation_id": fourier, **{key: fourier for key in ("window_generation_id", "fundamental_generation_id", "harmonic_generation_id", "phase_generation_id", "trace_generation_id", "result_generation_id")},
            "analysis_window_s": [0.01, 0.02], "result_analysis_window_s": [0.01, 0.02],
            "fundamental_hz": 1000.0, "result_fundamental_hz": 1000.0,
            "harmonic_orders": [1, 2, 3], "result_harmonic_orders": [1, 2, 3],
            "harmonic_amplitudes": [1.0, 0.1, 0.02], "result_harmonic_amplitudes": [1.0, 0.1, 0.02],
            "harmonic_phase_deg": [0.0, -90.0, 45.0], "result_harmonic_phase_deg": [0.0, -90.0, 45.0],
            "phase_reference": "cosine_at_window_start", "result_phase_reference": "cosine_at_window_start",
            "trace_owner": "trace:test", "result_trace_owner": "trace:test",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def test_v48_public_replay_identity():
    assert validate_ltspice_v48_identity(_positive()) is True


def test_v48_public_rejects_behavioral_integration_mutation():
    value = deepcopy(_positive())
    value[BEHAVIOR].update({"result_discontinuity_events_s": [1.5e-6], "result_accepted_timesteps_s": [0.0, 1.5e-6], "result_reltol": 1.0e-2, "result_waveform_owner": "waveform:old"})
    assert validate_ltspice_v48_identity(value) is False


def test_v48_public_rejects_fourier_basis_mutation():
    value = deepcopy(_positive())
    value[FOURIER].update({"result_analysis_window_s": [0.0, 0.01], "result_fundamental_hz": 2000.0, "result_harmonic_orders": [2, 1, 3], "result_phase_reference": "sine_at_zero", "result_trace_owner": "trace:old"})
    assert validate_ltspice_v48_identity(value) is False


def test_v48_public_rejects_self_consistent_nonintegral_fourier_window():
    value = deepcopy(_positive())
    value[FOURIER]["analysis_window_s"] = value[FOURIER]["result_analysis_window_s"] = [0.01, 0.0105]
    assert validate_ltspice_v48_identity(value) is False
