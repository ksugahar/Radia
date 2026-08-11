from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v36 import _v36

_PROMOTED_CASE_IDS = (
    "v37_public_oscillator_startup_barkhausen_frequency_amplitude_limitcycle_energy_owner_mismatch",
    "v37_public_conducted_emi_lisn_spectrum_window_band_limit_power_owner_mismatch",
)


def _v37():
    payload = _v36()
    positive = payload["metrics"]["positive"]

    generation = "oscillator-closure-237"
    positive[
        "oscillator_startup_barkhausen_frequency_amplitude_limitcycle_energy_owner_result_identity"
    ] = {
        "oscillator_generation_id": generation,
        **{
            key: generation
            for key in (
                "startup_oscillator_generation_id",
                "barkhausen_oscillator_generation_id",
                "frequency_oscillator_generation_id",
                "amplitude_oscillator_generation_id",
                "limitcycle_oscillator_generation_id",
                "energy_oscillator_generation_id",
                "waveform_oscillator_generation_id",
                "circuit_oscillator_generation_id",
                "result_oscillator_generation_id",
            )
        },
        "startup_time_s": [0.0, 1.0e-4, 2.0e-4, 3.0e-4],
        "result_startup_time_s": [0.0, 1.0e-4, 2.0e-4, 3.0e-4],
        "startup_amplitude_v": [0.01, 0.02, 0.04, 0.08],
        "result_startup_amplitude_v": [0.01, 0.02, 0.04, 0.08],
        "barkhausen_frequency_hz": 1.0e5,
        "result_barkhausen_frequency_hz": 1.0e5,
        "loop_gain_magnitude": 1.0,
        "result_loop_gain_magnitude": 1.0,
        "loop_phase_deg": 0.0,
        "result_loop_phase_deg": 0.0,
        "steady_frequency_hz": 1.0e5,
        "result_steady_frequency_hz": 1.0e5,
        "steady_amplitude_v": 2.5,
        "result_steady_amplitude_v": 2.5,
        "limitcycle_start_state": [2.5, 0.0],
        "limitcycle_end_state": [2.5, 0.0],
        "result_limitcycle_start_state": [2.5, 0.0],
        "result_limitcycle_end_state": [2.5, 0.0],
        "limitcycle_state_tolerance": 1.0e-9,
        "result_limitcycle_state_tolerance": 1.0e-9,
        "cycle_source_energy_j": 1.2e-6,
        "result_cycle_source_energy_j": 1.2e-6,
        "cycle_dissipated_energy_j": 1.2e-6,
        "result_cycle_dissipated_energy_j": 1.2e-6,
        "waveform_owner": "oscillator/waveform-237",
        "accepted_waveform_owner": "oscillator/waveform-237",
        "waveform_sha256": "a" * 64,
        "accepted_waveform_sha256": "a" * 64,
        "circuit_sha256": "b" * 64,
        "accepted_circuit_sha256": "b" * 64,
        "result_sha256": "c" * 64,
        "accepted_result_sha256": "c" * 64,
    }

    generation = "conducted-emi-237"
    positive[
        "conducted_emi_lisn_spectrum_window_band_limit_power_owner_result_identity"
    ] = {
        "emi_generation_id": generation,
        **{
            key: generation
            for key in (
                "lisn_emi_generation_id",
                "waveform_emi_generation_id",
                "fft_emi_generation_id",
                "detector_emi_generation_id",
                "limit_emi_generation_id",
                "power_emi_generation_id",
                "circuit_emi_generation_id",
                "result_emi_generation_id",
            )
        },
        "lisn_impedance_ohm": 50.0,
        "result_lisn_impedance_ohm": 50.0,
        "waveform_window": "hann",
        "result_waveform_window": "hann",
        "fft_normalization": "rms_single_sided",
        "result_fft_normalization": "rms_single_sided",
        "detector": "quasi_peak",
        "result_detector": "quasi_peak",
        "resolution_bandwidth_hz": 9.0e3,
        "result_resolution_bandwidth_hz": 9.0e3,
        "frequency_hz": [1.5e5, 1.0e6, 3.0e7],
        "result_frequency_hz": [1.5e5, 1.0e6, 3.0e7],
        "spectrum_dbuv": [45.0, 48.0, 50.0],
        "result_spectrum_dbuv": [45.0, 48.0, 50.0],
        "limit_dbuv": [50.0, 53.0, 55.0],
        "result_limit_dbuv": [50.0, 53.0, 55.0],
        "margin_db": [5.0, 5.0, 5.0],
        "result_margin_db": [5.0, 5.0, 5.0],
        "source_input_power_w": 100.0,
        "result_source_input_power_w": 100.0,
        "load_power_w": 80.0,
        "result_load_power_w": 80.0,
        "lisn_loss_w": 5.0,
        "result_lisn_loss_w": 5.0,
        "switching_loss_w": 15.0,
        "result_switching_loss_w": 15.0,
        "circuit_owner": "emi/circuit-237",
        "accepted_circuit_owner": "emi/circuit-237",
        "circuit_sha256": "d" * 64,
        "accepted_circuit_sha256": "d" * 64,
        "result_sha256": "e" * 64,
        "accepted_result_sha256": "e" * 64,
    }
    return payload


def test_v37_public_positive_oscillator_and_conducted_emi_closure():
    assert ideal_transformer_identity_gate(_v37())["status"] == "ok"


def test_v37_public_oscillator_startup_barkhausen_frequency_amplitude_limitcycle_energy_owner_mismatch():
    payload = _v37()
    identity = payload["metrics"]["positive"][
        "oscillator_startup_barkhausen_frequency_amplitude_limitcycle_energy_owner_result_identity"
    ]
    identity.update(
        {
            "startup_oscillator_generation_id": "oscillator-closure-236",
            "result_startup_amplitude_v": [0.08, 0.04, 0.02, 0.01],
            "result_loop_gain_magnitude": 1.4,
            "result_loop_phase_deg": 90.0,
            "result_steady_frequency_hz": 8.0e4,
            "result_steady_amplitude_v": -2.5,
            "result_limitcycle_end_state": [-2.5, 1.0],
            "result_cycle_dissipated_energy_j": 4.0e-7,
            "accepted_waveform_owner": "oscillator/old-waveform",
            "accepted_result_sha256": "5" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "oscillators_use_current_startup_barkhausen_frequency_amplitude_limitcycle_energy_waveform_and_result"
    ]


def test_v37_public_conducted_emi_lisn_spectrum_window_band_limit_power_owner_mismatch():
    payload = _v37()
    identity = payload["metrics"]["positive"][
        "conducted_emi_lisn_spectrum_window_band_limit_power_owner_result_identity"
    ]
    identity.update(
        {
            "fft_emi_generation_id": "conducted-emi-236",
            "result_lisn_impedance_ohm": 75.0,
            "result_waveform_window": "rectangular",
            "result_fft_normalization": "peak_double_sided",
            "result_detector": "average",
            "result_resolution_bandwidth_hz": 1.0e6,
            "result_frequency_hz": [3.0e7, 1.5e5],
            "result_spectrum_dbuv": [70.0, 70.0, 70.0],
            "result_margin_db": [-20.0, -17.0, -15.0],
            "result_source_input_power_w": 80.0,
            "accepted_circuit_owner": "emi/old-circuit",
            "accepted_result_sha256": "6" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "conducted_emi_uses_current_lisn_window_fft_detector_band_limit_power_circuit_and_result"
    ]


def test_v37_public_rejects_self_consistent_nonexponential_startup():
    payload = _v37()
    identity = payload["metrics"]["positive"][
        "oscillator_startup_barkhausen_frequency_amplitude_limitcycle_energy_owner_result_identity"
    ]
    identity["startup_amplitude_v"] = identity["result_startup_amplitude_v"] = [
        0.01,
        0.02,
        0.025,
        0.08,
    ]
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_emi_limit_violation():
    payload = _v37()
    identity = payload["metrics"]["positive"][
        "conducted_emi_lisn_spectrum_window_band_limit_power_owner_result_identity"
    ]
    identity["spectrum_dbuv"] = identity["result_spectrum_dbuv"] = [60.0, 60.0, 60.0]
    identity["margin_db"] = identity["result_margin_db"] = [-10.0, -7.0, -5.0]
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
