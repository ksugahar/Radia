from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v26 import _v26


_PROMOTED_CASE_IDS = (
    "v27_public_noise_input_output_source_normalization_psd_sidedness_integration_grid_mismatch",
    "v27_public_switch_hysteresis_event_order_max_timestep_measure_window_waveform_generation_mismatch",
)


def _v27():
    summary = _v26()
    positive = summary["metrics"]["positive"]
    generation = "noise-141"
    positive[
        "noise_input_output_source_normalization_psd_sidedness_integration_grid_generation_identity"
    ] = {
        "noise_generation_id": generation,
        "input_source_noise_generation_id": generation,
        "output_source_noise_generation_id": generation,
        "normalization_noise_generation_id": generation,
        "psd_noise_generation_id": generation,
        "integration_grid_noise_generation_id": generation,
        "result_noise_generation_id": generation,
        "input_source_id": "V1",
        "result_input_source_id": "V1",
        "input_node": "in",
        "result_input_node": "in",
        "output_node": "out",
        "result_output_node": "out",
        "normalization": "input_referred_voltage_density",
        "result_normalization": "input_referred_voltage_density",
        "psd_sidedness": "one_sided",
        "result_psd_sidedness": "one_sided",
        "psd_unit": "V^2/Hz",
        "result_psd_unit": "V^2/Hz",
        "frequency_grid_hz": [10.0, 100.0, 1000.0, 10000.0],
        "result_frequency_grid_hz": [10.0, 100.0, 1000.0, 10000.0],
        "psd_values": [1.0e-18, 8.0e-19, 6.0e-19, 5.0e-19],
        "result_psd_values": [1.0e-18, 8.0e-19, 6.0e-19, 5.0e-19],
        "integration_rule": "log_frequency_trapezoid",
        "result_integration_rule": "log_frequency_trapezoid",
        "integrated_noise_v_rms": 9.0e-8,
        "result_integrated_noise_v_rms": 9.0e-8,
        "noise_result_sha256": "1" * 64,
        "accepted_noise_result_sha256": "1" * 64,
    }
    generation = "switch-event-141"
    positive[
        "switch_hysteresis_event_order_max_timestep_measure_window_waveform_generation_identity"
    ] = {
        "switch_generation_id": generation,
        "hysteresis_switch_generation_id": generation,
        "event_order_switch_generation_id": generation,
        "timestep_switch_generation_id": generation,
        "measure_window_switch_generation_id": generation,
        "waveform_switch_generation_id": generation,
        "result_switch_generation_id": generation,
        "switch_model_id": "voltage_hysteretic_switch",
        "result_switch_model_id": "voltage_hysteretic_switch",
        "threshold_high_v": 3.0,
        "result_threshold_high_v": 3.0,
        "threshold_low_v": 2.0,
        "result_threshold_low_v": 2.0,
        "event_order": ["rising_on", "falling_off", "rising_on"],
        "result_event_order": ["rising_on", "falling_off", "rising_on"],
        "max_timestep_s": 1.0e-7,
        "result_max_timestep_s": 1.0e-7,
        "measure_window_s": [1.0e-6, 20.0e-6],
        "result_measure_window_s": [1.0e-6, 20.0e-6],
        "waveform_sha256": "2" * 64,
        "result_waveform_sha256": "2" * 64,
        "measure_table_sha256": "3" * 64,
        "accepted_measure_table_sha256": "3" * 64,
    }
    return summary


def test_v27_positive_noise_and_switch_identities():
    assert ideal_transformer_identity_gate(_v27())["status"] == "ok"


def test_v27_public_noise_identity_mismatch():
    summary = _v27()
    contract = summary["metrics"]["positive"][
        "noise_input_output_source_normalization_psd_sidedness_integration_grid_generation_identity"
    ]
    contract.update(
        {
            "input_source_noise_generation_id": "noise-140",
            "psd_noise_generation_id": "noise-139",
            "integration_grid_noise_generation_id": "noise-138",
            "result_input_source_id": "I1",
            "result_input_node": "old_in",
            "result_output_node": "old_out",
            "result_normalization": "output_referred_current_density",
            "result_psd_sidedness": "two_sided",
            "result_psd_unit": "A^2/Hz",
            "result_frequency_grid_hz": [10.0, 500.0, 10000.0],
            "result_psd_values": [1.0e-12],
            "result_integration_rule": "linear_rectangle",
            "result_integrated_noise_v_rms": 2.0e-3,
            "accepted_noise_result_sha256": "7" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "noise_integration_uses_current_sources_normalization_sidedness_grid_units_and_result"
    ]


def test_v27_public_switch_event_identity_mismatch():
    summary = _v27()
    contract = summary["metrics"]["positive"][
        "switch_hysteresis_event_order_max_timestep_measure_window_waveform_generation_identity"
    ]
    contract.update(
        {
            "hysteresis_switch_generation_id": "switch-event-140",
            "event_order_switch_generation_id": "switch-event-139",
            "waveform_switch_generation_id": "switch-event-138",
            "result_switch_model_id": "ideal_switch_no_hysteresis",
            "result_threshold_high_v": 2.0,
            "result_threshold_low_v": 3.0,
            "result_event_order": ["falling_off", "rising_on"],
            "result_max_timestep_s": 1.0e-5,
            "result_measure_window_s": [0.0, 0.5e-6],
            "result_waveform_sha256": "8" * 64,
            "accepted_measure_table_sha256": "9" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "switch_timing_uses_current_hysteresis_events_timestep_window_waveform_and_measures"
    ]
