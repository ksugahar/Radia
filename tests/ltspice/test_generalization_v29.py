from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v28 import _v28


_PROMOTED_CASE_IDS = (
    "v29_public_mosfet_switching_loss_gate_charge_overlap_deadtime_event_grid_temperature_mismatch",
    "v29_public_step_response_rise_settling_overshoot_initial_state_threshold_window_mismatch",
)


def _v29():
    summary = _v28()
    positive = summary["metrics"]["positive"]
    generation = "switching-loss-161"
    positive[
        "mosfet_switching_loss_gate_charge_overlap_deadtime_event_grid_temperature_cycle_result_generation_identity"
    ] = {
        "switching_generation_id": generation,
        "gate_charge_switching_generation_id": generation,
        "overlap_switching_generation_id": generation,
        "deadtime_switching_generation_id": generation,
        "event_grid_switching_generation_id": generation,
        "temperature_switching_generation_id": generation,
        "cycle_switching_generation_id": generation,
        "result_switching_generation_id": generation,
        "gate_charge_trace_id": "Qgate(M1)",
        "result_gate_charge_trace_id": "Qgate(M1)",
        "overlap_power_trace_id": "Vds(M1)*Id(M1)",
        "result_overlap_power_trace_id": "Vds(M1)*Id(M1)",
        "deadtime_s": 5.0e-8,
        "result_deadtime_s": 5.0e-8,
        "event_times_s": [0.0040001, 0.0040051, 0.0040101, 0.0040151],
        "result_event_times_s": [0.0040001, 0.0040051, 0.0040101, 0.0040151],
        "event_grid_rule": "edge-aligned-local-refinement",
        "result_event_grid_rule": "edge-aligned-local-refinement",
        "junction_temperature_c": 100.0,
        "result_junction_temperature_c": 100.0,
        "cycle_window_s": [0.004, 0.005],
        "result_cycle_window_s": [0.004, 0.005],
        "turn_on_energy_j": 2.0e-6,
        "result_turn_on_energy_j": 2.0e-6,
        "turn_off_energy_j": 3.0e-6,
        "result_turn_off_energy_j": 3.0e-6,
        "switching_frequency_hz": 1.0e5,
        "result_switching_frequency_hz": 1.0e5,
        "switching_loss_w": 0.5,
        "result_switching_loss_w": 0.5,
        "event_grid_sha256": "1" * 64,
        "result_event_grid_sha256": "1" * 64,
        "switching_waveform_sha256": "2" * 64,
        "result_switching_waveform_sha256": "2" * 64,
        "switching_loss_result_sha256": "3" * 64,
        "accepted_switching_loss_result_sha256": "3" * 64,
    }
    generation = "step-response-161"
    positive[
        "step_response_initial_final_rise_threshold_settling_band_overshoot_window_waveform_result_generation_identity"
    ] = {
        "step_generation_id": generation,
        "initial_step_generation_id": generation,
        "final_step_generation_id": generation,
        "rise_step_generation_id": generation,
        "settling_step_generation_id": generation,
        "overshoot_step_generation_id": generation,
        "window_step_generation_id": generation,
        "waveform_step_generation_id": generation,
        "result_step_generation_id": generation,
        "initial_value": 0.0,
        "result_initial_value": 0.0,
        "final_value": 5.0,
        "result_final_value": 5.0,
        "rise_threshold_fractions": [0.1, 0.9],
        "result_rise_threshold_fractions": [0.1, 0.9],
        "rise_crossing_times_s": [1.0e-4, 4.0e-4],
        "result_rise_crossing_times_s": [1.0e-4, 4.0e-4],
        "rise_time_s": 3.0e-4,
        "result_rise_time_s": 3.0e-4,
        "settling_band_fraction": 0.02,
        "result_settling_band_fraction": 0.02,
        "settling_time_s": 1.5e-3,
        "result_settling_time_s": 1.5e-3,
        "overshoot_peak": 5.5,
        "result_overshoot_peak": 5.5,
        "overshoot_fraction": 0.1,
        "result_overshoot_fraction": 0.1,
        "measurement_window_s": [0.0, 0.01],
        "result_measurement_window_s": [0.0, 0.01],
        "waveform_sha256": "4" * 64,
        "result_waveform_sha256": "4" * 64,
        "step_result_sha256": "5" * 64,
        "accepted_step_result_sha256": "5" * 64,
    }
    return summary


def test_v29_positive_switching_loss_and_step_response_identities():
    assert ideal_transformer_identity_gate(_v29())["status"] == "ok"


def test_v29_public_mosfet_switching_loss_gate_charge_overlap_deadtime_event_grid_temperature_mismatch():
    summary = _v29()
    contract = summary["metrics"]["positive"][
        "mosfet_switching_loss_gate_charge_overlap_deadtime_event_grid_temperature_cycle_result_generation_identity"
    ]
    contract.update(
        {
            "gate_charge_switching_generation_id": "switching-loss-160",
            "event_grid_switching_generation_id": "switching-loss-159",
            "result_switching_generation_id": "switching-loss-158",
            "result_gate_charge_trace_id": "V(gate)",
            "result_overlap_power_trace_id": "Vds(M1)",
            "result_deadtime_s": 0.0,
            "result_event_times_s": [0.0, 1.0e-3],
            "result_event_grid_rule": "uniform-coarse",
            "result_junction_temperature_c": 25.0,
            "result_cycle_window_s": [0.0, 1.0e-4],
            "result_turn_on_energy_j": 5.0e-6,
            "result_turn_off_energy_j": 6.0e-6,
            "result_switching_frequency_hz": 5.0e4,
            "result_switching_loss_w": 0.1,
            "result_event_grid_sha256": "b" * 64,
            "result_switching_waveform_sha256": "c" * 64,
            "accepted_switching_loss_result_sha256": "d" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mosfet_switching_loss_uses_current_gate_charge_overlap_deadtime_events_temperature_cycle_and_result"
    ]


def test_v29_public_step_response_rise_settling_overshoot_initial_state_threshold_window_mismatch():
    summary = _v29()
    contract = summary["metrics"]["positive"][
        "step_response_initial_final_rise_threshold_settling_band_overshoot_window_waveform_result_generation_identity"
    ]
    contract.update(
        {
            "initial_step_generation_id": "step-response-160",
            "window_step_generation_id": "step-response-159",
            "result_step_generation_id": "step-response-158",
            "result_initial_value": 1.0,
            "result_final_value": 4.0,
            "result_rise_threshold_fractions": [0.2, 0.8],
            "result_rise_crossing_times_s": [4.0e-4, 1.0e-4],
            "result_rise_time_s": -3.0e-4,
            "result_settling_band_fraction": 0.05,
            "result_settling_time_s": 0.02,
            "result_overshoot_peak": 4.5,
            "result_overshoot_fraction": -0.1,
            "result_measurement_window_s": [0.002, 0.001],
            "result_waveform_sha256": "e" * 64,
            "accepted_step_result_sha256": "f" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_response_uses_current_initial_final_rise_settling_overshoot_window_and_waveform"
    ]
