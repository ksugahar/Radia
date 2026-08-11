from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v27 import _v27


_PROMOTED_CASE_IDS = (
    "v28_public_smps_efficiency_energy_integration_steady_cycle_window_source_load_waveform_mismatch",
    "v28_public_loop_gain_break_injection_sign_phase_unwrap_crossover_margin_frequency_grid_mismatch",
)


def _v28():
    summary = _v27()
    positive = summary["metrics"]["positive"]
    generation = "smps-efficiency-151"
    positive[
        "smps_efficiency_source_load_steady_cycle_energy_integration_switching_waveform_timestep_result_generation_identity"
    ] = {
        "efficiency_generation_id": generation,
        "source_trace_efficiency_generation_id": generation,
        "load_trace_efficiency_generation_id": generation,
        "cycle_window_efficiency_generation_id": generation,
        "integration_efficiency_generation_id": generation,
        "waveform_efficiency_generation_id": generation,
        "timestep_efficiency_generation_id": generation,
        "result_efficiency_generation_id": generation,
        "source_trace_id": "V(in)*-I(VIN)",
        "result_source_trace_id": "V(in)*-I(VIN)",
        "load_trace_id": "V(out)*I(RLOAD)",
        "result_load_trace_id": "V(out)*I(RLOAD)",
        "steady_cycle_window_s": [0.004, 0.005],
        "result_steady_cycle_window_s": [0.004, 0.005],
        "energy_integration_rule": "trapezoid_power_over_time",
        "result_energy_integration_rule": "trapezoid_power_over_time",
        "source_energy_j": 0.010,
        "result_source_energy_j": 0.010,
        "load_energy_j": 0.009,
        "result_load_energy_j": 0.009,
        "efficiency": 0.9,
        "result_efficiency": 0.9,
        "max_timestep_s": 1.0e-7,
        "result_max_timestep_s": 1.0e-7,
        "switching_waveform_sha256": "1" * 64,
        "result_switching_waveform_sha256": "1" * 64,
        "efficiency_result_sha256": "2" * 64,
        "accepted_efficiency_result_sha256": "2" * 64,
    }
    generation = "loop-gain-151"
    positive[
        "loop_gain_break_injection_sign_phase_unwrap_crossover_margin_frequency_grid_result_generation_identity"
    ] = {
        "loop_gain_generation_id": generation,
        "break_loop_gain_generation_id": generation,
        "injection_loop_gain_generation_id": generation,
        "phase_loop_gain_generation_id": generation,
        "crossover_loop_gain_generation_id": generation,
        "frequency_loop_gain_generation_id": generation,
        "result_loop_gain_generation_id": generation,
        "loop_break_element": "VLOOP",
        "result_loop_break_element": "VLOOP",
        "loop_break_nodes": ["fb", "comp"],
        "result_loop_break_nodes": ["fb", "comp"],
        "injection_sign": 1,
        "result_injection_sign": 1,
        "phase_unwrap_rule": "continuous_negative_180",
        "result_phase_unwrap_rule": "continuous_negative_180",
        "crossover_interpolation": "log_frequency_linear_db",
        "result_crossover_interpolation": "log_frequency_linear_db",
        "frequency_grid_hz": [10.0, 100.0, 1000.0, 10000.0, 100000.0],
        "result_frequency_grid_hz": [10.0, 100.0, 1000.0, 10000.0, 100000.0],
        "loop_gain_db": [40.0, 20.0, 0.0, -20.0, -40.0],
        "result_loop_gain_db": [40.0, 20.0, 0.0, -20.0, -40.0],
        "phase_deg": [-90.0, -120.0, -150.0, -190.0, -240.0],
        "result_phase_deg": [-90.0, -120.0, -150.0, -190.0, -240.0],
        "gain_crossover_hz": 1000.0,
        "result_gain_crossover_hz": 1000.0,
        "phase_margin_deg": 30.0,
        "result_phase_margin_deg": 30.0,
        "gain_margin_db": 10.0,
        "result_gain_margin_db": 10.0,
        "loop_gain_result_sha256": "3" * 64,
        "accepted_loop_gain_result_sha256": "3" * 64,
    }
    return summary


def test_v28_positive_smps_efficiency_and_loop_gain_identities():
    assert ideal_transformer_identity_gate(_v28())["status"] == "ok"


def test_v28_public_smps_efficiency_energy_integration_steady_cycle_window_source_load_waveform_mismatch():
    summary = _v28()
    contract = summary["metrics"]["positive"][
        "smps_efficiency_source_load_steady_cycle_energy_integration_switching_waveform_timestep_result_generation_identity"
    ]
    contract.update(
        {
            "source_trace_efficiency_generation_id": "smps-efficiency-150",
            "cycle_window_efficiency_generation_id": "smps-efficiency-149",
            "result_efficiency_generation_id": "smps-efficiency-148",
            "result_source_trace_id": "V(in)*I(VIN)",
            "result_load_trace_id": "V(out)*-I(RLOAD)",
            "result_steady_cycle_window_s": [0.0, 0.0002],
            "result_energy_integration_rule": "sample_average_ratio",
            "result_source_energy_j": 0.020,
            "result_load_energy_j": 0.012,
            "result_efficiency": 1.2,
            "result_max_timestep_s": 1.0e-4,
            "result_switching_waveform_sha256": "9" * 64,
            "accepted_efficiency_result_sha256": "a" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "smps_efficiency_uses_current_traces_cycle_window_integration_waveform_timestep_and_result"
    ]


def test_v28_public_loop_gain_break_injection_sign_phase_unwrap_crossover_margin_frequency_grid_mismatch():
    summary = _v28()
    contract = summary["metrics"]["positive"][
        "loop_gain_break_injection_sign_phase_unwrap_crossover_margin_frequency_grid_result_generation_identity"
    ]
    contract.update(
        {
            "break_loop_gain_generation_id": "loop-gain-150",
            "phase_loop_gain_generation_id": "loop-gain-149",
            "result_loop_gain_generation_id": "loop-gain-148",
            "result_loop_break_element": "VOLD",
            "result_loop_break_nodes": ["out", "0"],
            "result_injection_sign": -1,
            "result_phase_unwrap_rule": "principal_value",
            "result_crossover_interpolation": "linear_frequency_nearest",
            "result_frequency_grid_hz": [10.0, 500.0, 50000.0],
            "result_loop_gain_db": [20.0, -10.0],
            "result_phase_deg": [90.0, 120.0],
            "result_gain_crossover_hz": 500.0,
            "result_phase_margin_deg": 270.0,
            "result_gain_margin_db": -5.0,
            "accepted_loop_gain_result_sha256": "b" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "loop_gain_uses_current_break_injection_phase_grid_crossover_margins_and_result"
    ]
