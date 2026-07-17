from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v25 import _v25


_PROMOTED_CASE_IDS = (
    "v26_public_ac_sweep_mode_frequency_grid_complex_phase_unwrap_measure_generation_mismatch",
    "v26_public_electrothermal_device_power_temperature_model_thermal_network_timestep_mismatch",
)


def _v26():
    summary = _v25()
    positive = summary["metrics"]["positive"]
    generation = "ac-sweep-131"
    positive[
        "ac_sweep_mode_frequency_grid_complex_phase_unwrap_measure_generation_identity"
    ] = {
        "ac_generation_id": generation,
        "sweep_mode_ac_generation_id": generation,
        "frequency_grid_ac_generation_id": generation,
        "complex_basis_ac_generation_id": generation,
        "phase_unwrap_ac_generation_id": generation,
        "measure_row_ac_generation_id": generation,
        "result_ac_generation_id": generation,
        "sweep_mode": "decade",
        "result_sweep_mode": "decade",
        "points_per_decade": 20,
        "result_points_per_decade": 20,
        "frequency_grid_hz": [10.0, 31.6227766, 100.0, 316.227766, 1000.0],
        "result_frequency_grid_hz": [10.0, 31.6227766, 100.0, 316.227766, 1000.0],
        "complex_basis": "real_imaginary",
        "result_complex_basis": "real_imaginary",
        "phase_unwrap": "continuous_radians",
        "result_phase_unwrap": "continuous_radians",
        "measure_row_ids": ["gain_10", "gain_100", "gain_1000"],
        "result_measure_row_ids": ["gain_10", "gain_100", "gain_1000"],
        "measure_table_sha256": "1" * 64,
        "result_measure_table_sha256": "1" * 64,
    }
    generation = "electrothermal-131"
    positive[
        "electrothermal_device_power_temperature_model_thermal_network_timestep_generation_identity"
    ] = {
        "electrothermal_generation_id": generation,
        "device_power_electrothermal_generation_id": generation,
        "temperature_model_electrothermal_generation_id": generation,
        "thermal_network_electrothermal_generation_id": generation,
        "timestep_electrothermal_generation_id": generation,
        "result_electrothermal_generation_id": generation,
        "device_power_trace_ids": ["P(M1)", "P(D1)"],
        "result_device_power_trace_ids": ["P(M1)", "P(D1)"],
        "device_power_sha256": "2" * 64,
        "result_device_power_sha256": "2" * 64,
        "temperature_model_id": "linear_temperature_coefficients",
        "result_temperature_model_id": "linear_temperature_coefficients",
        "temperature_model_sha256": "3" * 64,
        "result_temperature_model_sha256": "3" * 64,
        "thermal_network_id": "foster_rc_3stage",
        "result_thermal_network_id": "foster_rc_3stage",
        "thermal_network_sha256": "4" * 64,
        "result_thermal_network_sha256": "4" * 64,
        "time_step_s": 1.0e-5,
        "result_time_step_s": 1.0e-5,
        "time_grid_s": [0.0, 1.0e-5, 2.0e-5, 3.0e-5],
        "result_time_grid_s": [0.0, 1.0e-5, 2.0e-5, 3.0e-5],
        "temperature_waveform_sha256": "5" * 64,
        "result_temperature_waveform_sha256": "5" * 64,
    }
    return summary


def test_v26_positive():
    assert ideal_transformer_identity_gate(_v26())["status"] == "ok"


def test_v26_public_ac_sweep_identity_mismatch():
    summary = _v26()
    contract = summary["metrics"]["positive"][
        "ac_sweep_mode_frequency_grid_complex_phase_unwrap_measure_generation_identity"
    ]
    contract.update(
        {
            "sweep_mode_ac_generation_id": "ac-sweep-130",
            "result_sweep_mode": "linear",
            "result_frequency_grid_hz": [10.0, 100.0, 1000.0],
            "result_complex_basis": "magnitude_phase_degrees",
            "result_phase_unwrap": "wrapped_degrees",
            "result_measure_row_ids": ["gain_old"],
            "result_measure_table_sha256": "8" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "ac_measures_use_current_sweep_grid_complex_basis_phase_unwrap_rows_and_result"
    ]


def test_v26_public_electrothermal_identity_mismatch():
    summary = _v26()
    contract = summary["metrics"]["positive"][
        "electrothermal_device_power_temperature_model_thermal_network_timestep_generation_identity"
    ]
    contract.update(
        {
            "device_power_electrothermal_generation_id": "electrothermal-130",
            "result_device_power_trace_ids": ["P(M0)"],
            "result_device_power_sha256": "9" * 64,
            "result_temperature_model_id": "constant_ambient",
            "result_thermal_network_id": "cauer_rc_previous",
            "result_time_step_s": 2.0e-5,
            "result_time_grid_s": [0.0, 2.0e-5, 4.0e-5],
            "result_temperature_waveform_sha256": "c" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "electrothermal_waveforms_use_current_power_temperature_model_network_timestep_and_result"
    ]
