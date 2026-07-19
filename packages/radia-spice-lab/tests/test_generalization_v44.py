from copy import deepcopy
from math import sqrt

from ltspice_converter.ltspice_v44_gates import validate_ltspice_v44_public_identity


BUCK_CASE = "v44_public_buck_transient_startup_inductorcurrent_outputripple_efficiency_energy_waveform_mismatch"
NOISE_CASE = "v44_public_noise_ac_transfer_inputreferred_psd_bandwidth_correlation_measure_owner_mismatch"


def _contract(generation, fields, owner_key, owner):
    return {
        "buck_generation_id": generation,
        **{name: generation for name in ("startup_buck_generation_id", "inductor_current_buck_generation_id", "ripple_buck_generation_id", "efficiency_buck_generation_id", "energy_buck_generation_id", "waveform_buck_generation_id", "result_buck_generation_id")},
        **fields, **{f"result_{name}": value for name, value in fields.items()},
        owner_key: owner, f"accepted_{owner_key}": owner,
        "waveform_sha256": "1" * 64, "accepted_waveform_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }


def _positive():
    buck = _contract("buck-transient-844", {
        "startup_window_start_s": 1e-4, "startup_window_stop_s": 3e-4,
        "inductor_current_avg_a": 2.0, "inductor_current_ripple_a": 0.4,
        "output_voltage_v": 5.0, "output_ripple_v": 0.05,
        "input_power_w": 10.0, "output_power_w": 9.0, "efficiency_fraction": 0.9,
        "energy_start_j": 1e-6, "energy_end_j": 1.5e-6, "energy_balance_residual_j": 0.0,
    }, "waveform_owner", "buck/waveform-844")
    noise_fields = {
        "frequency_hz": 1000.0, "transfer_magnitude_v_per_v": 2.0,
        "input_referred_psd_v2_per_hz": 4e-12, "bandwidth_hz": 1000.0,
        "integrated_noise_v_rms": sqrt(4e-9), "correlation_coefficient": 0.25,
        "output_noise_v_rms": 2.0 * sqrt(4e-9),
    }
    noise = {
        "noise_generation_id": "noise-ac-844",
        **{name: "noise-ac-844" for name in ("transfer_noise_generation_id", "input_psd_noise_generation_id", "bandwidth_noise_generation_id", "correlation_noise_generation_id", "measure_noise_generation_id", "waveform_noise_generation_id", "result_noise_generation_id")},
        **noise_fields, **{f"result_{name}": value for name, value in noise_fields.items()},
        "measure_owner": "noise/measure-844", "accepted_measure_owner": "noise/measure-844",
        "waveform_sha256": "3" * 64, "accepted_waveform_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return {"buck_transient_startup_inductorcurrent_outputripple_efficiency_energy_waveform_result_identity": buck, "noise_ac_transfer_inputreferred_psd_bandwidth_correlation_measure_owner_identity": noise}


def test_v44_public_positive_transient_and_noise_closure():
    checks = validate_ltspice_v44_public_identity(_positive())
    assert checks == {"buck_v44_transient_identity": True, "noise_v44_ac_identity": True}


def test_v44_public_rejects_buck_result_mismatch():
    positive = _positive()
    positive["buck_transient_startup_inductorcurrent_outputripple_efficiency_energy_waveform_result_identity"]["result_output_power_w"] = 8.0
    assert validate_ltspice_v44_public_identity(positive)["buck_v44_transient_identity"] is False


def test_v44_public_rejects_noise_bandwidth_mismatch():
    positive = _positive()
    positive["noise_ac_transfer_inputreferred_psd_bandwidth_correlation_measure_owner_identity"]["result_bandwidth_hz"] = 900.0
    assert validate_ltspice_v44_public_identity(positive)["noise_v44_ac_identity"] is False


assert BUCK_CASE and NOISE_CASE
