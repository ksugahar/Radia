from copy import deepcopy
import math

from ltspice_converter.ltspice_v55_gates import NOISE, POWER, validate_ltspice_v55_identity


PROMOTED_CASE_IDS = {
    "v55_public_acnoise_input_output_density_transfergain_bandwidth_owner_mismatch",
    "v55_public_switch_averagepower_input_output_loss_efficiency_window_owner_mismatch",
}


def _positive() -> dict[str, object]:
    noise_generation = "noise-public-v55"
    power_generation = "power-public-v55"
    frequency = [100.0, 1000.0, 10000.0]
    input_density = [1.0e-9, 2.0e-9, 1.5e-9]
    gain = [10.0, 8.0, 5.0]
    output_density = [source * transfer for source, transfer in zip(input_density, gain)]
    variance = sum(0.5 * (output_density[index] ** 2 + output_density[index + 1] ** 2) * (frequency[index + 1] - frequency[index]) for index in range(2))
    window = {"start_s": 0.002, "stop_s": 0.004}
    return {
        NOISE: {
            "generation_id": noise_generation, **{field: noise_generation for field in ("density_generation_id", "gain_generation_id", "bandwidth_generation_id", "grid_generation_id", "owner_generation_id", "result_generation_id")},
            "frequency_hz": frequency, "result_frequency_hz": frequency,
            "input_noise_density_v_per_sqrt_hz": input_density, "result_input_noise_density_v_per_sqrt_hz": input_density,
            "transfer_gain_v_per_v": gain, "result_transfer_gain_v_per_v": gain,
            "output_noise_density_v_per_sqrt_hz": output_density, "result_output_noise_density_v_per_sqrt_hz": output_density,
            "integration_bandwidth_hz": 9900.0, "result_integration_bandwidth_hz": 9900.0,
            "integrated_output_noise_v_rms": math.sqrt(variance), "result_integrated_output_noise_v_rms": math.sqrt(variance),
            "trace_owner": "trace:noise-v55", "result_trace_owner": "trace:noise-v55",
            "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
        },
        POWER: {
            "generation_id": power_generation, **{field: power_generation for field in ("input_generation_id", "output_generation_id", "loss_generation_id", "efficiency_generation_id", "window_generation_id", "owner_generation_id", "result_generation_id")},
            "average_input_power_w": 12.0, "result_average_input_power_w": 12.0,
            "average_output_power_w": 9.0, "result_average_output_power_w": 9.0,
            "average_loss_power_w": 3.0, "result_average_loss_power_w": 3.0,
            "efficiency": 0.75, "result_efficiency": 0.75,
            "measurement_window": window, "result_measurement_window": window,
            "waveform_owner": "waveform:power-v55", "result_waveform_owner": "waveform:power-v55",
            "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64,
        },
    }


def test_v55_positive_identity_is_accepted() -> None:
    assert validate_ltspice_v55_identity(_positive()) is True


def test_v55_frozen_public_counterfactuals_are_rejected() -> None:
    value = deepcopy(_positive())
    value[NOISE].update({"result_frequency_hz": [100.0], "result_trace_owner": "trace:stale"})
    value[POWER].update({"result_average_output_power_w": 11.0, "result_waveform_owner": "waveform:stale"})
    assert validate_ltspice_v55_identity(value) is False


def test_v55_self_consistent_noise_and_power_contradictions_are_rejected() -> None:
    value = deepcopy(_positive())
    bad_output = [1.0e-9, 1.0e-9, 1.0e-9]
    value[NOISE]["output_noise_density_v_per_sqrt_hz"] = value[NOISE]["result_output_noise_density_v_per_sqrt_hz"] = bad_output
    value[POWER]["average_loss_power_w"] = value[POWER]["result_average_loss_power_w"] = 1.0
    assert validate_ltspice_v55_identity(value) is False


def test_v55_malformed_values_reject_without_raising() -> None:
    value = deepcopy(_positive())
    value[NOISE]["frequency_hz"] = [[100.0]]
    value[POWER]["measurement_window"] = {"start_s": [0.0], "stop_s": 1.0}
    assert validate_ltspice_v55_identity(value) is False
