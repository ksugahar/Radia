import math
from copy import deepcopy

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v40 import _v40


_PROMOTED_CASE_IDS = (
    "v41_public_boost_duty_inductorripple_ccmboundary_efficiency_switchstress_energy_mismatch",
    "v41_public_activefilter_polezero_q_gain_noise_slew_saturation_power_mismatch",
)
_BOOST_KEY = (
    "boost_duty_inductor_ripple_ccm_output_ripple_efficiency_stress_cycle_"
    "energy_waveform_result_identity"
)
_FILTER_KEY = (
    "active_filter_pole_zero_q_gain_noise_slew_saturation_power_circuit_"
    "result_identity"
)


def _v41():
    payload = deepcopy(_v40())
    positive = payload["metrics"]["positive"]
    generation = "boost-731"
    vin, vout, frequency = 12.0, 24.0, 100_000.0
    inductance, capacitance, output_current = 100.0e-6, 100.0e-6, 1.0
    duty, efficiency = 1.0 - vin / vout, 0.90
    ripple = vin * duty / (inductance * frequency)
    output_power = vout * output_current
    input_power = output_power / efficiency
    average_current = input_power / vin
    mirrored = {
        "input_voltage_v": vin, "output_voltage_v": vout, "duty_ratio": duty,
        "switching_frequency_hz": frequency, "inductance_h": inductance,
        "inductor_average_current_a": average_current,
        "inductor_ripple_peak_to_peak_a": ripple,
        "ccm_boundary_current_a": ripple / 2.0,
        "ccm_condition_met": True,
        "output_capacitance_f": capacitance, "output_current_a": output_current,
        "output_voltage_ripple_peak_to_peak_v": output_current * duty / (capacitance * frequency),
        "efficiency": efficiency, "switch_voltage_stress_v": vout,
        "diode_reverse_voltage_stress_v": vout,
        "switch_peak_current_a": average_current + ripple / 2.0,
        "input_energy_per_cycle_j": input_power / frequency,
        "output_energy_per_cycle_j": output_power / frequency,
        "loss_energy_per_cycle_j": (input_power - output_power) / frequency,
        "cycle_energy_residual_j": 0.0,
    }
    positive[_BOOST_KEY] = {
        "boost_generation_id": generation,
        **{key: generation for key in (
            "duty_boost_generation_id", "ripple_boost_generation_id",
            "ccm_boost_generation_id", "output_boost_generation_id",
            "efficiency_boost_generation_id", "stress_boost_generation_id",
            "energy_boost_generation_id", "waveform_boost_generation_id",
            "result_boost_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "waveform_owner": "boost/waveform-731",
        "accepted_waveform_owner": "boost/waveform-731",
        "waveform_sha256": "1" * 64, "accepted_waveform_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }

    generation = "active-filter-731"
    natural_frequency = 2.0 * math.pi * 1_000.0
    quality = 1.0 / math.sqrt(2.0)
    pole_real = -natural_frequency / (2.0 * quality)
    pole_imag = natural_frequency * math.sqrt(1.0 - 1.0 / (4.0 * quality**2))
    noise_density, noise_bandwidth = 10.0e-9, 2_500.0
    mirrored = {
        "pole_locations_rad_s": [[pole_real, pole_imag], [pole_real, -pole_imag]],
        "zero_locations_rad_s": [], "natural_frequency_rad_s": natural_frequency,
        "quality_factor": quality, "passband_gain_v_per_v": 2.0,
        "input_noise_density_v_per_sqrt_hz": noise_density,
        "noise_bandwidth_hz": noise_bandwidth,
        "integrated_output_noise_v_rms": noise_density * math.sqrt(noise_bandwidth),
        "input_frequency_hz": 1_000.0, "input_amplitude_v_peak": 0.5,
        "output_amplitude_v_peak": 1.0,
        "slew_rate_demand_v_per_s": 2.0 * math.pi * 1_000.0,
        "available_slew_rate_v_per_s": 5.0e6,
        "positive_supply_v": 12.0, "negative_supply_v": -12.0,
        "saturation_margin_v": 11.0, "quiescent_supply_current_a": 0.002,
        "supply_power_w": 0.048,
    }
    positive[_FILTER_KEY] = {
        "filter_generation_id": generation,
        **{key: generation for key in (
            "pole_filter_generation_id", "zero_filter_generation_id",
            "q_filter_generation_id", "gain_filter_generation_id",
            "noise_filter_generation_id", "slew_filter_generation_id",
            "saturation_filter_generation_id", "power_filter_generation_id",
            "circuit_filter_generation_id", "result_filter_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "circuit_owner": "active-filter/circuit-731",
        "accepted_circuit_owner": "active-filter/circuit-731",
        "circuit_sha256": "3" * 64, "accepted_circuit_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v41_public_positive_boost_and_active_filter_closure():
    assert ideal_transformer_identity_gate(_v41())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v41_public_boost_duty_inductorripple_ccmboundary_efficiency_switchstress_energy_mismatch():
    payload = _v41()
    identity = payload["metrics"]["positive"][_BOOST_KEY]
    identity.update({
        "duty_boost_generation_id": "boost-730",
        "result_duty_ratio": 0.9, "result_inductor_ripple_peak_to_peak_a": -1.0,
        "result_ccm_condition_met": False, "result_efficiency": 1.5,
        "result_switch_voltage_stress_v": -1.0,
        "result_cycle_energy_residual_j": 1.0,
        "accepted_waveform_owner": "boost/old", "accepted_result_sha256": "a" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "boost_converters_use_current_duty_ripple_ccm_output_efficiency_stress_energy_waveform_and_result"
    ]


def test_v41_public_activefilter_polezero_q_gain_noise_slew_saturation_power_mismatch():
    payload = _v41()
    identity = payload["metrics"]["positive"][_FILTER_KEY]
    identity.update({
        "pole_filter_generation_id": "active-filter-730",
        "result_pole_locations_rad_s": [[1.0, 0.0]],
        "result_quality_factor": -1.0, "result_passband_gain_v_per_v": -1.0,
        "result_integrated_output_noise_v_rms": -1.0,
        "result_slew_rate_demand_v_per_s": 1.0e9,
        "result_saturation_margin_v": -1.0, "result_supply_power_w": -1.0,
        "accepted_circuit_owner": "active-filter/old",
        "accepted_result_sha256": "b" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "active_filters_use_current_poles_zeros_q_gain_noise_slew_saturation_power_circuit_and_result"
    ]


def test_v41_public_rejects_self_consistent_boost_energy_creation():
    payload = _v41()
    identity = payload["metrics"]["positive"][_BOOST_KEY]
    identity["input_energy_per_cycle_j"] = identity["result_input_energy_per_cycle_j"] = 0.0
    identity["cycle_energy_residual_j"] = identity["result_cycle_energy_residual_j"] = 0.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_unstable_filter_poles():
    payload = _v41()
    identity = payload["metrics"]["positive"][_FILTER_KEY]
    poles = [[1.0, 2.0], [1.0, -2.0]]
    identity["pole_locations_rad_s"] = identity["result_pole_locations_rad_s"] = poles
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
