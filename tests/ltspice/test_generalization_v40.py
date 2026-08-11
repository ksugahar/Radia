import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v39 import _v39


_PROMOTED_CASE_IDS = (
    "v40_public_llc_resonant_gain_frequency_zvs_current_loss_cycle_energy_mismatch",
    "v40_public_bjt_bias_gm_gain_pole_noise_distortion_thermal_power_mismatch",
)
_LLC_KEY = (
    "llc_resonant_elements_frequency_gain_zvs_current_loss_cycle_energy_"
    "waveform_result_identity"
)
_BJT_KEY = (
    "bjt_bias_gm_gain_pole_noise_distortion_thermal_power_circuit_result_identity"
)


def _v40():
    payload = _v39()
    positive = payload["metrics"]["positive"]
    generation = "llc-resonant-724"
    inductance, capacitance, magnetizing = 10.0e-6, 100.0e-9, 60.0e-6
    resonance = 1.0 / (2.0 * math.pi * math.sqrt(inductance * capacitance))
    frequency, voltage, turns, gain = 150_000.0, 400.0, 0.1, 0.95
    current, resistance = 4.0, 0.08
    conduction, switching, tank = current**2 * resistance, 0.2, 2.0
    device = conduction + switching
    dead_time, coss = 200.0e-9, 200.0e-12
    output_power = 380.0
    mirrored = {
        "resonant_inductance_h": inductance,
        "resonant_capacitance_f": capacitance,
        "magnetizing_inductance_h": magnetizing,
        "resonant_frequency_hz": resonance,
        "switching_frequency_hz": frequency,
        "normalized_switching_frequency": frequency / resonance,
        "input_voltage_v": voltage,
        "transformer_turns_ratio": turns,
        "conversion_gain_v_per_v": gain,
        "output_voltage_v": voltage * turns * gain,
        "dead_time_s": dead_time,
        "switch_output_capacitance_f": coss,
        "magnetizing_current_at_transition_a": 1.2,
        "required_zvs_current_a": 2.0 * coss * voltage / dead_time,
        "zvs_condition_met": True,
        "circulating_current_rms_a": current,
        "switch_conduction_resistance_ohm": resistance,
        "conduction_loss_w": conduction,
        "switching_loss_w": switching,
        "device_loss_w": device,
        "tank_loss_w": tank,
        "input_energy_per_cycle_j": (output_power + device + tank) / frequency,
        "output_energy_per_cycle_j": output_power / frequency,
        "device_loss_energy_per_cycle_j": device / frequency,
        "tank_loss_energy_per_cycle_j": tank / frequency,
        "cycle_energy_residual_j": 0.0,
    }
    positive[_LLC_KEY] = {
        "llc_generation_id": generation,
        **{key: generation for key in (
            "element_llc_generation_id", "frequency_llc_generation_id",
            "gain_llc_generation_id", "zvs_llc_generation_id",
            "current_llc_generation_id", "loss_llc_generation_id",
            "energy_llc_generation_id", "waveform_llc_generation_id",
            "result_llc_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "waveform_owner": "llc/waveform-724",
        "accepted_waveform_owner": "llc/waveform-724",
        "waveform_sha256": "1" * 64, "accepted_waveform_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }

    generation = "bjt-amplifier-724"
    current, thermal_voltage = 1.0e-3, 25.85e-3
    gm = current / thermal_voltage
    collector_resistance, load_resistance = 4_700.0, 10_000.0
    parallel = 1.0 / (1.0 / collector_resistance + 1.0 / load_resistance)
    gain = -gm * parallel
    capacitance = 20.0e-12
    pole = 1.0 / (2.0 * math.pi * parallel * capacitance)
    noise_density = 4.0e-9
    second, third = 0.01, 0.005
    power = 5.0 * current
    mirrored = {
        "collector_current_a": current, "collector_emitter_voltage_v": 5.0,
        "thermal_voltage_v": thermal_voltage, "transconductance_s": gm,
        "collector_resistance_ohm": collector_resistance,
        "load_resistance_ohm": load_resistance, "voltage_gain_v_per_v": gain,
        "output_capacitance_f": capacitance, "dominant_pole_frequency_hz": pole,
        "input_noise_density_v_per_sqrt_hz": noise_density,
        "integrated_output_noise_v_rms": noise_density * abs(gain) * math.sqrt(pole),
        "fundamental_output_v": 1.0, "second_harmonic_v": second,
        "third_harmonic_v": third,
        "total_harmonic_distortion": math.sqrt(second**2 + third**2),
        "ambient_temperature_c": 25.0, "junction_to_ambient_k_per_w": 100.0,
        "device_power_w": power, "junction_temperature_c": 25.0 + 100.0 * power,
    }
    positive[_BJT_KEY] = {
        "bjt_generation_id": generation,
        **{key: generation for key in (
            "bias_bjt_generation_id", "gm_bjt_generation_id", "gain_bjt_generation_id",
            "pole_bjt_generation_id", "noise_bjt_generation_id",
            "distortion_bjt_generation_id", "thermal_bjt_generation_id",
            "circuit_bjt_generation_id", "result_bjt_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "circuit_owner": "bjt/circuit-724",
        "accepted_circuit_owner": "bjt/circuit-724",
        "circuit_sha256": "3" * 64, "accepted_circuit_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v40_public_positive_llc_and_bjt_closure():
    assert ideal_transformer_identity_gate(_v40())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v40_public_llc_resonant_gain_frequency_zvs_current_loss_cycle_energy_mismatch():
    payload = _v40()
    identity = payload["metrics"]["positive"][_LLC_KEY]
    identity.update({
        "frequency_llc_generation_id": "llc-resonant-723",
        "result_resonant_frequency_hz": -1.0,
        "result_conversion_gain_v_per_v": 9.0,
        "result_zvs_condition_met": False,
        "result_device_loss_w": -1.0,
        "result_cycle_energy_residual_j": 9.0,
        "accepted_waveform_owner": "llc/old",
        "accepted_result_sha256": "a" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["llc_resonant_converters_use_current_elements_frequency_gain_zvs_current_loss_cycle_energy_waveform_and_result"]


def test_v40_public_bjt_bias_gm_gain_pole_noise_distortion_thermal_power_mismatch():
    payload = _v40()
    identity = payload["metrics"]["positive"][_BJT_KEY]
    identity.update({
        "bias_bjt_generation_id": "bjt-amplifier-723",
        "result_collector_current_a": -1.0,
        "result_transconductance_s": -1.0,
        "result_voltage_gain_v_per_v": 9.0,
        "result_dominant_pole_frequency_hz": -1.0,
        "result_total_harmonic_distortion": -1.0,
        "result_junction_temperature_c": -273.15,
        "accepted_circuit_owner": "bjt/old",
        "accepted_result_sha256": "b" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["bjt_amplifiers_use_current_bias_gm_gain_pole_noise_distortion_thermal_power_circuit_and_result"]


def test_v40_public_rejects_self_consistent_llc_energy_creation():
    payload = _v40()
    identity = payload["metrics"]["positive"][_LLC_KEY]
    identity["input_energy_per_cycle_j"] = identity["result_input_energy_per_cycle_j"] = 0.0
    identity["cycle_energy_residual_j"] = identity["result_cycle_energy_residual_j"] = 0.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_wrong_bjt_transconductance():
    payload = _v40()
    identity = payload["metrics"]["positive"][_BJT_KEY]
    identity["transconductance_s"] = identity["result_transconductance_s"] = 1.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
