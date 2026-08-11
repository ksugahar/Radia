import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v38 import _v38


_PROMOTED_CASE_IDS = (
    "v39_public_flyback_magnetizing_leakage_snubber_demag_flux_current_energy_mismatch",
    "v39_public_transimpedance_photodiode_capacitance_noise_gain_stability_step_mismatch",
)
_FLYBACK_KEY = (
    "flyback_magnetizing_leakage_snubber_demag_flux_current_cycle_energy_"
    "waveform_result_identity"
)
_TIA_KEY = (
    "transimpedance_photodiode_capacitance_noise_gain_bandwidth_stability_"
    "noise_step_circuit_result_identity"
)


def _v39():
    payload = _v38()
    positive = payload["metrics"]["positive"]
    generation = "flyback-719"
    frequency, magnetizing, leakage = 100_000.0, 100.0e-6, 2.0e-6
    period, voltage, duty = 1.0 / frequency, 24.0, 0.4
    peak_current = voltage * duty * period / magnetizing
    magnetizing_energy = 0.5 * magnetizing * peak_current**2
    leakage_energy = 0.5 * leakage * peak_current**2
    mirrored = {
        "switching_frequency_hz": frequency, "switching_period_s": period,
        "input_voltage_v": voltage, "duty_cycle": duty,
        "magnetizing_inductance_h": magnetizing, "leakage_inductance_h": leakage,
        "peak_switch_current_a": peak_current, "demagnetization_interval_s": 3.0e-6,
        "peak_core_flux_density_t": 0.2,
        "snubber_energy_per_cycle_j": leakage_energy,
        "average_snubber_dissipation_w": leakage_energy * frequency,
        "magnetizing_energy_per_cycle_j": magnetizing_energy,
        "delivered_energy_per_cycle_j": magnetizing_energy,
        "source_energy_per_cycle_j": magnetizing_energy + leakage_energy,
        "cycle_energy_residual_j": 0.0,
    }
    positive[_FLYBACK_KEY] = {
        "flyback_generation_id": generation,
        **{key: generation for key in (
            "magnetizing_flyback_generation_id", "leakage_flyback_generation_id",
            "snubber_flyback_generation_id", "demag_flyback_generation_id",
            "flux_flyback_generation_id", "current_flyback_generation_id",
            "energy_flyback_generation_id", "waveform_flyback_generation_id",
            "result_flyback_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "waveform_owner": "flyback/waveform-719",
        "accepted_waveform_owner": "flyback/waveform-719",
        "waveform_sha256": "1" * 64, "accepted_waveform_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }

    generation = "transimpedance-719"
    resistance, feedback_capacitance = 100_000.0, 2.0e-12
    detector_capacitance, input_capacitance = 5.0e-12, 2.0e-12
    bandwidth = 1.0 / (2.0 * math.pi * resistance * feedback_capacitance)
    noise_gain = 1.0 + (detector_capacitance + input_capacitance) / feedback_capacitance
    mirrored = {
        "feedback_resistance_ohm": resistance,
        "feedback_capacitance_f": feedback_capacitance,
        "photodiode_capacitance_f": detector_capacitance,
        "opamp_input_capacitance_f": input_capacitance,
        "transimpedance_dc_gain_v_per_a": resistance,
        "high_frequency_noise_gain": noise_gain,
        "closed_loop_bandwidth_hz": bandwidth,
        "opamp_gain_bandwidth_hz": 20.0e6, "phase_margin_deg": 60.0,
        "input_referred_current_noise_a_per_sqrt_hz": 2.0e-12,
        "integrated_output_noise_v_rms": 50.0e-6,
        "step_rise_time_s": 0.44 / bandwidth, "step_overshoot_fraction": 0.1,
        "step_settling_time_s": 5.0 / bandwidth,
    }
    positive[_TIA_KEY] = {
        "transimpedance_generation_id": generation,
        **{key: generation for key in (
            "capacitance_transimpedance_generation_id",
            "noisegain_transimpedance_generation_id",
            "bandwidth_transimpedance_generation_id",
            "stability_transimpedance_generation_id",
            "noise_transimpedance_generation_id", "step_transimpedance_generation_id",
            "circuit_transimpedance_generation_id", "result_transimpedance_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "circuit_owner": "transimpedance/circuit-719",
        "accepted_circuit_owner": "transimpedance/circuit-719",
        "circuit_sha256": "3" * 64, "accepted_circuit_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v39_public_positive_flyback_and_transimpedance_closure():
    assert ideal_transformer_identity_gate(_v39())["status"] == "ok"


def test_v39_public_flyback_magnetizing_leakage_snubber_demag_flux_current_energy_mismatch():
    payload = _v39()
    identity = payload["metrics"]["positive"][_FLYBACK_KEY]
    identity.update({
        "leakage_flyback_generation_id": "flyback-718",
        "energy_flyback_generation_id": "flyback-717",
        "result_flyback_generation_id": "flyback-716",
        "result_magnetizing_inductance_h": -1.0,
        "result_leakage_inductance_h": 9.0,
        "result_average_snubber_dissipation_w": -1.0,
        "result_demagnetization_interval_s": -1.0,
        "result_peak_core_flux_density_t": 9.0,
        "result_peak_switch_current_a": -1.0,
        "result_source_energy_per_cycle_j": 0.0,
        "accepted_waveform_owner": "flyback/old",
        "accepted_result_sha256": "a" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["flyback_uses_current_magnetizing_leakage_snubber_demag_flux_current_cycle_energy_waveform_and_result"]


def test_v39_public_transimpedance_photodiode_capacitance_noise_gain_stability_step_mismatch():
    payload = _v39()
    identity = payload["metrics"]["positive"][_TIA_KEY]
    identity.update({
        "capacitance_transimpedance_generation_id": "transimpedance-718",
        "result_transimpedance_generation_id": "transimpedance-717",
        "result_photodiode_capacitance_f": -1.0,
        "result_high_frequency_noise_gain": -1.0,
        "result_closed_loop_bandwidth_hz": -1.0,
        "result_phase_margin_deg": -10.0,
        "result_input_referred_current_noise_a_per_sqrt_hz": -1.0,
        "result_integrated_output_noise_v_rms": -1.0,
        "result_step_rise_time_s": -1.0,
        "accepted_circuit_owner": "transimpedance/old",
        "accepted_result_sha256": "b" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["transimpedance_uses_current_capacitance_noise_gain_bandwidth_stability_noise_step_circuit_and_result"]


def test_v39_public_rejects_self_consistent_flyback_energy_creation():
    payload = _v39()
    identity = payload["metrics"]["positive"][_FLYBACK_KEY]
    identity["source_energy_per_cycle_j"] = identity["result_source_energy_per_cycle_j"] = 0.0
    identity["cycle_energy_residual_j"] = identity["result_cycle_energy_residual_j"] = 0.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_wrong_tia_noise_gain():
    payload = _v39()
    identity = payload["metrics"]["positive"][_TIA_KEY]
    identity["high_frequency_noise_gain"] = identity["result_high_frequency_noise_gain"] = 1.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
