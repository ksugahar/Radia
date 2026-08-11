import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v37 import _v37


_PROMOTED_CASE_IDS = (
    "v38_public_sallenkey_poles_q_gain_opamp_gbw_phase_noise_step_owner_mismatch",
    "v38_public_bridge_rectifier_inrush_ripple_charge_conduction_loss_load_cycle_energy_mismatch",
)


def _v38():
    payload = _v37()
    positive = payload["metrics"]["positive"]
    generation = "sallen-key-341"
    resistance = 10_000.0
    capacitance = 15.915494309189533e-9
    mirrored = {
        "r1_ohm": resistance, "r2_ohm": resistance,
        "c1_f": capacitance, "c2_f": capacitance,
        "noninverting_gain": 1.0,
        "pole_frequency_hz": 1.0 / (2.0 * math.pi * resistance * capacitance),
        "quality_factor": 0.5, "dc_gain": 1.0,
        "opamp_gbw_hz": 1.0e6, "phase_margin_deg": 60.0,
        "integrated_output_noise_v_rms": 20.0e-6,
        "step_overshoot_fraction": 0.0, "step_settling_time_s": 5.0e-3,
    }
    positive["sallenkey_poles_q_gain_opamp_gbw_phase_noise_step_circuit_owner_result_identity"] = {
        "filter_generation_id": generation,
        **{key: generation for key in (
            "pole_filter_generation_id", "q_filter_generation_id",
            "gain_filter_generation_id", "opamp_filter_generation_id",
            "noise_filter_generation_id", "step_filter_generation_id",
            "circuit_filter_generation_id", "result_filter_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "circuit_owner": "filter/sallen-key-341",
        "accepted_circuit_owner": "filter/sallen-key-341",
        "circuit_sha256": "1" * 64, "accepted_circuit_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }

    generation = "bridge-rectifier-341"
    line_frequency, capacitance = 60.0, 1.0e-3
    ripple_frequency, voltage, resistance = 120.0, 12.0, 1.0e3
    current = voltage / resistance
    load_power = voltage * current
    diode_loss = 2.0 * 0.7 * current
    period = 1.0 / ripple_frequency
    mirrored = {
        "line_frequency_hz": line_frequency, "ripple_frequency_hz": ripple_frequency,
        "capacitor_f": capacitance, "average_output_voltage_v": voltage,
        "load_resistance_ohm": resistance, "average_load_current_a": current,
        "ripple_voltage_pp_v": current / (ripple_frequency * capacitance),
        "inrush_peak_current_a": 10.0,
        "capacitor_charge_per_cycle_c": current / ripple_frequency,
        "rectifier_delivered_charge_per_cycle_c": current / ripple_frequency,
        "diode_conduction_fraction": 0.08, "diode_forward_drop_v": 0.7,
        "average_diode_loss_w": diode_loss, "average_load_power_w": load_power,
        "ripple_cycle_period_s": period,
        "source_energy_per_cycle_j": (load_power + diode_loss) * period,
        "load_energy_per_cycle_j": load_power * period,
        "diode_energy_per_cycle_j": diode_loss * period,
    }
    positive["bridge_rectifier_inrush_ripple_charge_conduction_loss_load_cycle_energy_waveform_result_identity"] = {
        "rectifier_generation_id": generation,
        **{key: generation for key in (
            "inrush_rectifier_generation_id", "ripple_rectifier_generation_id",
            "charge_rectifier_generation_id", "conduction_rectifier_generation_id",
            "loss_rectifier_generation_id", "load_rectifier_generation_id",
            "energy_rectifier_generation_id", "waveform_rectifier_generation_id",
            "result_rectifier_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "waveform_owner": "rectifier/waveform-341",
        "accepted_waveform_owner": "rectifier/waveform-341",
        "waveform_sha256": "3" * 64, "accepted_waveform_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v38_public_positive_sallen_key_and_bridge_rectifier_closure():
    assert ideal_transformer_identity_gate(_v38())["status"] == "ok"


def test_v38_public_sallenkey_poles_q_gain_opamp_gbw_phase_noise_step_owner_mismatch():
    payload = _v38()
    identity = payload["metrics"]["positive"]["sallenkey_poles_q_gain_opamp_gbw_phase_noise_step_circuit_owner_result_identity"]
    identity.update({
        "pole_filter_generation_id": "sallen-key-340",
        "result_filter_generation_id": "sallen-key-339",
        "result_pole_frequency_hz": 2000.0, "result_quality_factor": 5.0,
        "result_dc_gain": -1.0, "result_opamp_gbw_hz": 100.0,
        "result_phase_margin_deg": -10.0,
        "result_integrated_output_noise_v_rms": -1.0,
        "result_step_overshoot_fraction": 2.0,
        "result_step_settling_time_s": -1.0,
        "accepted_circuit_owner": "filter/old", "accepted_result_sha256": "9" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["sallen_key_uses_current_poles_q_gain_opamp_gbw_phase_noise_step_circuit_owner_and_result"]


def test_v38_public_bridge_rectifier_inrush_ripple_charge_conduction_loss_load_cycle_energy_mismatch():
    payload = _v38()
    identity = payload["metrics"]["positive"]["bridge_rectifier_inrush_ripple_charge_conduction_loss_load_cycle_energy_waveform_result_identity"]
    identity.update({
        "charge_rectifier_generation_id": "bridge-rectifier-340",
        "result_rectifier_generation_id": "bridge-rectifier-339",
        "result_inrush_peak_current_a": -1.0, "result_ripple_frequency_hz": 60.0,
        "result_ripple_voltage_pp_v": 2.0,
        "result_capacitor_charge_per_cycle_c": -1.0,
        "result_rectifier_delivered_charge_per_cycle_c": 2.0,
        "result_diode_conduction_fraction": 2.0,
        "result_average_diode_loss_w": -1.0,
        "result_average_load_power_w": -1.0,
        "result_source_energy_per_cycle_j": 0.0,
        "accepted_waveform_owner": "rectifier/old",
        "accepted_result_sha256": "a" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["bridge_rectifier_uses_current_inrush_ripple_charge_conduction_loss_load_cycle_energy_waveform_and_result"]


def test_v38_public_rejects_self_consistent_wrong_sallen_key_q():
    payload = _v38()
    identity = payload["metrics"]["positive"]["sallenkey_poles_q_gain_opamp_gbw_phase_noise_step_circuit_owner_result_identity"]
    identity["quality_factor"] = identity["result_quality_factor"] = 2.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_rectifier_energy_creation():
    payload = _v38()
    identity = payload["metrics"]["positive"]["bridge_rectifier_inrush_ripple_charge_conduction_loss_load_cycle_energy_waveform_result_identity"]
    identity["source_energy_per_cycle_j"] = identity["result_source_energy_per_cycle_j"] = 0.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
