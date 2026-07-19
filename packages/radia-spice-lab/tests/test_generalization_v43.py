from copy import deepcopy

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v42 import _v42


_PROMOTED_CASE_IDS = (
    "v43_public_flyback_turnsratio_magnetizing_leakage_clamp_stress_loss_power_energy_mismatch",
    "v43_public_opamp_loopgain_phasemargin_crossover_step_overshoot_slew_power_mismatch",
)


def _v43():
    payload = deepcopy(_v42())
    positive = payload["metrics"]["positive"]
    generation = "flyback-843"
    values = {
        "primary_turns": 40.0, "secondary_turns": 10.0, "turns_ratio": 4.0,
        "magnetizing_inductance_h": 2.0e-3, "leakage_inductance_h": 20.0e-6,
        "clamp_voltage_v": 80.0, "switching_frequency_hz": 100_000.0,
        "input_voltage_v": 24.0, "duty_ratio": 0.4,
        "magnetizing_current_peak_a": 0.048, "leakage_current_peak_a": 1.6,
        "switch_voltage_peak_v": 104.0, "switch_current_peak_a": 1.6,
        "copper_loss_w": 0.32, "core_loss_w": 0.48, "clamp_loss_w": 0.16,
        "input_power_w": 10.96, "output_power_w": 10.0,
        "cycle_energy_j": 1.096e-4,
    }
    key = "flyback_turnsratio_magnetizing_leakage_clamp_stress_loss_power_energy_waveform_result_identity"
    positive[key] = {
        "flyback_generation_id": generation,
        **{name: generation for name in (
            "turnsratio_flyback_generation_id", "magnetizing_flyback_generation_id",
            "leakage_flyback_generation_id", "clamp_flyback_generation_id",
            "stress_flyback_generation_id", "loss_flyback_generation_id",
            "power_flyback_generation_id", "energy_flyback_generation_id",
            "waveform_flyback_generation_id", "result_flyback_generation_id",
        )},
        **values, **{f"result_{name}": value for name, value in values.items()},
        "waveform_owner": "flyback/waveform-843", "accepted_waveform_owner": "flyback/waveform-843",
        "waveform_sha256": "9" * 64, "accepted_waveform_sha256": "9" * 64,
        "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
    }
    generation = "opamp-stability-843"
    values = {
        "loop_gain_dc": 100_000.0, "crossover_frequency_hz": 1.0e6,
        "phase_margin_deg": 60.0, "input_step_v": 0.01, "output_step_v": 1.0,
        "overshoot_fraction": 0.10, "settling_time_s": 4.0e-6,
        "slew_rate_v_per_s": 400_000.0, "supply_voltage_v": 5.0,
        "supply_current_a": 2.0e-3, "supply_power_w": 10.0e-3,
    }
    key = "opamp_loopgain_phasemargin_crossover_step_overshoot_slew_power_result_identity"
    positive[key] = {
        "opamp_generation_id": generation,
        **{name: generation for name in (
            "loopgain_opamp_generation_id", "phasemargin_opamp_generation_id",
            "crossover_opamp_generation_id", "step_opamp_generation_id",
            "overshoot_opamp_generation_id", "slew_opamp_generation_id",
            "power_opamp_generation_id", "result_opamp_generation_id",
        )},
        **values, **{f"result_{name}": value for name, value in values.items()},
        "waveform_owner": "opamp/stability-843", "accepted_waveform_owner": "opamp/stability-843",
        "waveform_sha256": "b" * 64, "accepted_waveform_sha256": "b" * 64,
        "result_sha256": "c" * 64, "accepted_result_sha256": "c" * 64,
    }
    return payload


def test_v43_public_positive_flyback_and_opamp_closure():
    assert ideal_transformer_identity_gate(_v43())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v43_public_rejects_flyback_contract_mismatch():
    payload = _v43()
    identity = payload["metrics"]["positive"][
        "flyback_turnsratio_magnetizing_leakage_clamp_stress_loss_power_energy_waveform_result_identity"
    ]
    identity["result_input_power_w"] = 2.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v43_public_rejects_opamp_contract_mismatch():
    payload = _v43()
    identity = payload["metrics"]["positive"][
        "opamp_loopgain_phasemargin_crossover_step_overshoot_slew_power_result_identity"
    ]
    identity["result_phase_margin_deg"] = -5.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
