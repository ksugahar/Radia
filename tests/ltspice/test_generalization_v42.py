import math
from copy import deepcopy

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v41 import _v41


_PROMOTED_CASE_IDS = (
    "v42_public_mosfet_gatedrive_charge_current_deadtime_loss_temperature_energy_mismatch",
    "v42_public_instrumentationamp_gain_cmrr_inputrange_noise_outputheadroom_power_mismatch",
)
_GATE_KEY = (
    "mosfet_gatedrive_charge_current_deadtime_switching_conduction_loss_"
    "temperature_cycle_energy_waveform_result_identity"
)
_INSTRUMENTATION_KEY = (
    "instrumentation_amplifier_gain_cmrr_inputrange_noise_output_headroom_"
    "power_circuit_result_identity"
)


def _v42():
    payload = deepcopy(_v41())
    positive = payload["metrics"]["positive"]
    generation = "mosfet-gatedrive-842"
    charge, source_current, sink_current = 60.0e-9, 2.0, 2.4
    rise, fall, frequency = charge / source_current, charge / sink_current, 100_000.0
    voltage, current, duty, resistance = 400.0, 5.0, 0.5, 0.08
    switching_loss = 0.5 * voltage * current * (rise + fall) * frequency
    conduction_loss = current**2 * resistance * duty
    total_loss = switching_loss + conduction_loss
    mirrored = {
        "gate_charge_c": charge,
        "driver_source_current_a": source_current,
        "driver_sink_current_a": sink_current,
        "rise_time_s": rise,
        "fall_time_s": fall,
        "dead_time_s": 100.0e-9,
        "switching_frequency_hz": frequency,
        "drain_voltage_v": voltage,
        "drain_current_a": current,
        "duty_ratio": duty,
        "rds_on_ohm": resistance,
        "switching_loss_w": switching_loss,
        "conduction_loss_w": conduction_loss,
        "total_loss_w": total_loss,
        "ambient_temperature_c": 25.0,
        "junction_to_ambient_k_per_w": 2.0,
        "junction_temperature_c": 25.0 + 2.0 * total_loss,
        "cycle_loss_energy_j": total_loss / frequency,
    }
    positive[_GATE_KEY] = {
        "gatedrive_generation_id": generation,
        **{key: generation for key in (
            "charge_gatedrive_generation_id", "current_gatedrive_generation_id",
            "timing_gatedrive_generation_id", "loss_gatedrive_generation_id",
            "temperature_gatedrive_generation_id", "energy_gatedrive_generation_id",
            "waveform_gatedrive_generation_id", "result_gatedrive_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "waveform_owner": "gate-drive/waveform-842",
        "accepted_waveform_owner": "gate-drive/waveform-842",
        "waveform_sha256": "1" * 64, "accepted_waveform_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }

    generation = "instrumentation-amplifier-842"
    feedback, gain_resistance = 49_400.0, 499.0
    gain = 1.0 + feedback / gain_resistance
    common_mode, cmrr = 2.5, 100.0
    common_mode_gain = gain / (10.0 ** (cmrr / 20.0))
    signal_output = gain * 0.010
    noise_density, noise_bandwidth = 10.0e-9, 10_000.0
    mirrored = {
        "feedback_resistance_ohm": feedback,
        "gain_resistance_ohm": gain_resistance,
        "differential_gain_v_per_v": gain,
        "differential_input_v": 0.010,
        "common_mode_input_v": common_mode,
        "input_common_mode_min_v": 0.1,
        "input_common_mode_max_v": 4.9,
        "cmrr_db": cmrr,
        "common_mode_gain_v_per_v": common_mode_gain,
        "signal_output_v": signal_output,
        "common_mode_output_error_v": common_mode_gain * common_mode,
        "input_noise_density_v_per_sqrt_hz": noise_density,
        "noise_bandwidth_hz": noise_bandwidth,
        "integrated_output_noise_v_rms": noise_density * math.sqrt(noise_bandwidth),
        "output_low_limit_v": 0.1,
        "output_high_limit_v": 4.9,
        "output_headroom_v": min(signal_output - 0.1, 4.9 - signal_output),
        "supply_voltage_v": 5.0,
        "quiescent_supply_current_a": 1.0e-3,
        "supply_power_w": 5.0e-3,
    }
    positive[_INSTRUMENTATION_KEY] = {
        "instrumentation_generation_id": generation,
        **{key: generation for key in (
            "gain_instrumentation_generation_id", "cmrr_instrumentation_generation_id",
            "inputrange_instrumentation_generation_id", "noise_instrumentation_generation_id",
            "headroom_instrumentation_generation_id", "power_instrumentation_generation_id",
            "circuit_instrumentation_generation_id", "result_instrumentation_generation_id",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "circuit_owner": "instrumentation-amplifier/circuit-842",
        "accepted_circuit_owner": "instrumentation-amplifier/circuit-842",
        "circuit_sha256": "3" * 64, "accepted_circuit_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v42_public_positive_gate_drive_and_instrumentation_amplifier_closure():
    assert ideal_transformer_identity_gate(_v42())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v42_public_mosfet_gatedrive_charge_current_deadtime_loss_temperature_energy_mismatch():
    payload = _v42()
    identity = payload["metrics"]["positive"][_GATE_KEY]
    identity.update({
        "current_gatedrive_generation_id": "mosfet-gatedrive-841",
        "loss_gatedrive_generation_id": "mosfet-gatedrive-840",
        "result_rise_time_s": -1.0, "result_dead_time_s": 0.0,
        "result_switching_loss_w": -1.0, "result_conduction_loss_w": -1.0,
        "result_junction_temperature_c": 500.0,
        "result_cycle_loss_energy_j": -1.0,
        "accepted_waveform_owner": "gate-drive/old",
        "accepted_result_sha256": "a" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mosfet_gate_drives_use_current_charge_current_deadtime_losses_temperature_cycle_energy_waveform_and_result"
    ]


def test_v42_public_instrumentationamp_gain_cmrr_inputrange_noise_outputheadroom_power_mismatch():
    payload = _v42()
    identity = payload["metrics"]["positive"][_INSTRUMENTATION_KEY]
    identity.update({
        "gain_instrumentation_generation_id": "instrumentation-amplifier-841",
        "noise_instrumentation_generation_id": "instrumentation-amplifier-840",
        "result_differential_gain_v_per_v": -1.0, "result_cmrr_db": -10.0,
        "result_common_mode_input_v": 9.0,
        "result_integrated_output_noise_v_rms": -1.0,
        "result_output_headroom_v": -1.0, "result_supply_power_w": -1.0,
        "accepted_circuit_owner": "instrumentation-amplifier/old",
        "accepted_result_sha256": "b" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "instrumentation_amplifiers_use_current_gain_cmrr_input_range_noise_headroom_power_circuit_and_result"
    ]


def test_v42_public_rejects_self_consistent_short_deadtime():
    payload = _v42()
    identity = payload["metrics"]["positive"][_GATE_KEY]
    identity["dead_time_s"] = identity["result_dead_time_s"] = 1.0e-9
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_common_mode_outside_input_range():
    payload = _v42()
    identity = payload["metrics"]["positive"][_INSTRUMENTATION_KEY]
    identity["common_mode_input_v"] = identity["result_common_mode_input_v"] = 9.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
