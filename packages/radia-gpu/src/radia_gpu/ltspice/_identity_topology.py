"""Internal topology evidence checks for LTspice identities."""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from .ltspice_v43_gates import flyback_identity_ok, opamp_stability_identity_ok
from ._identity_common import (
    _is_sha256,
)


def _sallen_key_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "sallenkey_poles_q_gain_opamp_gbw_phase_noise_step_circuit_owner_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("filter_generation_id") or "")
    fields = (
        "r1_ohm", "r2_ohm", "c1_f", "c2_f", "noninverting_gain",
        "pole_frequency_hz", "quality_factor", "dc_gain", "opamp_gbw_hz",
        "phase_margin_deg", "integrated_output_noise_v_rms",
        "step_overshoot_fraction", "step_settling_time_s",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    r1, r2 = values["r1_ohm"], values["r2_ohm"]
    c1, c2 = values["c1_f"], values["c2_f"]
    gain = values["noninverting_gain"]
    scale = math.sqrt(r1 * r2 * c1 * c2) if min(r1, r2, c1, c2) > 0.0 else math.nan
    expected_frequency = 1.0 / (2.0 * math.pi * scale) if scale > 0.0 else math.nan
    q_denominator = c2 * (r1 + r2) + (1.0 - gain) * r1 * c1
    expected_q = scale / q_denominator if q_denominator > 0.0 else math.nan
    expected_overshoot = 0.0
    if math.isfinite(expected_q) and expected_q > 0.5:
        damping = 1.0 / (2.0 * expected_q)
        expected_overshoot = math.exp(
            -math.pi * damping / math.sqrt(1.0 - damping * damping)
        )
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "pole_filter_generation_id", "q_filter_generation_id",
            "gain_filter_generation_id", "opamp_filter_generation_id",
            "noise_filter_generation_id", "step_filter_generation_id",
            "circuit_filter_generation_id", "result_filter_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(r1, r2, c1, c2) > 0.0
        and 0.0 < gain < 3.0
        and math.isclose(values["pole_frequency_hz"], expected_frequency, rel_tol=1.0e-9, abs_tol=0.0)
        and math.isclose(values["quality_factor"], expected_q, rel_tol=1.0e-9, abs_tol=0.0)
        and math.isclose(values["dc_gain"], gain, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and values["opamp_gbw_hz"] >= 100.0 * expected_frequency
        and 45.0 <= values["phase_margin_deg"] <= 180.0
        and values["integrated_output_noise_v_rms"] >= 0.0
        and 0.0 <= values["step_overshoot_fraction"] < 1.0
        and math.isclose(values["step_overshoot_fraction"], expected_overshoot, rel_tol=1.0e-6, abs_tol=1.0e-9)
        and 0.0 < values["step_settling_time_s"] <= 20.0 / expected_frequency
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-15) for field in fields)
        and bool(str(contract.get("circuit_owner") or ""))
        and contract.get("accepted_circuit_owner") == contract.get("circuit_owner")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _bridge_rectifier_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "bridge_rectifier_inrush_ripple_charge_conduction_loss_load_cycle_energy_waveform_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("rectifier_generation_id") or "")
    fields = (
        "line_frequency_hz", "ripple_frequency_hz", "capacitor_f",
        "average_output_voltage_v", "load_resistance_ohm",
        "average_load_current_a", "ripple_voltage_pp_v",
        "inrush_peak_current_a", "capacitor_charge_per_cycle_c",
        "rectifier_delivered_charge_per_cycle_c", "diode_conduction_fraction",
        "diode_forward_drop_v", "average_diode_loss_w",
        "average_load_power_w", "ripple_cycle_period_s",
        "source_energy_per_cycle_j", "load_energy_per_cycle_j",
        "diode_energy_per_cycle_j",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    line_frequency = values["line_frequency_hz"]
    ripple_frequency = values["ripple_frequency_hz"]
    capacitance = values["capacitor_f"]
    voltage = values["average_output_voltage_v"]
    resistance = values["load_resistance_ohm"]
    expected_current = voltage / resistance if resistance > 0.0 else math.nan
    expected_power = voltage * expected_current
    expected_ripple = expected_current / (ripple_frequency * capacitance) if min(ripple_frequency, capacitance) > 0.0 else math.nan
    expected_charge = expected_current / ripple_frequency if ripple_frequency > 0.0 else math.nan
    expected_diode_loss = 2.0 * values["diode_forward_drop_v"] * expected_current
    expected_period = 1.0 / ripple_frequency if ripple_frequency > 0.0 else math.nan
    expected_load_energy = expected_power * expected_period
    expected_diode_energy = expected_diode_loss * expected_period
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "inrush_rectifier_generation_id", "ripple_rectifier_generation_id",
            "charge_rectifier_generation_id", "conduction_rectifier_generation_id",
            "loss_rectifier_generation_id", "load_rectifier_generation_id",
            "energy_rectifier_generation_id", "waveform_rectifier_generation_id",
            "result_rectifier_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(line_frequency, capacitance, voltage, resistance) > 0.0
        and math.isclose(ripple_frequency, 2.0 * line_frequency, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["average_load_current_a"], expected_current, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["average_load_power_w"], expected_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["ripple_voltage_pp_v"], expected_ripple, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and values["inrush_peak_current_a"] >= expected_current
        and math.isclose(values["capacitor_charge_per_cycle_c"], expected_charge, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["rectifier_delivered_charge_per_cycle_c"], expected_charge, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and 0.0 < values["diode_conduction_fraction"] < 1.0
        and values["diode_forward_drop_v"] >= 0.0
        and math.isclose(values["average_diode_loss_w"], expected_diode_loss, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["ripple_cycle_period_s"], expected_period, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["load_energy_per_cycle_j"], expected_load_energy, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["diode_energy_per_cycle_j"], expected_diode_energy, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["source_energy_per_cycle_j"], expected_load_energy + expected_diode_energy, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-15) for field in fields)
        and bool(str(contract.get("waveform_owner") or ""))
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _flyback_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "flyback_magnetizing_leakage_snubber_demag_flux_current_cycle_energy_waveform_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("flyback_generation_id") or "")
    fields = (
        "switching_frequency_hz", "switching_period_s", "input_voltage_v",
        "duty_cycle", "magnetizing_inductance_h", "leakage_inductance_h",
        "peak_switch_current_a", "demagnetization_interval_s",
        "peak_core_flux_density_t", "snubber_energy_per_cycle_j",
        "average_snubber_dissipation_w", "magnetizing_energy_per_cycle_j",
        "delivered_energy_per_cycle_j", "source_energy_per_cycle_j",
        "cycle_energy_residual_j",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    frequency = values["switching_frequency_hz"]
    period = values["switching_period_s"]
    voltage = values["input_voltage_v"]
    duty = values["duty_cycle"]
    magnetizing = values["magnetizing_inductance_h"]
    leakage = values["leakage_inductance_h"]
    peak_current = values["peak_switch_current_a"]
    expected_current = voltage * duty * period / magnetizing if magnetizing > 0.0 else math.nan
    expected_magnetizing_energy = 0.5 * magnetizing * peak_current**2
    expected_leakage_energy = 0.5 * leakage * peak_current**2
    expected_residual = (
        values["source_energy_per_cycle_j"]
        - values["delivered_energy_per_cycle_j"]
        - values["snubber_energy_per_cycle_j"]
    )
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "magnetizing_flyback_generation_id", "leakage_flyback_generation_id",
            "snubber_flyback_generation_id", "demag_flyback_generation_id",
            "flux_flyback_generation_id", "current_flyback_generation_id",
            "energy_flyback_generation_id", "waveform_flyback_generation_id",
            "result_flyback_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(frequency, period, voltage, magnetizing, leakage) > 0.0
        and 0.0 < duty < 1.0
        and math.isclose(period, 1.0 / frequency, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isclose(peak_current, expected_current, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and 0.0 < values["demagnetization_interval_s"] <= period * (1.0 - duty)
        and 0.0 < values["peak_core_flux_density_t"] <= 2.5
        and math.isclose(values["snubber_energy_per_cycle_j"], expected_leakage_energy, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isclose(values["average_snubber_dissipation_w"], expected_leakage_energy * frequency, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["magnetizing_energy_per_cycle_j"], expected_magnetizing_energy, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isclose(values["delivered_energy_per_cycle_j"], expected_magnetizing_energy, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isclose(values["cycle_energy_residual_j"], expected_residual, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and abs(expected_residual) <= 1.0e-12 * max(values["source_energy_per_cycle_j"], 1.0e-18)
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-18) for field in fields)
        and bool(str(contract.get("waveform_owner") or ""))
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _transimpedance_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "transimpedance_photodiode_capacitance_noise_gain_bandwidth_stability_noise_step_circuit_result_identity"
    )
    return _transimpedance_contract_ok(contract)


def _llc_resonant_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "llc_resonant_elements_frequency_gain_zvs_current_loss_cycle_energy_waveform_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("llc_generation_id") or "")
    fields = (
        "resonant_inductance_h", "resonant_capacitance_f",
        "magnetizing_inductance_h", "resonant_frequency_hz",
        "switching_frequency_hz", "normalized_switching_frequency",
        "input_voltage_v", "transformer_turns_ratio", "conversion_gain_v_per_v",
        "output_voltage_v", "dead_time_s", "switch_output_capacitance_f",
        "magnetizing_current_at_transition_a", "required_zvs_current_a",
        "circulating_current_rms_a", "switch_conduction_resistance_ohm",
        "conduction_loss_w", "switching_loss_w", "device_loss_w", "tank_loss_w",
        "input_energy_per_cycle_j", "output_energy_per_cycle_j",
        "device_loss_energy_per_cycle_j", "tank_loss_energy_per_cycle_j",
        "cycle_energy_residual_j",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    resonant_frequency = 1.0 / (
        2.0 * math.pi * math.sqrt(
            values["resonant_inductance_h"] * values["resonant_capacitance_f"]
        )
    ) if min(values["resonant_inductance_h"], values["resonant_capacitance_f"]) > 0.0 else math.nan
    required_zvs_current = (
        2.0 * values["switch_output_capacitance_f"] * values["input_voltage_v"]
        / values["dead_time_s"]
        if values["dead_time_s"] > 0.0
        else math.nan
    )
    expected_residual = (
        values["input_energy_per_cycle_j"]
        - values["output_energy_per_cycle_j"]
        - values["device_loss_energy_per_cycle_j"]
        - values["tank_loss_energy_per_cycle_j"]
    )
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "element_llc_generation_id", "frequency_llc_generation_id",
            "gain_llc_generation_id", "zvs_llc_generation_id",
            "current_llc_generation_id", "loss_llc_generation_id",
            "energy_llc_generation_id", "waveform_llc_generation_id",
            "result_llc_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(
            values["resonant_inductance_h"], values["resonant_capacitance_f"],
            values["magnetizing_inductance_h"], values["switching_frequency_hz"],
            values["input_voltage_v"], values["transformer_turns_ratio"],
            values["conversion_gain_v_per_v"], values["dead_time_s"],
            values["switch_output_capacitance_f"],
            values["magnetizing_current_at_transition_a"],
            values["circulating_current_rms_a"],
            values["switch_conduction_resistance_ohm"],
        ) > 0.0
        and values["magnetizing_inductance_h"] > values["resonant_inductance_h"]
        and math.isclose(values["resonant_frequency_hz"], resonant_frequency, rel_tol=1.0e-12, abs_tol=1.0e-6)
        and math.isclose(values["normalized_switching_frequency"], values["switching_frequency_hz"] / resonant_frequency, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["output_voltage_v"], values["input_voltage_v"] * values["transformer_turns_ratio"] * values["conversion_gain_v_per_v"], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["required_zvs_current_a"], required_zvs_current, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and contract.get("zvs_condition_met") is True
        and contract.get("result_zvs_condition_met") is True
        and values["magnetizing_current_at_transition_a"] >= required_zvs_current
        and math.isclose(values["conduction_loss_w"], values["circulating_current_rms_a"] ** 2 * values["switch_conduction_resistance_ohm"], rel_tol=1.0e-12, abs_tol=1.0e-15)
        and values["switching_loss_w"] >= 0.0 and values["tank_loss_w"] >= 0.0
        and math.isclose(values["device_loss_w"], values["conduction_loss_w"] + values["switching_loss_w"], rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["device_loss_energy_per_cycle_j"], values["device_loss_w"] / values["switching_frequency_hz"], rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isclose(values["tank_loss_energy_per_cycle_j"], values["tank_loss_w"] / values["switching_frequency_hz"], rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isclose(values["cycle_energy_residual_j"], expected_residual, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and abs(expected_residual) <= 1.0e-12 * max(values["input_energy_per_cycle_j"], 1.0e-18)
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-18) for field in fields)
        and bool(str(contract.get("waveform_owner") or ""))
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _bjt_amplifier_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "bjt_bias_gm_gain_pole_noise_distortion_thermal_power_circuit_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("bjt_generation_id") or "")
    fields = (
        "collector_current_a", "collector_emitter_voltage_v", "thermal_voltage_v",
        "transconductance_s", "collector_resistance_ohm", "load_resistance_ohm",
        "voltage_gain_v_per_v", "output_capacitance_f",
        "dominant_pole_frequency_hz", "input_noise_density_v_per_sqrt_hz",
        "integrated_output_noise_v_rms", "fundamental_output_v",
        "second_harmonic_v", "third_harmonic_v", "total_harmonic_distortion",
        "ambient_temperature_c", "junction_to_ambient_k_per_w",
        "device_power_w", "junction_temperature_c",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    parallel_resistance = 1.0 / (
        1.0 / values["collector_resistance_ohm"] + 1.0 / values["load_resistance_ohm"]
    ) if min(values["collector_resistance_ohm"], values["load_resistance_ohm"]) > 0.0 else math.nan
    expected_gm = values["collector_current_a"] / values["thermal_voltage_v"] if values["thermal_voltage_v"] > 0.0 else math.nan
    expected_gain = -expected_gm * parallel_resistance
    expected_pole = 1.0 / (2.0 * math.pi * parallel_resistance * values["output_capacitance_f"]) if values["output_capacitance_f"] > 0.0 else math.nan
    expected_noise = values["input_noise_density_v_per_sqrt_hz"] * abs(expected_gain) * math.sqrt(expected_pole)
    expected_thd = math.sqrt(values["second_harmonic_v"] ** 2 + values["third_harmonic_v"] ** 2) / values["fundamental_output_v"] if values["fundamental_output_v"] > 0.0 else math.nan
    expected_power = values["collector_emitter_voltage_v"] * values["collector_current_a"]
    expected_temperature = values["ambient_temperature_c"] + values["junction_to_ambient_k_per_w"] * expected_power
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "bias_bjt_generation_id", "gm_bjt_generation_id", "gain_bjt_generation_id",
            "pole_bjt_generation_id", "noise_bjt_generation_id",
            "distortion_bjt_generation_id", "thermal_bjt_generation_id",
            "circuit_bjt_generation_id", "result_bjt_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(
            values["collector_current_a"], values["collector_emitter_voltage_v"],
            values["thermal_voltage_v"], values["collector_resistance_ohm"],
            values["load_resistance_ohm"], values["output_capacitance_f"],
            values["input_noise_density_v_per_sqrt_hz"], values["fundamental_output_v"],
            values["junction_to_ambient_k_per_w"],
        ) > 0.0
        and math.isclose(values["transconductance_s"], expected_gm, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["voltage_gain_v_per_v"], expected_gain, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["dominant_pole_frequency_hz"], expected_pole, rel_tol=1.0e-12, abs_tol=1.0e-6)
        and math.isclose(values["integrated_output_noise_v_rms"], expected_noise, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and values["second_harmonic_v"] >= 0.0 and values["third_harmonic_v"] >= 0.0
        and math.isclose(values["total_harmonic_distortion"], expected_thd, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["device_power_w"], expected_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["junction_temperature_c"], expected_temperature, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and values["junction_temperature_c"] > values["ambient_temperature_c"]
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-15) for field in fields)
        and bool(str(contract.get("circuit_owner") or ""))
        and contract.get("accepted_circuit_owner") == contract.get("circuit_owner")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _boost_converter_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "boost_duty_inductor_ripple_ccm_output_ripple_efficiency_stress_cycle_energy_waveform_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("boost_generation_id") or "")
    fields = (
        "input_voltage_v", "output_voltage_v", "duty_ratio",
        "switching_frequency_hz", "inductance_h",
        "inductor_average_current_a", "inductor_ripple_peak_to_peak_a",
        "ccm_boundary_current_a", "output_capacitance_f", "output_current_a",
        "output_voltage_ripple_peak_to_peak_v", "efficiency",
        "switch_voltage_stress_v", "diode_reverse_voltage_stress_v",
        "switch_peak_current_a", "input_energy_per_cycle_j",
        "output_energy_per_cycle_j", "loss_energy_per_cycle_j",
        "cycle_energy_residual_j",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    vin = values["input_voltage_v"]
    vout = values["output_voltage_v"]
    duty = values["duty_ratio"]
    frequency = values["switching_frequency_hz"]
    inductance = values["inductance_h"]
    capacitance = values["output_capacitance_f"]
    output_current = values["output_current_a"]
    efficiency = values["efficiency"]
    expected_duty = 1.0 - vin / vout if vout > 0.0 else math.nan
    expected_ripple = vin * duty / (inductance * frequency) if min(inductance, frequency) > 0.0 else math.nan
    output_power = vout * output_current
    input_power = output_power / efficiency if efficiency > 0.0 else math.nan
    expected_average = input_power / vin if vin > 0.0 else math.nan
    expected_output_ripple = output_current * duty / (capacitance * frequency) if min(capacitance, frequency) > 0.0 else math.nan
    expected_residual = (
        values["input_energy_per_cycle_j"]
        - values["output_energy_per_cycle_j"]
        - values["loss_energy_per_cycle_j"]
    )
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "duty_boost_generation_id", "ripple_boost_generation_id",
            "ccm_boost_generation_id", "output_boost_generation_id",
            "efficiency_boost_generation_id", "stress_boost_generation_id",
            "energy_boost_generation_id", "waveform_boost_generation_id",
            "result_boost_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(vin, vout, frequency, inductance, capacitance, output_current) > 0.0
        and vout > vin and 0.0 < duty < 1.0 and 0.0 < efficiency <= 1.0
        and math.isclose(duty, expected_duty, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["inductor_ripple_peak_to_peak_a"], expected_ripple, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["ccm_boundary_current_a"], expected_ripple / 2.0, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["inductor_average_current_a"], expected_average, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and contract.get("ccm_condition_met") is True
        and contract.get("result_ccm_condition_met") is True
        and expected_average > expected_ripple / 2.0
        and math.isclose(values["output_voltage_ripple_peak_to_peak_v"], expected_output_ripple, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and values["switch_voltage_stress_v"] >= vout
        and values["diode_reverse_voltage_stress_v"] >= vout
        and math.isclose(values["switch_peak_current_a"], expected_average + expected_ripple / 2.0, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["input_energy_per_cycle_j"], input_power / frequency, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isclose(values["output_energy_per_cycle_j"], output_power / frequency, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and values["loss_energy_per_cycle_j"] >= 0.0
        and math.isclose(values["cycle_energy_residual_j"], expected_residual, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and abs(expected_residual) <= 1.0e-12 * max(values["input_energy_per_cycle_j"], 1.0e-18)
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-18) for field in fields)
        and bool(str(contract.get("waveform_owner") or ""))
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _active_filter_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "active_filter_pole_zero_q_gain_noise_slew_saturation_power_circuit_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("filter_generation_id") or "")
    fields = (
        "natural_frequency_rad_s", "quality_factor", "passband_gain_v_per_v",
        "input_noise_density_v_per_sqrt_hz", "noise_bandwidth_hz",
        "integrated_output_noise_v_rms", "input_frequency_hz",
        "input_amplitude_v_peak", "output_amplitude_v_peak",
        "slew_rate_demand_v_per_s", "available_slew_rate_v_per_s",
        "positive_supply_v", "negative_supply_v", "saturation_margin_v",
        "quiescent_supply_current_a", "supply_power_w",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
        poles = tuple(complex(float(row[0]), float(row[1])) for row in contract["pole_locations_rad_s"])
        result_poles = tuple(complex(float(row[0]), float(row[1])) for row in contract["result_pole_locations_rad_s"])
        zeros = tuple(complex(float(row[0]), float(row[1])) for row in contract["zero_locations_rad_s"])
        result_zeros = tuple(complex(float(row[0]), float(row[1])) for row in contract["result_zero_locations_rad_s"])
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    natural_frequency = values["natural_frequency_rad_s"]
    quality = values["quality_factor"]
    gain = values["passband_gain_v_per_v"]
    output_peak = values["output_amplitude_v_peak"]
    expected_noise = values["input_noise_density_v_per_sqrt_hz"] * math.sqrt(values["noise_bandwidth_hz"])
    expected_slew = 2.0 * math.pi * values["input_frequency_hz"] * output_peak
    expected_margin = min(
        values["positive_supply_v"] - output_peak,
        -values["negative_supply_v"] - output_peak,
    )
    expected_power = (
        values["positive_supply_v"] - values["negative_supply_v"]
    ) * values["quiescent_supply_current_a"]
    pole_model_ok = (
        len(poles) == 2
        and poles[0].real < 0.0
        and math.isclose(poles[0].real, poles[1].real, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(poles[0].imag, -poles[1].imag, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(abs(poles[0]), natural_frequency, rel_tol=1.0e-12, abs_tol=1.0e-9)
        and math.isclose(natural_frequency / (-2.0 * poles[0].real), quality, rel_tol=1.0e-12, abs_tol=1.0e-15)
    )
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "pole_filter_generation_id", "zero_filter_generation_id",
            "q_filter_generation_id", "gain_filter_generation_id",
            "noise_filter_generation_id", "slew_filter_generation_id",
            "saturation_filter_generation_id", "power_filter_generation_id",
            "circuit_filter_generation_id", "result_filter_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(
            natural_frequency, quality, gain,
            values["input_noise_density_v_per_sqrt_hz"], values["noise_bandwidth_hz"],
            values["input_frequency_hz"], values["input_amplitude_v_peak"],
            values["available_slew_rate_v_per_s"],
            values["quiescent_supply_current_a"],
        ) > 0.0
        and pole_model_ok and poles == result_poles and zeros == result_zeros
        and math.isclose(output_peak, gain * values["input_amplitude_v_peak"], rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["integrated_output_noise_v_rms"], expected_noise, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isclose(values["slew_rate_demand_v_per_s"], expected_slew, rel_tol=1.0e-12, abs_tol=1.0e-9)
        and values["available_slew_rate_v_per_s"] > expected_slew
        and values["positive_supply_v"] > 0.0 > values["negative_supply_v"]
        and math.isclose(values["saturation_margin_v"], expected_margin, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and expected_margin > 0.0
        and math.isclose(values["supply_power_w"], expected_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-18) for field in fields)
        and bool(str(contract.get("circuit_owner") or ""))
        and contract.get("accepted_circuit_owner") == contract.get("circuit_owner")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _mosfet_gatedrive_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "mosfet_gatedrive_charge_current_deadtime_switching_conduction_loss_temperature_cycle_energy_waveform_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("gatedrive_generation_id") or "")
    fields = (
        "gate_charge_c", "driver_source_current_a", "driver_sink_current_a",
        "rise_time_s", "fall_time_s", "dead_time_s", "switching_frequency_hz",
        "drain_voltage_v", "drain_current_a", "duty_ratio", "rds_on_ohm",
        "switching_loss_w", "conduction_loss_w", "total_loss_w",
        "ambient_temperature_c", "junction_to_ambient_k_per_w",
        "junction_temperature_c", "cycle_loss_energy_j",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    charge = values["gate_charge_c"]
    source_current = values["driver_source_current_a"]
    sink_current = values["driver_sink_current_a"]
    frequency = values["switching_frequency_hz"]
    expected_rise = charge / source_current if source_current > 0.0 else math.nan
    expected_fall = charge / sink_current if sink_current > 0.0 else math.nan
    expected_switching = (
        0.5 * values["drain_voltage_v"] * values["drain_current_a"]
        * (values["rise_time_s"] + values["fall_time_s"]) * frequency
    )
    expected_conduction = (
        values["drain_current_a"] ** 2 * values["rds_on_ohm"] * values["duty_ratio"]
    )
    expected_total = expected_switching + expected_conduction
    expected_temperature = (
        values["ambient_temperature_c"]
        + expected_total * values["junction_to_ambient_k_per_w"]
    )
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "charge_gatedrive_generation_id", "current_gatedrive_generation_id",
            "timing_gatedrive_generation_id", "loss_gatedrive_generation_id",
            "temperature_gatedrive_generation_id", "energy_gatedrive_generation_id",
            "waveform_gatedrive_generation_id", "result_gatedrive_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(
            charge, source_current, sink_current, frequency,
            values["drain_voltage_v"], values["drain_current_a"],
            values["rds_on_ohm"], values["junction_to_ambient_k_per_w"],
        ) > 0.0
        and 0.0 < values["duty_ratio"] < 1.0
        and math.isclose(values["rise_time_s"], expected_rise, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isclose(values["fall_time_s"], expected_fall, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and values["dead_time_s"] >= max(expected_rise, expected_fall) > 0.0
        and math.isclose(values["switching_loss_w"], expected_switching, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["conduction_loss_w"], expected_conduction, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["total_loss_w"], expected_total, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["junction_temperature_c"], expected_temperature, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and values["junction_temperature_c"] > values["ambient_temperature_c"]
        and math.isclose(values["cycle_loss_energy_j"], expected_total / frequency, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-18) for field in fields)
        and str(contract.get("waveform_owner") or "").startswith("gate-drive/")
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _instrumentation_amplifier_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "instrumentation_amplifier_gain_cmrr_inputrange_noise_output_headroom_power_circuit_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("instrumentation_generation_id") or "")
    fields = (
        "feedback_resistance_ohm", "gain_resistance_ohm",
        "differential_gain_v_per_v", "differential_input_v",
        "common_mode_input_v", "input_common_mode_min_v",
        "input_common_mode_max_v", "cmrr_db", "common_mode_gain_v_per_v",
        "signal_output_v", "common_mode_output_error_v",
        "input_noise_density_v_per_sqrt_hz", "noise_bandwidth_hz",
        "integrated_output_noise_v_rms", "output_low_limit_v",
        "output_high_limit_v", "output_headroom_v", "supply_voltage_v",
        "quiescent_supply_current_a", "supply_power_w",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    gain = values["differential_gain_v_per_v"]
    expected_gain = 1.0 + values["feedback_resistance_ohm"] / values["gain_resistance_ohm"]
    expected_common_mode_gain = gain / (10.0 ** (values["cmrr_db"] / 20.0))
    expected_output = gain * values["differential_input_v"]
    expected_cm_error = expected_common_mode_gain * values["common_mode_input_v"]
    expected_noise = values["input_noise_density_v_per_sqrt_hz"] * math.sqrt(values["noise_bandwidth_hz"])
    expected_headroom = min(
        expected_output - values["output_low_limit_v"],
        values["output_high_limit_v"] - expected_output,
    )
    expected_power = values["supply_voltage_v"] * values["quiescent_supply_current_a"]
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "gain_instrumentation_generation_id", "cmrr_instrumentation_generation_id",
            "inputrange_instrumentation_generation_id", "noise_instrumentation_generation_id",
            "headroom_instrumentation_generation_id", "power_instrumentation_generation_id",
            "circuit_instrumentation_generation_id", "result_instrumentation_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(
            values["feedback_resistance_ohm"], values["gain_resistance_ohm"],
            gain, values["input_noise_density_v_per_sqrt_hz"],
            values["noise_bandwidth_hz"], values["supply_voltage_v"],
            values["quiescent_supply_current_a"],
        ) > 0.0
        and values["input_common_mode_min_v"] <= values["common_mode_input_v"] <= values["input_common_mode_max_v"]
        and values["cmrr_db"] > 0.0
        and math.isclose(gain, expected_gain, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["common_mode_gain_v_per_v"], expected_common_mode_gain, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["signal_output_v"], expected_output, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["common_mode_output_error_v"], expected_cm_error, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["integrated_output_noise_v_rms"], expected_noise, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and values["output_low_limit_v"] < expected_output < values["output_high_limit_v"]
        and math.isclose(values["output_headroom_v"], expected_headroom, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and expected_headroom > 0.0
        and math.isclose(values["supply_power_w"], expected_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-18) for field in fields)
        and str(contract.get("circuit_owner") or "").startswith("instrumentation-amplifier/")
        and contract.get("accepted_circuit_owner") == contract.get("circuit_owner")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _transimpedance_contract_ok(contract: object) -> bool:
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("transimpedance_generation_id") or "")
    fields = (
        "feedback_resistance_ohm", "feedback_capacitance_f",
        "photodiode_capacitance_f", "opamp_input_capacitance_f",
        "transimpedance_dc_gain_v_per_a", "high_frequency_noise_gain",
        "closed_loop_bandwidth_hz", "opamp_gain_bandwidth_hz",
        "phase_margin_deg", "input_referred_current_noise_a_per_sqrt_hz",
        "integrated_output_noise_v_rms", "step_rise_time_s",
        "step_overshoot_fraction", "step_settling_time_s",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
        results = {field: float(contract[f"result_{field}"]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    resistance = values["feedback_resistance_ohm"]
    feedback_capacitance = values["feedback_capacitance_f"]
    detector_capacitance = values["photodiode_capacitance_f"]
    input_capacitance = values["opamp_input_capacitance_f"]
    expected_noise_gain = 1.0 + (detector_capacitance + input_capacitance) / feedback_capacitance if feedback_capacitance > 0.0 else math.nan
    expected_bandwidth = 1.0 / (2.0 * math.pi * resistance * feedback_capacitance) if min(resistance, feedback_capacitance) > 0.0 else math.nan
    return (
        bool(generation)
        and all(contract.get(key) == generation for key in (
            "capacitance_transimpedance_generation_id",
            "noisegain_transimpedance_generation_id",
            "bandwidth_transimpedance_generation_id",
            "stability_transimpedance_generation_id",
            "noise_transimpedance_generation_id", "step_transimpedance_generation_id",
            "circuit_transimpedance_generation_id", "result_transimpedance_generation_id",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(resistance, feedback_capacitance) > 0.0
        and detector_capacitance >= 0.0 and input_capacitance >= 0.0
        and math.isclose(values["transimpedance_dc_gain_v_per_a"], resistance, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["high_frequency_noise_gain"], expected_noise_gain, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["closed_loop_bandwidth_hz"], expected_bandwidth, rel_tol=1.0e-12, abs_tol=1.0e-9)
        and values["opamp_gain_bandwidth_hz"] >= expected_noise_gain * expected_bandwidth
        and 45.0 <= values["phase_margin_deg"] <= 180.0
        and values["input_referred_current_noise_a_per_sqrt_hz"] >= 0.0
        and values["integrated_output_noise_v_rms"] >= 0.0
        and math.isclose(values["step_rise_time_s"], 0.44 / expected_bandwidth, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and 0.0 <= values["step_overshoot_fraction"] < 1.0
        and values["step_settling_time_s"] >= values["step_rise_time_s"]
        and all(math.isclose(results[field], values[field], rel_tol=1.0e-12, abs_tol=1.0e-18) for field in fields)
        and bool(str(contract.get("circuit_owner") or ""))
        and contract.get("accepted_circuit_owner") == contract.get("circuit_owner")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _flyback_v43_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "flyback_turnsratio_magnetizing_leakage_clamp_stress_loss_power_energy_waveform_result_identity"
    )
    return contract is None or (
        isinstance(contract, Mapping) and flyback_identity_ok(contract)
    )


def _opamp_stability_v43_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "opamp_loopgain_phasemargin_crossover_step_overshoot_slew_power_result_identity"
    )
    return contract is None or (
        isinstance(contract, Mapping) and opamp_stability_identity_ok(contract)
    )
