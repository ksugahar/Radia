"""Focused LTSpice v43 artifact gates.

These validators deliberately require both the declared contract and the
replayed result to agree.  They are small enough to teach the identity rules
without pretending to be a full waveform parser.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _sha256(value: object) -> bool:
    text = str(value or "")
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text.lower())


def _numbers_equal(contract: Mapping[str, object], fields: Sequence[str]) -> bool:
    try:
        return all(
            math.isclose(
                float(contract[field]),
                float(contract[f"result_{field}"]),
                rel_tol=1.0e-12,
                abs_tol=1.0e-15,
            )
            for field in fields
        )
    except (KeyError, TypeError, ValueError):
        return False


def flyback_identity_ok(contract: Mapping[str, object]) -> bool:
    """Check a flyback magnetic/stress/loss/power identity contract."""

    generation = str(contract.get("flyback_generation_id") or "")
    generation_fields = (
        "turnsratio_flyback_generation_id", "magnetizing_flyback_generation_id",
        "leakage_flyback_generation_id", "clamp_flyback_generation_id",
        "stress_flyback_generation_id", "loss_flyback_generation_id",
        "power_flyback_generation_id", "energy_flyback_generation_id",
        "waveform_flyback_generation_id", "result_flyback_generation_id",
    )
    fields = (
        "primary_turns", "secondary_turns", "turns_ratio",
        "magnetizing_inductance_h", "leakage_inductance_h", "clamp_voltage_v",
        "switching_frequency_hz", "input_voltage_v", "duty_ratio",
        "magnetizing_current_peak_a", "leakage_current_peak_a",
        "switch_voltage_peak_v", "switch_current_peak_a", "copper_loss_w",
        "core_loss_w", "clamp_loss_w", "input_power_w", "output_power_w",
        "cycle_energy_j",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    positive = (
        bool(generation)
        and all(contract.get(field) == generation for field in generation_fields)
        and all(math.isfinite(value) for value in values.values())
        and all(values[field] > 0.0 for field in (
            "primary_turns", "secondary_turns", "magnetizing_inductance_h",
            "leakage_inductance_h", "clamp_voltage_v", "switching_frequency_hz",
            "input_voltage_v", "magnetizing_current_peak_a",
            "leakage_current_peak_a", "switch_voltage_peak_v",
            "switch_current_peak_a", "copper_loss_w", "core_loss_w",
            "clamp_loss_w", "input_power_w", "output_power_w",
            "cycle_energy_j",
        ))
        and 0.0 < values["duty_ratio"] < 1.0
        and math.isclose(
            values["turns_ratio"],
            values["primary_turns"] / values["secondary_turns"],
            rel_tol=1.0e-12,
        )
        and math.isclose(
            values["magnetizing_current_peak_a"],
            values["input_voltage_v"] * values["duty_ratio"]
            / (values["magnetizing_inductance_h"] * values["switching_frequency_hz"]),
            rel_tol=1.0e-12,
        )
        and math.isclose(
            values["leakage_current_peak_a"],
            0.1 * values["clamp_voltage_v"] * values["duty_ratio"]
            / (values["leakage_inductance_h"] * values["switching_frequency_hz"]),
            rel_tol=1.0e-12,
        )
        and math.isclose(
            values["switch_voltage_peak_v"],
            values["input_voltage_v"] + values["clamp_voltage_v"],
            rel_tol=1.0e-12,
        )
        and values["switch_current_peak_a"] >= values["leakage_current_peak_a"]
        and math.isclose(
            values["input_power_w"],
            values["output_power_w"] + values["copper_loss_w"]
            + values["core_loss_w"] + values["clamp_loss_w"],
            rel_tol=1.0e-12,
        )
        and math.isclose(
            values["cycle_energy_j"],
            values["input_power_w"] / values["switching_frequency_hz"],
            rel_tol=1.0e-12,
        )
        and _numbers_equal(contract, fields)
        and bool(str(contract.get("waveform_owner") or ""))
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _sha256(contract.get("waveform_sha256"))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _sha256(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )
    return positive


def opamp_stability_identity_ok(contract: Mapping[str, object]) -> bool:
    """Check loop-gain, transient stability, slew, and supply-power closure."""

    generation = str(contract.get("opamp_generation_id") or "")
    generation_fields = (
        "loopgain_opamp_generation_id", "phasemargin_opamp_generation_id",
        "crossover_opamp_generation_id", "step_opamp_generation_id",
        "overshoot_opamp_generation_id", "slew_opamp_generation_id",
        "power_opamp_generation_id", "result_opamp_generation_id",
    )
    fields = (
        "loop_gain_dc", "crossover_frequency_hz", "phase_margin_deg",
        "input_step_v", "output_step_v", "overshoot_fraction",
        "settling_time_s", "slew_rate_v_per_s", "supply_voltage_v",
        "supply_current_a", "supply_power_w",
    )
    try:
        values = {field: float(contract[field]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(contract.get(field) == generation for field in generation_fields)
        and all(math.isfinite(value) for value in values.values())
        and values["loop_gain_dc"] > 0.0
        and values["crossover_frequency_hz"] > 0.0
        and 0.0 < values["phase_margin_deg"] < 180.0
        and values["input_step_v"] > 0.0
        and values["output_step_v"] > 0.0
        and 0.0 <= values["overshoot_fraction"] < 1.0
        and values["settling_time_s"] > 0.0
        and values["slew_rate_v_per_s"] >= values["output_step_v"] / values["settling_time_s"]
        and values["supply_voltage_v"] > 0.0
        and values["supply_current_a"] > 0.0
        and math.isclose(
            values["supply_power_w"],
            values["supply_voltage_v"] * values["supply_current_a"],
            rel_tol=1.0e-12,
        )
        and _numbers_equal(contract, fields)
        and bool(str(contract.get("waveform_owner") or ""))
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _sha256(contract.get("waveform_sha256"))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _sha256(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )
