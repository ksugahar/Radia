"""Solver-neutral LTSpice v44 transient and noise identity gates."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _sha256(value: object) -> bool:
    text = str(value or "")
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text.lower())


def _mirror_ok(contract: Mapping[str, object], fields: Sequence[str]) -> bool:
    try:
        return all(math.isclose(float(contract[field]), float(contract[f"result_{field}"]), rel_tol=1e-12, abs_tol=1e-18) for field in fields)
    except (KeyError, TypeError, ValueError):
        return False


def _lineage_ok(contract: Mapping[str, object], generation: str, fields: Sequence[str], owner: str) -> bool:
    return (
        bool(generation)
        and all(contract.get(name) == generation for name in fields)
        and bool(str(contract.get(owner) or ""))
        and contract.get(f"accepted_{owner}") == contract.get(owner)
        and _sha256(contract.get("waveform_sha256"))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _sha256(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def buck_transient_identity_ok(contract: Mapping[str, object]) -> bool:
    fields = ("startup_window_start_s", "startup_window_stop_s", "inductor_current_avg_a", "inductor_current_ripple_a", "output_voltage_v", "output_ripple_v", "input_power_w", "output_power_w", "efficiency_fraction", "energy_start_j", "energy_end_j", "energy_balance_residual_j")
    try:
        values = {name: float(contract[name]) for name in fields}
    except (KeyError, TypeError, ValueError):
        return False
    return (
        _lineage_ok(contract, str(contract.get("buck_generation_id") or ""), ("startup_buck_generation_id", "inductor_current_buck_generation_id", "ripple_buck_generation_id", "efficiency_buck_generation_id", "energy_buck_generation_id", "waveform_buck_generation_id", "result_buck_generation_id"), "waveform_owner")
        and all(math.isfinite(value) for value in values.values())
        and values["startup_window_stop_s"] > values["startup_window_start_s"] >= 0.0
        and values["inductor_current_avg_a"] > 0.0 and values["inductor_current_ripple_a"] >= 0.0
        and values["output_voltage_v"] > 0.0 and values["output_ripple_v"] >= 0.0
        and values["input_power_w"] > 0.0 and values["output_power_w"] > 0.0
        and 0.0 < values["efficiency_fraction"] <= 1.0
        and values["energy_start_j"] >= 0.0 and values["energy_end_j"] >= 0.0
        and abs(values["energy_balance_residual_j"]) <= 1e-12
        and math.isclose(values["output_power_w"], values["input_power_w"] * values["efficiency_fraction"], rel_tol=1e-12)
        and _mirror_ok(contract, fields)
    )


def noise_ac_identity_ok(contract: Mapping[str, object]) -> bool:
    fields = ("frequency_hz", "transfer_magnitude_v_per_v", "input_referred_psd_v2_per_hz", "bandwidth_hz", "integrated_noise_v_rms", "correlation_coefficient", "output_noise_v_rms")
    try:
        values = {name: float(contract[name]) for name in fields}
    except (KeyError, TypeError, ValueError):
        return False
    return (
        _lineage_ok(contract, str(contract.get("noise_generation_id") or ""), ("transfer_noise_generation_id", "input_psd_noise_generation_id", "bandwidth_noise_generation_id", "correlation_noise_generation_id", "measure_noise_generation_id", "waveform_noise_generation_id", "result_noise_generation_id"), "measure_owner")
        and all(math.isfinite(value) for value in values.values())
        and values["frequency_hz"] > 0.0 and values["transfer_magnitude_v_per_v"] >= 0.0
        and values["input_referred_psd_v2_per_hz"] >= 0.0 and values["bandwidth_hz"] > 0.0
        and values["integrated_noise_v_rms"] >= 0.0 and -1.0 <= values["correlation_coefficient"] <= 1.0
        and values["output_noise_v_rms"] >= 0.0
        and math.isclose(values["integrated_noise_v_rms"] ** 2, values["input_referred_psd_v2_per_hz"] * values["bandwidth_hz"], rel_tol=1e-12, abs_tol=1e-18)
        and math.isclose(values["output_noise_v_rms"], values["integrated_noise_v_rms"] * values["transfer_magnitude_v_per_v"], rel_tol=1e-12, abs_tol=1e-18)
        and _mirror_ok(contract, fields)
    )


def validate_ltspice_v44_public_identity(positive: Mapping[str, object]) -> dict[str, bool]:
    """Return optional v44 checks while preserving older artifact compatibility."""
    buck_key = "buck_transient_startup_inductorcurrent_outputripple_efficiency_energy_waveform_result_identity"
    noise_key = "noise_ac_transfer_inputreferred_psd_bandwidth_correlation_measure_owner_identity"
    buck = positive.get(buck_key)
    noise = positive.get(noise_key)
    return {
        "buck_v44_transient_identity": buck is None or (isinstance(buck, Mapping) and buck_transient_identity_ok(buck)),
        "noise_v44_ac_identity": noise is None or (isinstance(noise, Mapping) and noise_ac_identity_ok(noise)),
    }
