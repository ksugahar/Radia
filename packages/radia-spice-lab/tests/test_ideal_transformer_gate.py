from __future__ import annotations

import copy
import math

import pytest

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate


def _pair(value: complex) -> list[float]:
    return [value.real, value.imag]


def _summary() -> dict:
    n = 4.0
    source_peak = 170.0
    source_rms = source_peak / math.sqrt(2.0)
    source_resistance = 1.0
    load_resistance = 100.0
    reflected = n**2 * load_resistance
    primary_current = source_rms / (source_resistance + reflected)
    primary_voltage = primary_current * reflected
    secondary_voltage = -primary_voltage / n
    secondary_current = -n * primary_current
    source_voltage = primary_voltage + source_resistance * primary_current
    transformer_power = primary_voltage * primary_current
    source_power = source_voltage * primary_current
    series_loss = source_resistance * primary_current**2
    positive = {
        "point_count": 1042,
        "fit_point_count": 683,
        "fit_window_start_s": 1.0 / 60.0,
        "fit_window_stop_s": 3.0 / 60.0,
        "source_voltage_phasor_rms_v": _pair(complex(source_voltage)),
        "primary_voltage_phasor_rms_v": _pair(complex(primary_voltage)),
        "secondary_voltage_phasor_rms_v": _pair(complex(secondary_voltage)),
        "source_delivery_current_phasor_rms_a": _pair(complex(primary_current)),
        "primary_current_phasor_rms_a": _pair(complex(primary_current)),
        "secondary_current_phasor_rms_a": _pair(complex(secondary_current)),
        "source_voltage_relative_error": 0.0,
        "primary_voltage_relative_error": 0.0,
        "secondary_voltage_relative_error": 0.0,
        "primary_current_relative_error": 0.0,
        "secondary_current_relative_error": 0.0,
        "voltage_turns_identity_relative_error": 0.0,
        "current_turns_identity_relative_error": 0.0,
        "secondary_load_ohm_law_relative_error": 0.0,
        "reflected_load_relative_error": 0.0,
        "source_series_kvl_relative_error": 0.0,
        "source_primary_current_relative_error": 0.0,
        "transformer_complex_power_relative_error": 0.0,
        "source_power_balance_relative_error": 0.0,
        "instantaneous_transformer_power_relative_error": 0.0,
        "instantaneous_source_balance_relative_error": 0.0,
        "maximum_phasor_fit_relative_residual": 1.0e-8,
        "transformer_input_complex_power_va": _pair(complex(transformer_power)),
        "transformer_output_complex_power_va": _pair(complex(transformer_power)),
        "source_delivery_complex_power_va": _pair(complex(source_power)),
        "series_loss_w": series_loss,
    }
    return {
        "model_contract": {
            "topology": "ideal_two_winding_transformer_with_series_source_resistance",
            "turns_ratio_primary_to_secondary": n,
            "source_offset_v": 0.0,
            "source_peak_v": source_peak,
            "frequency_hz": 60.0,
            "series_resistance_ohm": source_resistance,
            "load_resistance_ohm": load_resistance,
            "expected_source_voltage_rms_v": source_rms,
            "expected_reflected_load_ohm": reflected,
            "expected_primary_voltage_rms_v": primary_voltage,
            "expected_secondary_voltage_rms_v": abs(secondary_voltage),
            "expected_primary_current_rms_a": primary_current,
            "expected_secondary_current_rms_a": abs(secondary_current),
            "expected_load_power_w": transformer_power,
            "expected_series_loss_w": series_loss,
        },
        "metrics": {
            "maximum_phasor_replay_relative_error": 0.0,
            "positive": positive,
        },
        "timing_breakdown_s": {
            "preflight": 0.01,
            "solve": 0.1,
            "analysis": 0.01,
            "serialization": 0.01,
        },
    }


def test_accepts_turns_reflection_network_power_and_replay() -> None:
    result = ideal_transformer_identity_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_rejects_current_gain_mismatch_even_when_voltage_ratio_survives() -> None:
    bad = copy.deepcopy(_summary())
    positive = bad["metrics"]["positive"]
    positive["primary_current_phasor_rms_a"][0] *= 0.8
    positive["primary_current_relative_error"] = 0.2
    positive["current_turns_identity_relative_error"] = 0.25
    positive["reflected_load_relative_error"] = 0.25
    positive["transformer_complex_power_relative_error"] = 0.2
    positive["source_power_balance_relative_error"] = 0.2
    positive["instantaneous_transformer_power_relative_error"] = 0.2
    positive["instantaneous_source_balance_relative_error"] = 0.25
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["voltage_and_current_turns_identities_close"] is False
    assert result["checks"]["load_reflection_and_series_source_network_close"] is False
    assert result["checks"]["complex_and_instantaneous_power_are_conserved"] is False


def test_rejects_stale_reported_error_and_short_fit_window() -> None:
    bad = copy.deepcopy(_summary())
    bad["metrics"]["positive"]["voltage_turns_identity_relative_error"] = 0.1
    bad["metrics"]["positive"]["fit_point_count"] = 20
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["reported_errors_and_complex_powers_match_recomputation"] is False
    assert result["checks"]["two_period_phasor_fit_is_dense_and_accurate"] is False


def test_rejects_non_mapping_input() -> None:
    try:
        ideal_transformer_identity_gate([])  # type: ignore[arg-type]
    except ValueError as exc:
        assert "object" in str(exc)
    else:
        raise AssertionError("non-mapping input was accepted")


def test_rejects_source_phasor_breaking_kvl_and_power_balance() -> None:
    bad = copy.deepcopy(_summary())
    bad["metrics"]["positive"]["source_voltage_phasor_rms_v"][0] *= 1.15
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["load_reflection_and_series_source_network_close"] is False
    assert result["checks"]["complex_and_instantaneous_power_are_conserved"] is False


@pytest.mark.parametrize(
    "case_id",
    ["topology", "reflected_load", "fit_density", "replay_error", "primary_current"],
)
def test_counterfactual_curriculum90_public(case_id: str) -> None:
    bad = copy.deepcopy(_summary())
    if case_id == "topology":
        bad["model_contract"]["topology"] = "unknown"
    elif case_id == "reflected_load":
        bad["model_contract"]["expected_reflected_load_ohm"] *= 1.1
    elif case_id == "fit_density":
        bad["metrics"]["positive"]["fit_point_count"] = 20
    elif case_id == "replay_error":
        bad["metrics"]["maximum_phasor_replay_relative_error"] = 0.1
    else:
        bad["metrics"]["positive"]["primary_current_phasor_rms_a"][0] *= 0.8
    assert ideal_transformer_identity_gate(bad)["status"] == "needs_attention"


def test_generalization_v3s_rejects_secondary_voltage_drift() -> None:
    bad = copy.deepcopy(_summary())
    bad["metrics"]["positive"]["secondary_voltage_phasor_rms_v"][0] *= 1.2
    assert ideal_transformer_identity_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v4_turns_ratio", "v4_source_offset", "v4_frequency", "v4_load_resistance", "v4_timing_stage"],
)
def test_counterfactual_curriculum90_v4_public(case_id: str) -> None:
    bad = copy.deepcopy(_summary())
    if case_id == "v4_turns_ratio":
        bad["model_contract"]["turns_ratio_primary_to_secondary"] = 0.0
    elif case_id == "v4_source_offset":
        bad["model_contract"]["source_offset_v"] = 1.0
    elif case_id == "v4_frequency":
        bad["model_contract"]["frequency_hz"] = 0.0
    elif case_id == "v4_load_resistance":
        bad["model_contract"]["load_resistance_ohm"] = -100.0
    else:
        bad["timing_breakdown_s"].pop("serialization")
    try:
        result = ideal_transformer_identity_gate(bad)
    except ValueError:
        return
    assert result["status"] == "needs_attention"


def test_generalization_v5_rejects_negative_series_resistance() -> None:
    bad = copy.deepcopy(_summary())
    bad["model_contract"]["series_resistance_ohm"] = -1.0
    with pytest.raises(ValueError):
        ideal_transformer_identity_gate(bad)
