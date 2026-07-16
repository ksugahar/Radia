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
        "phasor_fit_frequency_hz": 60.0,
        "phasor_replay_evidence": [
            {
                "fit_id": "positive-replay-1",
                "raw_generation_id": "raw-generation-41",
                "trace_group_generation_id": "raw-generation-41",
                "scalar_phasor_digest": "a" * 64,
            },
            {
                "fit_id": "positive-replay-2",
                "raw_generation_id": "raw-generation-42",
                "trace_group_generation_id": "raw-generation-42",
                "scalar_phasor_digest": "a" * 64,
            },
        ],
        "fit_window_segment_contract": {
            "fit_run_generation_id": "transient-generation-42",
            "sample_run_generation_ids": ["transient-generation-42"],
            "restart_discontinuities_s": [],
        },
        "current_role_contract": {
            "source_delivery_current_phasor_rms_a": "source_delivery_into_network",
            "primary_current_phasor_rms_a": "transformer_primary_absorption",
            "secondary_current_phasor_rms_a": "transformer_secondary_delivery_to_load",
        },
        "phasor_current_roles": {
            "source_delivery_current_phasor_rms_a": "source_delivery_into_network",
            "primary_current_phasor_rms_a": "transformer_primary_absorption",
            "secondary_current_phasor_rms_a": "transformer_secondary_delivery_to_load",
        },
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


@pytest.mark.parametrize(
    "case_id",
    ["v6_public_complex_power_drift", "v6_public_invalid_fit_window"],
)
def test_generalization_v6_public(case_id: str) -> None:
    bad = copy.deepcopy(_summary())
    if case_id == "v6_public_complex_power_drift":
        bad["metrics"]["positive"]["transformer_input_complex_power_va"][0] *= 1.10
    else:
        bad["metrics"]["positive"]["fit_window_start_s"] = bad["metrics"][
            "positive"
        ]["fit_window_stop_s"]
    assert ideal_transformer_identity_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_phasor_source_frequency_mismatch",
        "v7_public_current_role_power_false_closure",
    ],
)
def test_generalization_v7_public(case_id: str) -> None:
    bad = copy.deepcopy(_summary())
    positive = bad["metrics"]["positive"]
    if case_id == "v7_public_phasor_source_frequency_mismatch":
        positive["phasor_fit_frequency_hz"] = 61.0
    else:
        roles = positive["phasor_current_roles"]
        roles["source_delivery_current_phasor_rms_a"], roles[
            "secondary_current_phasor_rms_a"
        ] = (
            roles["secondary_current_phasor_rms_a"],
            roles["source_delivery_current_phasor_rms_a"],
        )
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    expected_check = (
        "phasor_fit_frequency_matches_source_contract"
        if case_id == "v7_public_phasor_source_frequency_mismatch"
        else "current_phasor_roles_match_terminal_contract"
    )
    assert result["checks"][expected_check] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v8_public_phasor_replay_mixed_raw_generation",
        "v8_public_fit_window_crosses_restart_discontinuity",
    ],
)
def test_generalization_v8_public(case_id: str) -> None:
    bad = copy.deepcopy(_summary())
    positive = bad["metrics"]["positive"]
    if case_id == "v8_public_phasor_replay_mixed_raw_generation":
        positive["phasor_replay_evidence"][1][
            "trace_group_generation_id"
        ] = "raw-generation-41"
        expected_check = "phasor_replays_bind_each_fit_to_one_raw_generation"
    else:
        positive["fit_window_segment_contract"].update(
            {
                "sample_run_generation_ids": [
                    "transient-generation-41",
                    "transient-generation-42",
                ],
                "restart_discontinuities_s": [0.025],
            }
        )
        expected_check = "phasor_fit_window_does_not_cross_a_restart_segment"
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected_check] is False


def _with_v9_phasor_identity(summary: dict) -> dict:
    positive = summary["metrics"]["positive"]
    names = (
        "source_voltage_phasor_rms_v",
        "primary_voltage_phasor_rms_v",
        "secondary_voltage_phasor_rms_v",
        "source_delivery_current_phasor_rms_a",
        "primary_current_phasor_rms_a",
        "secondary_current_phasor_rms_a",
    )
    positive["phasor_basis_contract"] = {
        "quantity_basis": {name: "rms" for name in names},
        "normalization_factor_to_rms": {name: 1.0 for name in names},
        "complex_power_formula": "rms_voltage_times_conjugate_rms_current",
    }
    positive["phase_unwrap_contract"] = {
        "phase_unit": "radian",
        "branch_period_rad": 2.0 * math.pi,
        "unwrap_convention": "continuous_signed_phase",
        "reference_trace": "source_voltage_phasor_rms_v",
        "fit_branch_sign_digest": "b" * 64,
        "replay_branch_sign_digest": "b" * 64,
    }
    return summary


def _with_v10_power_and_frequency_identity(summary: dict) -> dict:
    summary = _with_v9_phasor_identity(summary)
    positive = summary["metrics"]["positive"]
    positive["power_sign_convention"] = {
        "source_power_role": "delivered_positive",
        "passive_power_role": "absorbed_positive",
        "source_current_reference": "leaving_positive_terminal",
        "passive_current_reference": "entering_positive_terminal",
        "balance_equation": "source_delivered_equals_passive_absorbed",
        "recorded_sign_transform": "none",
    }
    positive["ac_frequency_interpolation_contract"] = {
        "source_frequency_coordinate": "log10_hz",
        "target_frequency_coordinate": "log10_hz",
        "interpolation_coordinate": "log10_hz",
        "frequency_unit": "Hz",
        "source_grid_generation": "ac-grid-42",
        "interpolated_trace_source_grid_generation": "ac-grid-42",
    }
    return summary


def _with_v11_transient_and_noise_identity(summary: dict) -> dict:
    summary = _with_v10_power_and_frequency_identity(summary)
    positive = summary["metrics"]["positive"]
    positive["transient_energy_window_event_identity"] = {
        "event_type": "switch_turn_on",
        "energy_window_start_event_index": 10,
        "reference_window_start_event_index": 10,
        "energy_window_duration_s": 2.0e-6,
        "reference_window_duration_s": 2.0e-6,
        "event_detection_generation": "switch-events-43",
        "energy_window_event_generation": "switch-events-43",
        "reference_window_event_generation": "switch-events-43",
    }
    positive["noise_density_band_identity"] = {
        "noise_density_unit": "V/sqrt(Hz)",
        "numerical_frequency_unit": "Hz",
        "band_limit_unit": "Hz",
        "numerical_frequency_scale_to_hz": 1.0,
        "band_limit_scale_to_hz": 1.0,
        "integrated_noise_unit": "V_rms",
        "frequency_grid_generation": "noise-grid-43",
        "band_integration_grid_generation": "noise-grid-43",
    }
    return summary


def _with_v12_cycle_and_seed_identity(summary: dict) -> dict:
    summary = _with_v11_transient_and_noise_identity(summary)
    positive = summary["metrics"]["positive"]
    positive["steady_cycle_average_identity"] = {
        "settled_cycle_generation_id": "steady-cycles-44",
        "waveform_cycle_generation_id": "steady-cycles-44",
        "average_window_cycle_generation_id": "steady-cycles-44",
        "first_settled_cycle_index": 18,
        "average_window_cycle_index": 20,
        "average_window_period_count": 1,
    }
    positive["monte_carlo_parameter_seed_identity"] = {
        "parameter_order": ["RPRI", "RLOAD", "CSTRAY"],
        "statistics_parameter_order": ["RPRI", "RLOAD", "CSTRAY"],
        "seed_policy": "one_seed_per_named_parameter",
        "seed_schedule_generation_id": "mc-seeds-44",
        "statistics_seed_schedule_generation_id": "mc-seeds-44",
        "parameter_seed_map_sha256": "d" * 64,
        "statistics_parameter_seed_map_sha256": "d" * 64,
    }
    return summary


def _with_v13_phase_and_derivative_identity(summary: dict) -> dict:
    summary = _with_v12_cycle_and_seed_identity(summary)
    positive = summary["metrics"]["positive"]
    positive["ac_phase_coordinate_identity"] = {
        "ac_sweep_generation_id": "ac-sweep-45",
        "phase_trace_sweep_generation_id": "ac-sweep-45",
        "reference_phase_trace_sweep_generation_id": "ac-sweep-45",
        "phase_coordinate_unit": "degree",
        "reference_phase_coordinate_unit": "degree",
        "phase_coordinate_convention": "principal_degree_minus180_180",
        "reference_phase_coordinate_convention": (
            "principal_degree_minus180_180"
        ),
        "phase_trace_value_sha256": "a" * 64,
        "reference_phase_trace_value_sha256": "a" * 64,
    }
    positive["transient_derivative_adaptive_history_identity"] = {
        "transient_generation_id": "transient-45",
        "accepted_step_generation_id": "transient-45",
        "derivative_history_generation_id": "transient-45",
        "current_sample_generation_id": "transient-45",
        "derivative_sample_generation_id": "transient-45",
        "accepted_time_grid_sha256": "b" * 64,
        "derivative_history_time_grid_sha256": "b" * 64,
        "accepted_step_count": 257,
        "derivative_history_step_count": 257,
        "derivative_scheme": "variable_step_bdf2",
        "history_derivative_scheme": "variable_step_bdf2",
    }
    return summary


def test_v9_public_peak_rms_phasor_basis_mixed() -> None:
    bad = _with_v9_phasor_identity(_summary())
    contract = bad["metrics"]["positive"]["phasor_basis_contract"]
    contract["quantity_basis"]["source_voltage_phasor_rms_v"] = "peak"
    contract["normalization_factor_to_rms"][
        "source_voltage_phasor_rms_v"
    ] = 2.0**-0.5
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_phasors_share_one_rms_normalization_basis"] is False


def test_v9_public_phase_unwrap_branch_sign_alias() -> None:
    bad = _with_v9_phasor_identity(_summary())
    bad["metrics"]["positive"]["phase_unwrap_contract"][
        "replay_branch_sign_digest"
    ] = "c" * 64
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["phase_replay_preserves_unwrap_branch_orientation"]
        is False
    )


def test_v10_public_source_power_passive_sign_mismatch() -> None:
    bad = _with_v10_power_and_frequency_identity(_summary())
    bad["metrics"]["positive"]["power_sign_convention"][
        "source_current_reference"
    ] = "entering_positive_terminal"
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["source_and_passive_power_use_recorded_sign_conventions"]
        is False
    )


def test_v10_public_ac_frequency_interpolation_scale_mismatch() -> None:
    bad = _with_v10_power_and_frequency_identity(_summary())
    bad["metrics"]["positive"]["ac_frequency_interpolation_contract"].update(
        {
            "source_frequency_coordinate": "linear_hz",
            "interpolated_trace_source_grid_generation": "ac-grid-41",
        }
    )
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["ac_traces_share_frequency_interpolation_coordinate"]
        is False
    )


def test_v11_public_transient_energy_window_event_alignment_mismatch() -> None:
    bad = _with_v11_transient_and_noise_identity(_summary())
    bad["metrics"]["positive"]["transient_energy_window_event_identity"].update(
        {
            "reference_window_start_event_index": 11,
            "reference_window_event_generation": "switch-events-42",
        }
    )
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["transient_energy_windows_share_switching_event_phase"]
        is False
    )


def test_v11_public_noise_density_integrated_band_unit_mismatch() -> None:
    bad = _with_v11_transient_and_noise_identity(_summary())
    bad["metrics"]["positive"]["noise_density_band_identity"].update(
        {"band_limit_unit": "kHz", "band_limit_scale_to_hz": 1.0}
    )
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "noise_density_and_band_limits_share_frequency_units"
        ]
        is False
    )


def test_v12_public_switched_average_steady_cycle_generation_mismatch() -> None:
    bad = _with_v12_cycle_and_seed_identity(_summary())
    bad["metrics"]["positive"]["steady_cycle_average_identity"].update(
        {
            "average_window_cycle_generation_id": "steady-cycles-43",
            "average_window_cycle_index": 12,
        }
    )
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["switched_averages_use_a_settled_cycle_generation"]
        is False
    )


def test_v12_public_monte_carlo_seed_parameter_order_mismatch() -> None:
    bad = _with_v12_cycle_and_seed_identity(_summary())
    bad["metrics"]["positive"]["monte_carlo_parameter_seed_identity"][
        "statistics_parameter_order"
    ] = ["CSTRAY", "RLOAD", "RPRI"]
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "monte_carlo_statistics_preserve_named_parameter_seed_order"
        ]
        is False
    )


def test_v13_public_ac_phase_degree_radian_coordinate_mismatch() -> None:
    bad = _with_v13_phase_and_derivative_identity(_summary())
    bad["metrics"]["positive"]["ac_phase_coordinate_identity"].update(
        {
            "reference_phase_coordinate_unit": "radian",
            "reference_phase_coordinate_convention": (
                "principal_radian_minuspi_pi"
            ),
        }
    )
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["ac_phase_traces_share_coordinate_unit_and_convention"]
        is False
    )


def test_v13_public_transient_derivative_adaptive_timestep_history_mismatch() -> None:
    bad = _with_v13_phase_and_derivative_identity(_summary())
    bad["metrics"]["positive"][
        "transient_derivative_adaptive_history_identity"
    ].update(
        {
            "derivative_history_generation_id": "transient-44",
            "derivative_history_time_grid_sha256": "c" * 64,
            "derivative_history_step_count": 241,
        }
    )
    result = ideal_transformer_identity_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "transient_derivatives_use_current_accepted_timestep_history"
        ]
        is False
    )
