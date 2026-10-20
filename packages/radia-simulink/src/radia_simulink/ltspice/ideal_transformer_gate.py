"""Stable public entry point for transformer and attached evidence checks.

Implementation ownership is split across private common, signal, circuit,
dynamic, and topology modules. Keep schema keys, tolerances, and public names
stable here; add new checks to the owning module rather than this facade.
The numbered ltspice_v*_gates interfaces remain compatibility contracts, not
LTspice product-version requirements. Retire them only with caller migration.
"""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from .ltspice_v44_gates import validate_ltspice_v44_public_identity
from .ltspice_v45_gates import validate_ltspice_v45_identity
from .ltspice_v46_gates import validate_ltspice_v46_identity
from .ltspice_v47_gates import validate_ltspice_v47_identity
from .ltspice_v48_gates import validate_ltspice_v48_identity
from .ltspice_v49_gates import validate_ltspice_v49_identity
from .ltspice_v50_gates import validate_ltspice_v50_identity
from typing import Any
from ._identity_common import (
    _complex_pair,
    _finite,
    _is_sha256,
    _relative_error,
)
from ._identity_signal import (
    _ac_frequency_interpolation_contract_ok,
    _ac_group_delay_phase_unwrap_grid_identity_ok,
    _ac_noise_integrated_density_bin_identity_ok,
    _ac_phase_coordinate_identity_ok,
    _fft_window_coherent_gain_amplitude_basis_identity_ok,
    _fft_window_harmonic_bin_identity_ok,
    _fit_window_stays_in_one_segment,
    _fourier_phase_reference_time_origin_identity_ok,
    _measure_crossing_interpolation_grid_identity_ok,
    _measure_trigger_target_crossing_identity_ok,
    _monte_carlo_parameter_seed_identity_ok,
    _monte_carlo_percentile_sample_filter_identity_ok,
    _monte_carlo_sample_trace_identity_ok,
    _noise_density_band_identity_ok,
    _noise_monte_carlo_psd_integration_identity_ok,
    _noise_spectral_density_sidedness_identity_ok,
    _phase_unwrap_contract_ok,
    _phasor_basis_contract_ok,
    _phasor_replay_generations_ok,
    _power_sign_convention_ok,
    _steady_cycle_average_identity_ok,
    _stepped_ac_parameter_tuple_grid_identity_ok,
    _stepped_parameter_interpolation_coordinate_identity_ok,
    _stepped_transient_measure_row_identity_ok,
    _switched_converter_cycle_measure_identity_ok,
    _transient_derivative_adaptive_history_identity_ok,
    _transient_energy_window_event_identity_ok,
    _transient_power_interpolation_grid_identity_ok,
    _transient_rms_average_event_window_identity_ok,
)
from ._identity_circuit import (
    _ac_noise_source_identity_ok,
    _ac_sweep_measure_identity_ok,
    _behavioral_switch_hysteresis_identity_ok,
    _electrothermal_generation_identity_ok,
    _electrothermal_waveform_closure_identity_ok,
    _hierarchical_step_identity_ok,
    _loop_gain_owner_identity_ok,
    _monte_carlo_subcircuit_identity_ok,
    _monte_carlo_yield_owner_identity_ok,
    _mosfet_soa_owner_identity_ok,
    _mosfet_switching_loss_owner_identity_ok,
    _noise_integration_owner_identity_ok,
    _smps_efficiency_owner_identity_ok,
    _step_response_owner_identity_ok,
    _stepped_monte_carlo_aggregation_identity_ok,
    _switch_event_timing_owner_identity_ok,
    _transient_startup_identity_ok,
)
from ._identity_dynamic import (
    _behavioral_source_event_owner_identity_ok,
    _conducted_emi_owner_identity_ok,
    _correlation_matrix_is_positive_semidefinite,
    _electrothermal_fixedpoint_owner_identity_ok,
    _noise_correlation_band_owner_identity_ok,
    _noise_referral_owner_identity_ok,
    _oscillator_closure_owner_identity_ok,
    _periodic_steady_state_owner_identity_ok,
    _pll_closure_owner_identity_ok,
    _sampled_loop_owner_identity_ok,
    _small_signal_consistency_owner_identity_ok,
    _smps_startup_owner_identity_ok,
    _switched_limit_cycle_owner_identity_ok,
    _touchstone_network_owner_identity_ok,
    _transmission_line_transient_owner_identity_ok,
)
from ._identity_topology import (
    _active_filter_owner_identity_ok,
    _bjt_amplifier_owner_identity_ok,
    _boost_converter_owner_identity_ok,
    _bridge_rectifier_owner_identity_ok,
    _flyback_owner_identity_ok,
    _flyback_v43_identity_ok,
    _instrumentation_amplifier_owner_identity_ok,
    _llc_resonant_owner_identity_ok,
    _mosfet_gatedrive_owner_identity_ok,
    _opamp_stability_v43_identity_ok,
    _sallen_key_owner_identity_ok,
    _transimpedance_contract_ok,
    _transimpedance_owner_identity_ok,
)
from .ltspice_v43_gates import flyback_identity_ok, opamp_stability_identity_ok


def ideal_transformer_identity_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate turns ratio, reflected impedance, network closure, power, and replay."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    model = summary.get("model_contract")
    metrics = summary.get("metrics")
    timing = summary.get("timing_breakdown_s")
    if not isinstance(model, Mapping) or not isinstance(metrics, Mapping):
        raise ValueError("model_contract and metrics must be objects")
    positive = metrics.get("positive")
    if not isinstance(positive, Mapping):
        raise ValueError("metrics.positive must be an object")

    turns_ratio = _finite(
        model.get("turns_ratio_primary_to_secondary"),
        "turns_ratio_primary_to_secondary",
        positive=True,
    )
    source_offset = _finite(model.get("source_offset_v"), "source_offset_v")
    source_peak = _finite(model.get("source_peak_v"), "source_peak_v", positive=True)
    frequency = _finite(model.get("frequency_hz"), "frequency_hz", positive=True)
    series_resistance = _finite(
        model.get("series_resistance_ohm"), "series_resistance_ohm", positive=True
    )
    load_resistance = _finite(
        model.get("load_resistance_ohm"), "load_resistance_ohm", positive=True
    )

    source_rms = source_peak / math.sqrt(2.0)
    reflected_load = turns_ratio**2 * load_resistance
    primary_current_expected = source_rms / (series_resistance + reflected_load)
    primary_voltage_expected = primary_current_expected * reflected_load
    secondary_voltage_expected = primary_voltage_expected / turns_ratio
    secondary_current_expected = secondary_voltage_expected / load_resistance
    load_power_expected = secondary_voltage_expected**2 / load_resistance
    series_loss_expected = primary_current_expected**2 * series_resistance
    derived = {
        "expected_source_voltage_rms_v": source_rms,
        "expected_reflected_load_ohm": reflected_load,
        "expected_primary_voltage_rms_v": primary_voltage_expected,
        "expected_secondary_voltage_rms_v": secondary_voltage_expected,
        "expected_primary_current_rms_a": primary_current_expected,
        "expected_secondary_current_rms_a": secondary_current_expected,
        "expected_load_power_w": load_power_expected,
        "expected_series_loss_w": series_loss_expected,
    }
    declared_expectation_error = max(
        _relative_error(_finite(model.get(name), name, positive=True), value)
        for name, value in derived.items()
    )

    v_source = _complex_pair(
        positive.get("source_voltage_phasor_rms_v"),
        "source_voltage_phasor_rms_v",
    )
    v_primary = _complex_pair(
        positive.get("primary_voltage_phasor_rms_v"),
        "primary_voltage_phasor_rms_v",
    )
    v_secondary = _complex_pair(
        positive.get("secondary_voltage_phasor_rms_v"),
        "secondary_voltage_phasor_rms_v",
    )
    i_source = _complex_pair(
        positive.get("source_delivery_current_phasor_rms_a"),
        "source_delivery_current_phasor_rms_a",
    )
    i_primary = _complex_pair(
        positive.get("primary_current_phasor_rms_a"),
        "primary_current_phasor_rms_a",
    )
    i_secondary = _complex_pair(
        positive.get("secondary_current_phasor_rms_a"),
        "secondary_current_phasor_rms_a",
    )
    if min(abs(v_source), abs(v_primary), abs(v_secondary), abs(i_primary), abs(i_secondary)) <= 1.0e-300:
        raise ValueError("transformer phasors must be nonzero")

    transformer_input_power = v_primary * i_primary.conjugate()
    transformer_output_power = v_secondary * i_secondary.conjugate()
    source_delivery_power = v_source * i_source.conjugate()
    series_loss = series_resistance * abs(i_primary) ** 2
    recomputed = {
        "source_voltage_relative_error": _relative_error(abs(v_source), source_rms),
        "primary_voltage_relative_error": _relative_error(
            abs(v_primary), primary_voltage_expected
        ),
        "secondary_voltage_relative_error": _relative_error(
            abs(v_secondary), secondary_voltage_expected
        ),
        "primary_current_relative_error": _relative_error(
            abs(i_primary), primary_current_expected
        ),
        "secondary_current_relative_error": _relative_error(
            abs(i_secondary), secondary_current_expected
        ),
        "voltage_turns_identity_relative_error": _relative_error(
            v_secondary, -v_primary / turns_ratio
        ),
        "current_turns_identity_relative_error": _relative_error(
            i_secondary, -turns_ratio * i_primary
        ),
        "secondary_load_ohm_law_relative_error": _relative_error(
            v_secondary, load_resistance * i_secondary
        ),
        "reflected_load_relative_error": _relative_error(
            v_primary / i_primary, reflected_load
        ),
        "source_series_kvl_relative_error": _relative_error(
            v_source, v_primary + series_resistance * i_primary
        ),
        "source_primary_current_relative_error": _relative_error(i_source, i_primary),
        "transformer_complex_power_relative_error": _relative_error(
            transformer_input_power, transformer_output_power
        ),
        "source_power_balance_relative_error": _relative_error(
            source_delivery_power, transformer_output_power + series_loss
        ),
    }
    reported_error_drift = max(
        abs(_finite(positive.get(name), name) - value)
        for name, value in recomputed.items()
    )

    reported_powers = {
        "transformer_input_complex_power_va": _complex_pair(
            positive.get("transformer_input_complex_power_va"),
            "transformer_input_complex_power_va",
        ),
        "transformer_output_complex_power_va": _complex_pair(
            positive.get("transformer_output_complex_power_va"),
            "transformer_output_complex_power_va",
        ),
        "source_delivery_complex_power_va": _complex_pair(
            positive.get("source_delivery_complex_power_va"),
            "source_delivery_complex_power_va",
        ),
    }
    power_report_drift = max(
        _relative_error(
            reported_powers["transformer_input_complex_power_va"],
            transformer_input_power,
        ),
        _relative_error(
            reported_powers["transformer_output_complex_power_va"],
            transformer_output_power,
        ),
        _relative_error(
            reported_powers["source_delivery_complex_power_va"],
            source_delivery_power,
        ),
        _relative_error(_finite(positive.get("series_loss_w"), "series_loss_w"), series_loss),
    )

    instantaneous_transformer_error = _finite(
        positive.get("instantaneous_transformer_power_relative_error"),
        "instantaneous_transformer_power_relative_error",
    )
    instantaneous_source_error = _finite(
        positive.get("instantaneous_source_balance_relative_error"),
        "instantaneous_source_balance_relative_error",
    )
    fit_residual = _finite(
        positive.get("maximum_phasor_fit_relative_residual"),
        "maximum_phasor_fit_relative_residual",
    )
    fit_frequency_matches_source = True
    if "phasor_fit_frequency_hz" in positive:
        fit_frequency = _finite(
            positive.get("phasor_fit_frequency_hz"),
            "phasor_fit_frequency_hz",
            positive=True,
        )
        fit_frequency_matches_source = (
            abs(fit_frequency - frequency) <= 1.0e-12 * max(frequency, 1.0)
        )

    expected_current_roles = {
        "source_delivery_current_phasor_rms_a": "source_delivery_into_network",
        "primary_current_phasor_rms_a": "transformer_primary_absorption",
        "secondary_current_phasor_rms_a": "transformer_secondary_delivery_to_load",
    }
    role_contract = positive.get("current_role_contract")
    observed_roles = positive.get("phasor_current_roles")
    role_evidence_present = role_contract is not None or observed_roles is not None
    current_roles_match_contract = not role_evidence_present or (
        isinstance(role_contract, Mapping)
        and isinstance(observed_roles, Mapping)
        and dict(role_contract) == expected_current_roles
        and dict(observed_roles) == expected_current_roles
    )
    point_count = int(_finite(positive.get("point_count"), "point_count", positive=True))
    fit_point_count = int(
        _finite(positive.get("fit_point_count"), "fit_point_count", positive=True)
    )
    fit_start = _finite(positive.get("fit_window_start_s"), "fit_window_start_s")
    fit_stop = _finite(
        positive.get("fit_window_stop_s"), "fit_window_stop_s", positive=True
    )
    replay_error = _finite(
        metrics.get("maximum_phasor_replay_relative_error"),
        "maximum_phasor_replay_relative_error",
    )

    timing_ok = False
    if isinstance(timing, Mapping) and len(timing) == 4:
        try:
            timing_ok = all(_finite(value, "timing") >= 0.0 for value in timing.values())
        except ValueError:
            timing_ok = False

    checks = {
        "ideal_two_winding_transformer_contract": model.get("topology")
        == "ideal_two_winding_transformer_with_series_source_resistance"
        and abs(source_offset) <= 1.0e-15,
        "declared_analytic_expectations_are_recomputed": declared_expectation_error
        <= 1.0e-12,
        "source_and_analytic_rms_values_close": max(
            recomputed["source_voltage_relative_error"],
            recomputed["primary_voltage_relative_error"],
            recomputed["secondary_voltage_relative_error"],
            recomputed["primary_current_relative_error"],
            recomputed["secondary_current_relative_error"],
        )
        <= 2.0e-5,
        "voltage_and_current_turns_identities_close": max(
            recomputed["voltage_turns_identity_relative_error"],
            recomputed["current_turns_identity_relative_error"],
        )
        <= 2.0e-6,
        "load_reflection_and_series_source_network_close": max(
            recomputed["secondary_load_ohm_law_relative_error"],
            recomputed["reflected_load_relative_error"],
            recomputed["source_series_kvl_relative_error"],
            recomputed["source_primary_current_relative_error"],
        )
        <= 2.0e-6,
        "complex_and_instantaneous_power_are_conserved": max(
            recomputed["transformer_complex_power_relative_error"],
            recomputed["source_power_balance_relative_error"],
            instantaneous_transformer_error,
            instantaneous_source_error,
        )
        <= 2.0e-6,
        "reported_errors_and_complex_powers_match_recomputation": max(
            reported_error_drift, power_report_drift
        )
        <= 1.0e-12,
        "two_period_phasor_fit_is_dense_and_accurate": point_count >= 1000
        and fit_point_count >= 500
        and fit_stop - fit_start >= 1.9 / frequency
        and fit_residual <= 2.0e-5,
        "phasor_fit_frequency_matches_source_contract": fit_frequency_matches_source,
        "current_phasor_roles_match_terminal_contract": current_roles_match_contract,
        "positive_phasor_replay_is_deterministic": replay_error <= 1.0e-12,
        "phasor_replays_bind_each_fit_to_one_raw_generation": (
            _phasor_replay_generations_ok(positive)
        ),
        "phasor_fit_window_does_not_cross_a_restart_segment": (
            _fit_window_stays_in_one_segment(positive, fit_start, fit_stop)
        ),
        "all_phasors_share_one_rms_normalization_basis": (
            _phasor_basis_contract_ok(positive)
        ),
        "phase_replay_preserves_unwrap_branch_orientation": (
            _phase_unwrap_contract_ok(positive)
        ),
        "source_and_passive_power_use_recorded_sign_conventions": (
            _power_sign_convention_ok(positive)
        ),
        "ac_traces_share_frequency_interpolation_coordinate": (
            _ac_frequency_interpolation_contract_ok(positive)
        ),
        "transient_energy_windows_share_switching_event_phase": (
            _transient_energy_window_event_identity_ok(positive)
        ),
        "noise_density_and_band_limits_share_frequency_units": (
            _noise_density_band_identity_ok(positive)
        ),
        "switched_averages_use_a_settled_cycle_generation": (
            _steady_cycle_average_identity_ok(positive)
        ),
        "monte_carlo_statistics_preserve_named_parameter_seed_order": (
            _monte_carlo_parameter_seed_identity_ok(positive)
        ),
        "ac_phase_traces_share_coordinate_unit_and_convention": (
            _ac_phase_coordinate_identity_ok(positive)
        ),
        "transient_derivatives_use_current_accepted_timestep_history": (
            _transient_derivative_adaptive_history_identity_ok(positive)
        ),
        "noise_density_sidedness_uses_required_amplitude_conversion": (
            _noise_spectral_density_sidedness_identity_ok(positive)
        ),
        "stepped_results_share_parameter_interpolation_coordinate": (
            _stepped_parameter_interpolation_coordinate_identity_ok(positive)
        ),
        "fft_amplitudes_use_current_window_gain_and_basis": (
            _fft_window_coherent_gain_amplitude_basis_identity_ok(positive)
        ),
        "monte_carlo_percentiles_use_current_filtered_sample_set": (
            _monte_carlo_percentile_sample_filter_identity_ok(positive)
        ),
        "measure_crossings_use_current_accepted_step_grid": (
            _measure_crossing_interpolation_grid_identity_ok(positive)
        ),
        "fourier_harmonic_phases_share_reference_time_origin": (
            _fourier_phase_reference_time_origin_identity_ok(positive)
        ),
        "ac_group_delay_uses_current_phase_unwrap_frequency_grid": (
            _ac_group_delay_phase_unwrap_grid_identity_ok(positive)
        ),
        "transient_rms_and_average_share_switching_event_window": (
            _transient_rms_average_event_window_identity_ok(positive)
        ),
        "ac_noise_integration_uses_current_sidedness_frequency_bins": (
            _ac_noise_integrated_density_bin_identity_ok(positive)
        ),
        "transient_power_uses_one_current_voltage_interpolation_grid": (
            _transient_power_interpolation_grid_identity_ok(positive)
        ),
        "stepped_ac_traces_use_current_parameter_tuple_and_frequency_grid": (
            _stepped_ac_parameter_tuple_grid_identity_ok(positive)
        ),
        "measure_trigger_target_use_current_crossing_counts_and_brackets": (
            _measure_trigger_target_crossing_identity_ok(positive)
        ),
        "monte_carlo_traces_use_current_seed_sample_tuple_row_order": (
            _monte_carlo_sample_trace_identity_ok(positive)
        ),
        "fft_harmonics_use_current_window_sample_rate_and_bin_table": (
            _fft_window_harmonic_bin_identity_ok(positive)
        ),
        "noise_monte_carlo_psd_integration_uses_current_filtered_samples_bins_and_sidedness": (
            _noise_monte_carlo_psd_integration_identity_ok(positive)
        ),
        "stepped_transient_measures_use_current_grid_parameter_tuples_and_row_order": (
            _stepped_transient_measure_row_identity_ok(positive)
        ),
        "switched_converter_measures_use_current_initial_state_topology_grid_and_cycle_window": (
            _switched_converter_cycle_measure_identity_ok(positive)
        ),
        "electrothermal_results_use_current_temperature_model_loss_network_and_timestep": (
            _electrothermal_generation_identity_ok(positive)
        ),
        "monte_carlo_statistics_use_current_seeds_models_parameters_and_raw_samples": (
            _monte_carlo_subcircuit_identity_ok(positive)
        ),
        "behavioral_switch_measures_use_current_hysteresis_events_timesteps_and_window": (
            _behavioral_switch_hysteresis_identity_ok(positive)
        ),
        "hierarchical_steps_use_current_scope_model_bin_temperature_and_sample_rows": (
            _hierarchical_step_identity_ok(positive)
        ),
        "ac_noise_uses_current_source_normalization_aliases_complex_axis_and_grid": (
            _ac_noise_source_identity_ok(positive)
        ),
        "transient_startup_uses_current_initial_conditions_uic_operating_point_grid_and_waveforms": (
            _transient_startup_identity_ok(positive)
        ),
        "stepped_monte_carlo_aggregation_uses_current_rows_seeds_filters_weights_and_rule": (
            _stepped_monte_carlo_aggregation_identity_ok(positive)
        ),
        "ac_measures_use_current_sweep_grid_complex_basis_phase_unwrap_rows_and_result": (
            _ac_sweep_measure_identity_ok(positive)
        ),
        "electrothermal_waveforms_use_current_power_temperature_model_network_timestep_and_result": (
            _electrothermal_waveform_closure_identity_ok(positive)
        ),
        "noise_integration_uses_current_sources_normalization_sidedness_grid_units_and_result": (
            _noise_integration_owner_identity_ok(positive)
        ),
        "switch_timing_uses_current_hysteresis_events_timestep_window_waveform_and_measures": (
            _switch_event_timing_owner_identity_ok(positive)
        ),
        "smps_efficiency_uses_current_traces_cycle_window_integration_waveform_timestep_and_result": (
            _smps_efficiency_owner_identity_ok(positive)
        ),
        "loop_gain_uses_current_break_injection_phase_grid_crossover_margins_and_result": (
            _loop_gain_owner_identity_ok(positive)
        ),
        "mosfet_switching_loss_uses_current_gate_charge_overlap_deadtime_events_temperature_cycle_and_result": (
            _mosfet_switching_loss_owner_identity_ok(positive)
        ),
        "step_response_uses_current_initial_final_rise_settling_overshoot_window_and_waveform": (
            _step_response_owner_identity_ok(positive)
        ),
        "mosfet_soa_uses_current_voltage_current_pulse_duty_temperature_model_waveform_and_result": (
            _mosfet_soa_owner_identity_ok(positive)
        ),
        "monte_carlo_yield_uses_current_distributions_tolerances_seeds_failure_samples_owner_and_result": (
            _monte_carlo_yield_owner_identity_ok(positive)
        ),
        "behavioral_sources_use_current_event_grid_derivative_charge_energy_initial_owner_and_result": (
            _behavioral_source_event_owner_identity_ok(positive)
        ),
        "touchstone_networks_use_current_impedance_units_parameters_ports_complex_passivity_file_and_result": (
            _touchstone_network_owner_identity_ok(positive)
        ),
        "smps_startup_uses_current_softstart_uvlo_switch_cycle_timestep_energy_waveform_and_result": (
            _smps_startup_owner_identity_ok(positive)
        ),
        "noise_bands_use_current_sources_correlation_psd_grid_bandwidth_transfer_integration_model_and_result": (
            _noise_correlation_band_owner_identity_ok(positive)
        ),
        "transmission_lines_use_current_z0_delay_reflections_arrivals_polarity_causality_energy_waveform_and_result": (
            _transmission_line_transient_owner_identity_ok(positive)
        ),
        "sampled_loops_use_current_sideband_injection_period_fft_phase_nyquist_crossover_waveform_and_result": (
            _sampled_loop_owner_identity_ok(positive)
        ),
        "periodic_steady_state_uses_current_charge_flux_energy_efficiency_phase_waveform_owner_and_result": (
            _periodic_steady_state_owner_identity_ok(positive)
        ),
        "small_signal_uses_current_bias_jacobian_ac_transient_poles_zeros_normalization_circuit_and_result": (
            _small_signal_consistency_owner_identity_ok(positive)
        ),
        "switched_limit_cycles_use_current_poincare_floquet_events_flux_charge_energy_step_waveform_and_result": (
            _switched_limit_cycle_owner_identity_ok(positive)
        ),
        "noise_referral_uses_current_input_output_density_correlation_band_gain_circuit_and_result": (
            _noise_referral_owner_identity_ok(positive)
        ),
        "pll_uses_current_lock_phase_error_noise_jitter_loop_gain_waveform_circuit_and_result": (
            _pll_closure_owner_identity_ok(positive)
        ),
        "electrothermal_uses_current_thermal_impedance_power_temperature_device_fixedpoint_stability_circuit_and_result": (
            _electrothermal_fixedpoint_owner_identity_ok(positive)
        ),
        "oscillators_use_current_startup_barkhausen_frequency_amplitude_limitcycle_energy_waveform_and_result": (
            _oscillator_closure_owner_identity_ok(positive)
        ),
        "conducted_emi_uses_current_lisn_window_fft_detector_band_limit_power_circuit_and_result": (
            _conducted_emi_owner_identity_ok(positive)
        ),
        "sallen_key_uses_current_poles_q_gain_opamp_gbw_phase_noise_step_circuit_owner_and_result": (
            _sallen_key_owner_identity_ok(positive)
        ),
        "bridge_rectifier_uses_current_inrush_ripple_charge_conduction_loss_load_cycle_energy_waveform_and_result": (
            _bridge_rectifier_owner_identity_ok(positive)
        ),
        "flyback_uses_current_magnetizing_leakage_snubber_demag_flux_current_cycle_energy_waveform_and_result": (
            _flyback_owner_identity_ok(positive)
        ),
        "transimpedance_uses_current_capacitance_noise_gain_bandwidth_stability_noise_step_circuit_and_result": (
            _transimpedance_owner_identity_ok(positive)
        ),
        "llc_resonant_converters_use_current_elements_frequency_gain_zvs_current_loss_cycle_energy_waveform_and_result": (
            _llc_resonant_owner_identity_ok(positive)
        ),
        "bjt_amplifiers_use_current_bias_gm_gain_pole_noise_distortion_thermal_power_circuit_and_result": (
            _bjt_amplifier_owner_identity_ok(positive)
        ),
        "boost_converters_use_current_duty_ripple_ccm_output_efficiency_stress_energy_waveform_and_result": (
            _boost_converter_owner_identity_ok(positive)
        ),
        "active_filters_use_current_poles_zeros_q_gain_noise_slew_saturation_power_circuit_and_result": (
            _active_filter_owner_identity_ok(positive)
        ),
        "mosfet_gate_drives_use_current_charge_current_deadtime_losses_temperature_cycle_energy_waveform_and_result": (
            _mosfet_gatedrive_owner_identity_ok(positive)
        ),
        "instrumentation_amplifiers_use_current_gain_cmrr_input_range_noise_headroom_power_circuit_and_result": (
            _instrumentation_amplifier_owner_identity_ok(positive)
        ),
        "flyback_converters_use_current_turnsratio_magnetizing_leakage_clamp_stress_loss_power_energy_waveform_and_result": (
            _flyback_v43_identity_ok(positive)
        ),
        "opamps_use_current_loopgain_phasemargin_crossover_step_overshoot_slew_power_waveform_and_result": (
            _opamp_stability_v43_identity_ok(positive)
        ),
        "v45_public_replays_bind_release_owner_and_digest": validate_ltspice_v45_identity(
            positive
        ),
        "v46_public_replays_bind_partial_state_seed_and_finite_artifacts": validate_ltspice_v46_identity(
            positive
        ),
        "v47_public_replays_bind_hierarchy_current_kcl_and_analysis_rows": validate_ltspice_v47_identity(
            positive
        ),
        "v48_public_replays_bind_behavioral_events_and_fourier_basis": validate_ltspice_v48_identity(
            positive
        ),
        "v49_public_replays_bind_monte_carlo_and_switch_state": validate_ltspice_v49_identity(
            positive
        ),
        "v50_public_replays_bind_step_rows_and_noise_traces": validate_ltspice_v50_identity(
            positive
        ),
        "exactly_four_timing_stages": timing_ok,
    }
    checks.update(validate_ltspice_v44_public_identity(positive))
    return {
        "schema": "radia-spice-lab.ideal-transformer-identity.v1",
        "policy": "ideal_transformer_identity_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            **derived,
            "maximum_declared_expectation_relative_error": declared_expectation_error,
            "maximum_recomputed_identity_relative_error": max(recomputed.values()),
            "maximum_reported_error_absolute_drift": reported_error_drift,
            "maximum_reported_power_relative_drift": power_report_drift,
            "transformer_input_active_power_w": transformer_input_power.real,
            "transformer_output_active_power_w": transformer_output_power.real,
            "source_delivery_active_power_w": source_delivery_power.real,
        },
        "notes": [
            "For the chosen dot/sign convention, V_secondary=-V_primary/N and I_secondary=-N*I_primary.",
            "The load reflected to the primary is N^2*R_load; include source resistance in the analytic RMS reference.",
            "Turns-ratio agreement alone is insufficient: require instantaneous and complex-power conservation plus deterministic replay.",
            "When phasor frequency or terminal-current role metadata is supplied, bind it to the source and sign contract before accepting scalar power closure.",
            "A matching scalar phasor is not replay evidence when its fitted traces mix RAW generations or cross a transient restart.",
            "Peak, RMS, and phasor quantities must declare one normalization basis before turns-ratio or complex-power closure is compared.",
            "Phase replay must retain the fitted unwrap branch orientation; equal magnitudes do not resolve a sign-changing branch alias.",
        ],
    }
