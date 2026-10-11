"""Internal signal evidence checks for LTspice identities."""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from ._identity_common import (
    _finite,
    _is_sha256,
)


def _phasor_replay_generations_ok(positive: Mapping[str, object]) -> bool:
    evidence = positive.get("phasor_replay_evidence")
    if evidence is None:
        return True
    if (
        not isinstance(evidence, Sequence)
        or isinstance(evidence, (str, bytes))
        or len(evidence) < 2
    ):
        return False
    rows = [row for row in evidence if isinstance(row, Mapping)]
    return len(rows) == len(evidence) and all(
        bool(row.get("fit_id"))
        and bool(row.get("raw_generation_id"))
        and row.get("trace_group_generation_id") == row.get("raw_generation_id")
        and len(str(row.get("scalar_phasor_digest") or "")) == 64
        for row in rows
    )


def _fit_window_stays_in_one_segment(
    positive: Mapping[str, object], fit_start: float, fit_stop: float
) -> bool:
    contract = positive.get("fit_window_segment_contract")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("fit_run_generation_id") or "")
    sample_generations = contract.get("sample_run_generation_ids")
    restarts = contract.get("restart_discontinuities_s")
    if (
        not generation
        or not isinstance(sample_generations, Sequence)
        or isinstance(sample_generations, (str, bytes))
        or not sample_generations
        or not isinstance(restarts, Sequence)
        or isinstance(restarts, (str, bytes))
        or any(str(item) != generation for item in sample_generations)
    ):
        return False
    try:
        restart_times = [float(value) for value in restarts]
    except (TypeError, ValueError):
        return False
    return all(
        math.isfinite(value) and not (fit_start < value < fit_stop)
        for value in restart_times
    )


def _phasor_basis_contract_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get("phasor_basis_contract")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    bases = contract.get("quantity_basis")
    factors = contract.get("normalization_factor_to_rms")
    names = {
        "source_voltage_phasor_rms_v",
        "primary_voltage_phasor_rms_v",
        "secondary_voltage_phasor_rms_v",
        "source_delivery_current_phasor_rms_a",
        "primary_current_phasor_rms_a",
        "secondary_current_phasor_rms_a",
    }
    if not isinstance(bases, Mapping) or not isinstance(factors, Mapping):
        return False
    try:
        factors_are_rms = all(
            math.isclose(float(factors.get(name)), 1.0, rel_tol=0.0, abs_tol=1.0e-15)
            for name in names
        )
    except (TypeError, ValueError):
        return False
    return (
        set(bases) == names
        and set(factors) == names
        and all(bases.get(name) == "rms" for name in names)
        and factors_are_rms
        and contract.get("complex_power_formula")
        == "rms_voltage_times_conjugate_rms_current"
    )


def _phase_unwrap_contract_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get("phase_unwrap_contract")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        period = float(contract.get("branch_period_rad"))
    except (TypeError, ValueError):
        return False
    fitted = str(contract.get("fit_branch_sign_digest") or "")
    replayed = str(contract.get("replay_branch_sign_digest") or "")
    return (
        contract.get("phase_unit") == "radian"
        and math.isclose(period, 2.0 * math.pi, rel_tol=0.0, abs_tol=1.0e-12)
        and contract.get("unwrap_convention") == "continuous_signed_phase"
        and contract.get("reference_trace") == "source_voltage_phasor_rms_v"
        and len(fitted) == 64
        and replayed == fitted
    )


def _power_sign_convention_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get("power_sign_convention")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    return (
        contract.get("source_power_role") == "delivered_positive"
        and contract.get("passive_power_role") == "absorbed_positive"
        and contract.get("source_current_reference")
        == "leaving_positive_terminal"
        and contract.get("passive_current_reference")
        == "entering_positive_terminal"
        and contract.get("balance_equation")
        == "source_delivered_equals_passive_absorbed"
        and contract.get("recorded_sign_transform") == "none"
    )


def _ac_frequency_interpolation_contract_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("ac_frequency_interpolation_contract")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    source_coordinate = str(contract.get("source_frequency_coordinate") or "")
    target_coordinate = str(contract.get("target_frequency_coordinate") or "")
    interpolation_coordinate = str(contract.get("interpolation_coordinate") or "")
    source_generation = str(contract.get("source_grid_generation") or "")
    trace_generation = str(
        contract.get("interpolated_trace_source_grid_generation") or ""
    )
    return (
        source_coordinate in {"linear_hz", "log10_hz", "natural_log_hz"}
        and source_coordinate == target_coordinate == interpolation_coordinate
        and contract.get("frequency_unit") == "Hz"
        and bool(source_generation)
        and trace_generation == source_generation
    )


def _transient_energy_window_event_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("transient_energy_window_event_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        energy_index = int(contract.get("energy_window_start_event_index"))
        reference_index = int(contract.get("reference_window_start_event_index"))
        energy_duration = float(contract.get("energy_window_duration_s"))
        reference_duration = float(contract.get("reference_window_duration_s"))
    except (TypeError, ValueError):
        return False
    event_generation = str(contract.get("event_detection_generation") or "")
    return (
        bool(str(contract.get("event_type") or ""))
        and energy_index >= 0
        and reference_index == energy_index
        and math.isfinite(energy_duration)
        and energy_duration > 0.0
        and math.isclose(
            reference_duration, energy_duration, rel_tol=0.0, abs_tol=0.0
        )
        and bool(event_generation)
        and contract.get("energy_window_event_generation") == event_generation
        and contract.get("reference_window_event_generation") == event_generation
    )


def _noise_density_band_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get("noise_density_band_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    scales = {"Hz": 1.0, "kHz": 1.0e3, "MHz": 1.0e6, "GHz": 1.0e9}
    numerical_unit = str(contract.get("numerical_frequency_unit") or "")
    band_unit = str(contract.get("band_limit_unit") or "")
    try:
        numerical_scale = float(contract.get("numerical_frequency_scale_to_hz"))
        band_scale = float(contract.get("band_limit_scale_to_hz"))
    except (TypeError, ValueError):
        return False
    grid_generation = str(contract.get("frequency_grid_generation") or "")
    return (
        contract.get("noise_density_unit") == "V/sqrt(Hz)"
        and numerical_unit in scales
        and band_unit in scales
        and math.isclose(
            numerical_scale, scales[numerical_unit], rel_tol=0.0, abs_tol=0.0
        )
        and math.isclose(band_scale, scales[band_unit], rel_tol=0.0, abs_tol=0.0)
        and contract.get("integrated_noise_unit") == "V_rms"
        and bool(grid_generation)
        and contract.get("band_integration_grid_generation") == grid_generation
    )


def _steady_cycle_average_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get("steady_cycle_average_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("settled_cycle_generation_id") or "")
    try:
        first_settled = int(contract.get("first_settled_cycle_index"))
        average_cycle = int(contract.get("average_window_cycle_index"))
        period_count = int(contract.get("average_window_period_count"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and contract.get("waveform_cycle_generation_id") == generation
        and contract.get("average_window_cycle_generation_id") == generation
        and first_settled >= 0
        and average_cycle >= first_settled
        and period_count >= 1
    )


def _monte_carlo_parameter_seed_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("monte_carlo_parameter_seed_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    parameter_order = contract.get("parameter_order")
    statistics_order = contract.get("statistics_parameter_order")
    if (
        not isinstance(parameter_order, Sequence)
        or isinstance(parameter_order, (str, bytes))
        or not isinstance(statistics_order, Sequence)
        or isinstance(statistics_order, (str, bytes))
    ):
        return False
    names = [str(name) for name in parameter_order]
    replayed_names = [str(name) for name in statistics_order]
    generation = str(contract.get("seed_schedule_generation_id") or "")
    seed_map_digest = str(contract.get("parameter_seed_map_sha256") or "")
    return (
        bool(names)
        and len(set(names)) == len(names)
        and replayed_names == names
        and contract.get("seed_policy") == "one_seed_per_named_parameter"
        and bool(generation)
        and contract.get("statistics_seed_schedule_generation_id") == generation
        and len(seed_map_digest) == 64
        and all(character in "0123456789abcdef" for character in seed_map_digest)
        and contract.get("statistics_parameter_seed_map_sha256")
        == seed_map_digest
    )


def _ac_phase_coordinate_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get("ac_phase_coordinate_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("ac_sweep_generation_id") or "")
    unit = str(contract.get("phase_coordinate_unit") or "")
    convention = str(contract.get("phase_coordinate_convention") or "")
    value_digest = str(contract.get("phase_trace_value_sha256") or "")
    expected_conventions = {
        "degree": "principal_degree_minus180_180",
        "radian": "principal_radian_minuspi_pi",
    }
    return (
        bool(generation)
        and contract.get("phase_trace_sweep_generation_id") == generation
        and contract.get("reference_phase_trace_sweep_generation_id")
        == generation
        and unit in expected_conventions
        and contract.get("reference_phase_coordinate_unit") == unit
        and convention == expected_conventions[unit]
        and contract.get("reference_phase_coordinate_convention") == convention
        and _is_sha256(value_digest)
        and contract.get("reference_phase_trace_value_sha256") == value_digest
    )


def _transient_derivative_adaptive_history_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("transient_derivative_adaptive_history_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("transient_generation_id") or "")
    time_grid_digest = str(contract.get("accepted_time_grid_sha256") or "")
    scheme = str(contract.get("derivative_scheme") or "")
    try:
        accepted_count = int(contract.get("accepted_step_count"))
        history_count = int(contract.get("derivative_history_step_count"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and contract.get("accepted_step_generation_id") == generation
        and contract.get("derivative_history_generation_id") == generation
        and contract.get("current_sample_generation_id") == generation
        and contract.get("derivative_sample_generation_id") == generation
        and _is_sha256(time_grid_digest)
        and contract.get("derivative_history_time_grid_sha256")
        == time_grid_digest
        and accepted_count > 2
        and history_count == accepted_count
        and scheme == "variable_step_bdf2"
        and contract.get("history_derivative_scheme") == scheme
    )


def _noise_spectral_density_sidedness_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("noise_spectral_density_sidedness_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("noise_generation_id") or "")
    target_basis = str(contract.get("density_sidedness_basis") or "")
    reference_basis = str(
        contract.get("reference_density_sidedness_basis") or ""
    )
    bases = {
        "one_sided_positive_frequency",
        "two_sided_full_frequency",
    }
    if target_basis not in bases or reference_basis not in bases:
        return False
    if target_basis == reference_basis:
        expected_scale = 1.0
    elif reference_basis == "two_sided_full_frequency":
        expected_scale = math.sqrt(2.0)
    else:
        expected_scale = 1.0 / math.sqrt(2.0)
    try:
        scale = _finite(
            contract.get("reference_to_density_amplitude_scale"),
            "reference_to_density_amplitude_scale",
            positive=True,
        )
    except ValueError:
        return False
    value_digest = str(contract.get("density_trace_value_sha256") or "")
    return (
        bool(generation)
        and contract.get("density_trace_generation_id") == generation
        and contract.get("reference_density_trace_generation_id") == generation
        and contract.get("density_quantity") == "amplitude_spectral_density"
        and contract.get("density_unit") == "V/sqrt(Hz)"
        and contract.get("integration_sidedness_basis") == target_basis
        and math.isclose(scale, expected_scale, rel_tol=1.0e-12, abs_tol=0.0)
        and _is_sha256(value_digest)
        and contract.get("reference_density_trace_value_sha256") == value_digest
    )


def _stepped_parameter_interpolation_coordinate_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("stepped_parameter_interpolation_coordinate_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    source_values = contract.get("source_parameter_values")
    interpolation_values = contract.get("interpolation_parameter_values")
    if (
        not isinstance(source_values, Sequence)
        or isinstance(source_values, (str, bytes))
        or not isinstance(interpolation_values, Sequence)
        or isinstance(interpolation_values, (str, bytes))
    ):
        return False
    try:
        source = [float(value) for value in source_values]
        interpolated = [float(value) for value in interpolation_values]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("step_generation_id") or "")
    coordinate = str(contract.get("source_parameter_coordinate") or "")
    axis_digest = str(contract.get("source_parameter_axis_sha256") or "")
    return (
        bool(generation)
        and contract.get("source_step_generation_id") == generation
        and contract.get("interpolated_result_step_generation_id") == generation
        and bool(str(contract.get("parameter_name") or ""))
        and bool(str(contract.get("parameter_unit") or ""))
        and coordinate in {"linear_value", "log10_value"}
        and contract.get("target_parameter_coordinate") == coordinate
        and contract.get("interpolation_parameter_coordinate") == coordinate
        and len(source) >= 2
        and len(interpolated) == len(source)
        and all(math.isfinite(value) and value > 0.0 for value in source)
        and all(right > left for left, right in zip(source, source[1:]))
        and interpolated == source
        and _is_sha256(axis_digest)
        and contract.get("interpolation_source_axis_sha256") == axis_digest
    )


def _fft_window_coherent_gain_amplitude_basis_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("fft_window_coherent_gain_amplitude_basis_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    trace_generation = str(contract.get("trace_generation_id") or "")
    window_generation = str(contract.get("window_generation_id") or "")
    input_basis = str(contract.get("input_amplitude_basis") or "")
    result_basis = str(contract.get("fft_result_amplitude_basis") or "")
    digest = str(contract.get("window_coefficients_sha256") or "")
    try:
        sample_count = int(contract.get("sample_count"))
        coherent_gain = _finite(
            contract.get("coherent_gain"), "coherent_gain", positive=True
        )
        correction = _finite(
            contract.get("coherent_gain_correction"),
            "coherent_gain_correction",
            positive=True,
        )
        conversion_count = int(contract.get("amplitude_basis_conversion_count"))
    except (TypeError, ValueError):
        return False
    bases = {"peak", "rms"}
    expected_conversions = 0 if input_basis == result_basis else 1
    return (
        bool(trace_generation)
        and contract.get("fft_input_trace_generation_id") == trace_generation
        and bool(window_generation)
        and contract.get("coherent_gain_window_generation_id")
        == window_generation
        and contract.get("fft_result_window_generation_id") == window_generation
        and contract.get("window_definition") == "periodic_hann"
        and sample_count >= 4
        and math.isclose(coherent_gain, 0.5, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(
            correction, 1.0 / coherent_gain, rel_tol=1.0e-15, abs_tol=0.0
        )
        and input_basis in bases
        and result_basis in bases
        and conversion_count == expected_conversions
        and _is_sha256(digest)
        and contract.get("fft_window_coefficients_sha256") == digest
    )


def _monte_carlo_percentile_sample_filter_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("monte_carlo_percentile_sample_filter_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    excluded = contract.get("excluded_sample_ids")
    if (
        not isinstance(excluded, Sequence)
        or isinstance(excluded, (str, bytes))
    ):
        return False
    try:
        raw_count = int(contract.get("raw_sample_count"))
        included_count = int(contract.get("included_sample_count"))
        percentile = _finite(contract.get("percentile"), "percentile")
        excluded_ids = [int(sample_id) for sample_id in excluded]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("statistics_generation_id") or "")
    digest = str(contract.get("sample_filter_sha256") or "")
    return (
        bool(generation)
        and contract.get("raw_sample_statistics_generation_id") == generation
        and contract.get("sample_filter_statistics_generation_id") == generation
        and contract.get("percentile_statistics_generation_id") == generation
        and raw_count > 0
        and 0 < included_count <= raw_count
        and len(excluded_ids) == raw_count - included_count
        and len(set(excluded_ids)) == len(excluded_ids)
        and all(1 <= sample_id <= raw_count for sample_id in excluded_ids)
        and 0.0 <= percentile <= 100.0
        and contract.get("sample_filter_policy") == "finite_converged_only"
        and contract.get("percentile_sample_filter_policy")
        == "finite_converged_only"
        and _is_sha256(digest)
        and contract.get("percentile_sample_filter_sha256") == digest
    )


def _measure_crossing_interpolation_grid_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "measure_crossing_interpolation_time_grid_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    bracket = contract.get("bracket_sample_indices")
    interpolation_bracket = contract.get(
        "interpolation_bracket_sample_indices"
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in (bracket, interpolation_bracket)
    ):
        return False
    try:
        ordinal = int(contract.get("crossing_ordinal"))
        crossing_time = _finite(contract.get("crossing_time_s"), "crossing_time_s")
        reported_time = _finite(
            contract.get("reported_crossing_time_s"),
            "reported_crossing_time_s",
        )
        bracket_indices = [int(index) for index in bracket]
        interpolation_indices = [int(index) for index in interpolation_bracket]
    except (TypeError, ValueError):
        return False
    transient_generation = str(contract.get("transient_generation_id") or "")
    grid_generation = str(
        contract.get("accepted_step_grid_generation_id") or ""
    )
    grid_digest = str(contract.get("accepted_step_grid_sha256") or "")
    return (
        bool(transient_generation)
        and contract.get("measure_generation_id") == transient_generation
        and bool(grid_generation)
        and contract.get("interpolation_grid_generation_id") == grid_generation
        and contract.get("crossing_bracket_grid_generation_id") == grid_generation
        and ordinal > 0
        and contract.get("crossing_direction") in {"rising", "falling", "either"}
        and contract.get("interpolation_method") == "linear"
        and crossing_time >= 0.0
        and math.isclose(reported_time, crossing_time, rel_tol=1.0e-12, abs_tol=0.0)
        and len(bracket_indices) == 2
        and bracket_indices[0] >= 0
        and bracket_indices[1] == bracket_indices[0] + 1
        and interpolation_indices == bracket_indices
        and _is_sha256(grid_digest)
        and contract.get("interpolation_grid_sha256") == grid_digest
    )


def _fourier_phase_reference_time_origin_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "fourier_harmonic_phase_reference_time_origin_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        frequency = _finite(
            contract.get("fundamental_frequency_hz"),
            "fundamental_frequency_hz",
            positive=True,
        )
        harmonic = int(contract.get("harmonic_number"))
        origin = _finite(
            contract.get("reference_time_origin_s"), "reference_time_origin_s"
        )
        fourier_origin = _finite(
            contract.get("fourier_reference_time_origin_s"),
            "fourier_reference_time_origin_s",
        )
        comparison_origin = _finite(
            contract.get("comparison_reference_time_origin_s"),
            "comparison_reference_time_origin_s",
        )
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("waveform_generation_id") or "")
    phase_basis = str(contract.get("phase_basis") or "")
    time_digest = str(contract.get("waveform_time_axis_sha256") or "")
    return (
        bool(generation)
        and contract.get("fourier_result_waveform_generation_id") == generation
        and contract.get("harmonic_table_waveform_generation_id") == generation
        and frequency > 0.0
        and harmonic > 0
        and phase_basis in {"cosine", "sine"}
        and contract.get("reported_phase_basis") == phase_basis
        and math.isclose(fourier_origin, origin, rel_tol=0.0, abs_tol=1.0e-18)
        and math.isclose(comparison_origin, origin, rel_tol=0.0, abs_tol=1.0e-18)
        and contract.get("reference_time_origin_convention")
        == "absolute_transient_time"
        and contract.get("comparison_time_origin_convention")
        == "absolute_transient_time"
        and _is_sha256(time_digest)
        and contract.get("fourier_time_axis_sha256") == time_digest
    )


def _ac_group_delay_phase_unwrap_grid_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "ac_group_delay_phase_unwrap_frequency_grid_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    sequences = (
        contract.get("frequency_hz"),
        contract.get("phase_unwrapped_rad"),
        contract.get("group_delay_s"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in sequences
    ):
        return False
    try:
        frequencies = [_finite(value, "frequency_hz", positive=True) for value in sequences[0]]
        phases = [_finite(value, "phase_unwrapped_rad") for value in sequences[1]]
        delays = [_finite(value, "group_delay_s") for value in sequences[2]]
        anchor = _finite(
            contract.get("phase_unwrap_branch_anchor_rad"),
            "phase_unwrap_branch_anchor_rad",
        )
        delay_anchor = _finite(
            contract.get("group_delay_phase_unwrap_branch_anchor_rad"),
            "group_delay_phase_unwrap_branch_anchor_rad",
        )
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("frequency_grid_generation_id") or "")
    digest = str(contract.get("frequency_grid_sha256") or "")
    computed_delays = [
        -(right_phase - left_phase) / (2.0 * math.pi * (right_f - left_f))
        for left_f, right_f, left_phase, right_phase in zip(
            frequencies, frequencies[1:], phases, phases[1:]
        )
    ]
    return (
        bool(str(contract.get("ac_sweep_generation_id") or ""))
        and bool(generation)
        and contract.get("phase_sample_frequency_grid_generation_id") == generation
        and contract.get("phase_unwrap_frequency_grid_generation_id") == generation
        and contract.get("group_delay_frequency_grid_generation_id") == generation
        and len(frequencies) >= 2
        and len(phases) == len(frequencies)
        and len(delays) == len(frequencies) - 1
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and all(
            math.isclose(reported, computed, rel_tol=1.0e-8, abs_tol=1.0e-15)
            for reported, computed in zip(delays, computed_delays)
        )
        and math.isclose(delay_anchor, anchor, rel_tol=0.0, abs_tol=1.0e-15)
        and contract.get("phase_unwrap_method") == "continuous_minimum_jump"
        and contract.get("group_delay_phase_unwrap_method")
        == "continuous_minimum_jump"
        and _is_sha256(digest)
        and contract.get("phase_unwrap_frequency_grid_sha256") == digest
    )


def _transient_rms_average_event_window_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("transient_rms_average_event_window_generation_identity")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        start = _finite(contract.get("window_start_s"), "window_start_s")
        rms_start = _finite(contract.get("rms_window_start_s"), "rms_window_start_s")
        average_start = _finite(
            contract.get("average_window_start_s"), "average_window_start_s"
        )
        end = _finite(contract.get("window_end_s"), "window_end_s", positive=True)
        rms_end = _finite(contract.get("rms_window_end_s"), "rms_window_end_s")
        average_end = _finite(
            contract.get("average_window_end_s"), "average_window_end_s"
        )
    except (TypeError, ValueError):
        return False
    event_generation = str(contract.get("switching_event_generation_id") or "")
    start_event = str(contract.get("window_start_event_id") or "")
    end_event = str(contract.get("window_end_event_id") or "")
    digest = str(contract.get("event_table_sha256") or "")
    return (
        bool(str(contract.get("transient_generation_id") or ""))
        and bool(event_generation)
        and contract.get("rms_window_event_generation_id") == event_generation
        and contract.get("average_window_event_generation_id") == event_generation
        and bool(start_event)
        and contract.get("rms_window_start_event_id") == start_event
        and contract.get("average_window_start_event_id") == start_event
        and bool(end_event)
        and end_event != start_event
        and contract.get("rms_window_end_event_id") == end_event
        and contract.get("average_window_end_event_id") == end_event
        and start < end
        and math.isclose(rms_start, start, rel_tol=0.0, abs_tol=1.0e-18)
        and math.isclose(average_start, start, rel_tol=0.0, abs_tol=1.0e-18)
        and math.isclose(rms_end, end, rel_tol=0.0, abs_tol=1.0e-18)
        and math.isclose(average_end, end, rel_tol=0.0, abs_tol=1.0e-18)
        and _is_sha256(digest)
        and contract.get("rms_event_table_sha256") == digest
        and contract.get("average_event_table_sha256") == digest
    )


def _ac_noise_integrated_density_bin_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "ac_noise_integrated_density_sidedness_bin_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    rows = (
        contract.get("frequency_hz"),
        contract.get("frequency_bin_width_hz"),
        contract.get("integration_frequency_bin_width_hz"),
        contract.get("noise_density_v_per_sqrt_hz"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in rows
    ):
        return False
    try:
        frequencies = [_finite(value, "frequency_hz", positive=True) for value in rows[0]]
        widths = [
            _finite(value, "frequency_bin_width_hz", positive=True)
            for value in rows[1]
        ]
        integration_widths = [
            _finite(value, "integration_frequency_bin_width_hz", positive=True)
            for value in rows[2]
        ]
        densities = [
            _finite(value, "noise_density_v_per_sqrt_hz", positive=True)
            for value in rows[3]
        ]
        amplitude_factor = _finite(
            contract.get("density_to_integration_amplitude_factor"),
            "density_to_integration_amplitude_factor",
            positive=True,
        )
        integrated = _finite(
            contract.get("integrated_noise_rms_v"),
            "integrated_noise_rms_v",
            positive=True,
        )
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("frequency_bin_generation_id") or "")
    digest = str(contract.get("frequency_bin_table_sha256") or "")
    recomputed = math.sqrt(
        sum(
            (amplitude_factor * density) ** 2 * width
            for density, width in zip(densities, widths)
        )
    )
    return (
        bool(str(contract.get("noise_generation_id") or ""))
        and bool(generation)
        and contract.get("density_frequency_bin_generation_id") == generation
        and contract.get("integration_frequency_bin_generation_id") == generation
        and contract.get("sidedness_conversion_frequency_bin_generation_id")
        == generation
        and len(frequencies) == len(widths) == len(integration_widths) == len(densities)
        and len(frequencies) >= 2
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and integration_widths == widths
        and contract.get("density_sidedness_basis")
        == "one_sided_positive_frequency"
        and contract.get("integration_sidedness_basis")
        == "one_sided_positive_frequency"
        and math.isclose(amplitude_factor, 1.0, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(integrated, recomputed, rel_tol=1.0e-12, abs_tol=1.0e-24)
        and _is_sha256(digest)
        and contract.get("integration_frequency_bin_table_sha256") == digest
    )


def _transient_power_interpolation_grid_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "transient_power_voltage_current_interpolation_grid_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    rows = (
        contract.get("voltage_sample_time_s"),
        contract.get("current_sample_time_s"),
        contract.get("power_interpolation_time_s"),
        contract.get("interpolated_voltage_v"),
        contract.get("interpolated_current_a"),
        contract.get("instantaneous_power_w"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in rows
    ):
        return False
    try:
        voltage_time = [_finite(value, "voltage_sample_time_s") for value in rows[0]]
        current_time = [_finite(value, "current_sample_time_s") for value in rows[1]]
        target_time = [_finite(value, "power_interpolation_time_s") for value in rows[2]]
        voltage = [_finite(value, "interpolated_voltage_v") for value in rows[3]]
        current = [_finite(value, "interpolated_current_a") for value in rows[4]]
        power = [_finite(value, "instantaneous_power_w") for value in rows[5]]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("voltage_sample_grid_generation_id") or "")
    digest = str(contract.get("voltage_current_grid_sha256") or "")
    expected_power = [volts * amps for volts, amps in zip(voltage, current)]
    return (
        bool(str(contract.get("transient_generation_id") or ""))
        and bool(generation)
        and contract.get("current_sample_grid_generation_id") == generation
        and contract.get("power_interpolation_grid_generation_id") == generation
        and contract.get("integration_grid_generation_id") == generation
        and len(voltage_time) == len(current_time) == len(target_time)
        == len(voltage) == len(current) == len(power)
        and len(target_time) >= 2
        and all(right > left for left, right in zip(target_time, target_time[1:]))
        and voltage_time == target_time
        and current_time == target_time
        and all(
            math.isclose(reported, expected, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for reported, expected in zip(power, expected_power)
        )
        and contract.get("power_sign_convention") == "passive_absorbed_positive"
        and _is_sha256(digest)
        and contract.get("power_interpolation_grid_sha256") == digest
    )


def _stepped_ac_parameter_tuple_grid_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "stepped_ac_trace_parameter_tuple_interpolation_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    names = contract.get("parameter_names")
    tuples = contract.get("parameter_tuples")
    trace_tuples = contract.get("trace_parameter_tuples")
    frequency = contract.get("frequency_hz")
    trace_frequency = contract.get("trace_frequency_hz")
    interpolation_frequency = contract.get("interpolation_frequency_hz")
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in (
            names,
            tuples,
            trace_tuples,
            frequency,
            trace_frequency,
            interpolation_frequency,
        )
    ):
        return False
    name_rows = [str(value) for value in names]
    try:
        tuple_rows = [
            [_finite(value, "parameter_tuple") for value in row]
            for row in tuples
            if isinstance(row, Sequence) and not isinstance(row, (str, bytes))
        ]
        trace_tuple_rows = [
            [_finite(value, "trace_parameter_tuple") for value in row]
            for row in trace_tuples
            if isinstance(row, Sequence) and not isinstance(row, (str, bytes))
        ]
        frequency_rows = [
            _finite(value, "frequency_hz", positive=True) for value in frequency
        ]
        trace_frequency_rows = [
            _finite(value, "trace_frequency_hz", positive=True)
            for value in trace_frequency
        ]
        interpolation_frequency_rows = [
            _finite(value, "interpolation_frequency_hz", positive=True)
            for value in interpolation_frequency
        ]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("sweep_generation_id") or "")
    tuple_digest = str(contract.get("parameter_tuple_table_sha256") or "")
    grid_digest = str(contract.get("trace_frequency_grid_sha256") or "")
    return (
        bool(generation)
        and contract.get("parameter_tuple_sweep_generation_id") == generation
        and contract.get("trace_grid_sweep_generation_id") == generation
        and contract.get("interpolator_sweep_generation_id") == generation
        and bool(name_rows)
        and all(name_rows)
        and len(set(name_rows)) == len(name_rows)
        and len(tuple_rows) == len(tuples) == len(trace_tuple_rows)
        and bool(tuple_rows)
        and all(len(row) == len(name_rows) for row in tuple_rows)
        and trace_tuple_rows == tuple_rows
        and len(frequency_rows) >= 2
        and all(
            right > left for left, right in zip(frequency_rows, frequency_rows[1:])
        )
        and trace_frequency_rows == frequency_rows
        and interpolation_frequency_rows == frequency_rows
        and _is_sha256(tuple_digest)
        and contract.get("trace_parameter_tuple_table_sha256") == tuple_digest
        and _is_sha256(grid_digest)
        and contract.get("interpolation_frequency_grid_sha256") == grid_digest
    )


def _measure_trigger_target_crossing_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "measure_trigger_target_crossing_edge_count_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    trigger_bracket = contract.get("trigger_bracket_time_s")
    target_bracket = contract.get("target_bracket_time_s")
    if not all(
        isinstance(values, Sequence)
        and not isinstance(values, (str, bytes))
        and len(values) == 2
        for values in (trigger_bracket, target_bracket)
    ):
        return False
    try:
        trigger_bounds = [
            _finite(value, "trigger_bracket_time_s") for value in trigger_bracket
        ]
        target_bounds = [
            _finite(value, "target_bracket_time_s") for value in target_bracket
        ]
        trigger_time = _finite(
            contract.get("trigger_crossing_time_s"), "trigger_crossing_time_s"
        )
        target_time = _finite(
            contract.get("target_crossing_time_s"), "target_crossing_time_s"
        )
        trigger_count = int(
            _finite(
                contract.get("trigger_crossing_count"),
                "trigger_crossing_count",
                positive=True,
            )
        )
        target_count = int(
            _finite(
                contract.get("target_crossing_count"),
                "target_crossing_count",
                positive=True,
            )
        )
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("accepted_grid_generation_id") or "")
    digest = str(contract.get("accepted_grid_sha256") or "")
    return (
        bool(str(contract.get("transient_generation_id") or ""))
        and bool(generation)
        and contract.get("trigger_crossing_grid_generation_id") == generation
        and contract.get("target_crossing_grid_generation_id") == generation
        and contract.get("measure_interpolator_grid_generation_id") == generation
        and contract.get("trigger_edge") in {"rise", "fall", "cross"}
        and contract.get("target_edge") in {"rise", "fall", "cross"}
        and trigger_count == contract.get("trigger_crossing_count")
        and target_count == contract.get("target_crossing_count")
        and trigger_bounds[0] < trigger_bounds[1]
        and target_bounds[0] < target_bounds[1]
        and trigger_bounds[0] <= trigger_time <= trigger_bounds[1]
        and target_bounds[0] <= target_time <= target_bounds[1]
        and target_time > trigger_time
        and _is_sha256(digest)
        and contract.get("trigger_crossing_table_sha256") == digest
        and contract.get("target_crossing_table_sha256") == digest
    )


def _monte_carlo_sample_trace_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "monte_carlo_seed_sample_tuple_trace_row_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    rows = (
        contract.get("sample_ids"),
        contract.get("seed_order"),
        contract.get("trace_seed_order"),
        contract.get("parameter_names"),
        contract.get("sample_parameter_tuples"),
        contract.get("trace_sample_ids"),
        contract.get("trace_parameter_tuples"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in rows
    ):
        return False
    try:
        sample_ids = [int(value) for value in rows[0]]
        seeds = [int(value) for value in rows[1]]
        trace_seeds = [int(value) for value in rows[2]]
        names = [str(value) for value in rows[3]]
        samples = [[_finite(value, "sample") for value in row] for row in rows[4]]
        trace_ids = [int(value) for value in rows[5]]
        trace_samples = [
            [_finite(value, "trace_sample") for value in row] for row in rows[6]
        ]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("monte_carlo_generation_id") or "")
    digest = str(contract.get("sample_table_sha256") or "")
    return (
        bool(generation)
        and contract.get("seed_order_generation_id") == generation
        and contract.get("sample_tuple_generation_id") == generation
        and contract.get("trace_row_generation_id") == generation
        and bool(sample_ids)
        and all(value > 0 for value in sample_ids)
        and len(set(sample_ids)) == len(sample_ids)
        and len(seeds) == len(sample_ids)
        and all(value > 0 for value in seeds)
        and len(set(seeds)) == len(seeds)
        and trace_seeds == seeds
        and bool(names)
        and all(names)
        and len(set(names)) == len(names)
        and len(samples) == len(sample_ids)
        and all(len(row) == len(names) for row in samples)
        and trace_ids == sample_ids
        and trace_samples == samples
        and _is_sha256(digest)
        and contract.get("trace_sample_table_sha256") == digest
    )


def _fft_window_harmonic_bin_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "fft_window_sample_rate_harmonic_bin_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        sample_count = int(_finite(contract.get("sample_count"), "sample_count", positive=True))
        window_count = int(
            _finite(contract.get("window_sample_count"), "window_sample_count", positive=True)
        )
        sample_rate = _finite(contract.get("sample_rate_hz"), "sample_rate_hz", positive=True)
        fft_rate = _finite(contract.get("fft_sample_rate_hz"), "fft_sample_rate_hz", positive=True)
        gain = _finite(contract.get("coherent_gain"), "coherent_gain", positive=True)
        applied_gain = _finite(
            contract.get("applied_coherent_gain"), "applied_coherent_gain", positive=True
        )
        bins = [int(value) for value in contract["harmonic_bin_indices"]]
        fft_bins = [int(value) for value in contract["fft_harmonic_bin_indices"]]
        frequencies = [float(value) for value in contract["harmonic_frequencies_hz"]]
        fft_frequencies = [
            float(value) for value in contract["fft_harmonic_frequencies_hz"]
        ]
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(contract.get("transient_generation_id") or "")
    digest = str(contract.get("fft_contract_sha256") or "")
    expected_frequencies = [value * sample_rate / sample_count for value in bins]
    return (
        bool(generation)
        and contract.get("window_transient_generation_id") == generation
        and contract.get("sample_rate_transient_generation_id") == generation
        and contract.get("harmonic_bin_transient_generation_id") == generation
        and contract.get("fft_result_transient_generation_id") == generation
        and sample_count == contract.get("sample_count") == window_count
        and math.isclose(fft_rate, sample_rate, rel_tol=1.0e-12)
        and contract.get("window_name") in {"hann", "hamming", "blackman", "rectangular"}
        and 0.0 < gain <= 1.0
        and math.isclose(applied_gain, gain, rel_tol=1.0e-12)
        and bool(bins)
        and all(0 < value < sample_count // 2 for value in bins)
        and len(set(bins)) == len(bins)
        and fft_bins == bins
        and len(frequencies) == len(bins)
        and all(
            math.isclose(value, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for value, expected in zip(frequencies, expected_frequencies)
        )
        and fft_frequencies == frequencies
        and _is_sha256(digest)
        and contract.get("result_fft_contract_sha256") == digest
    )


def _noise_monte_carlo_psd_integration_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "noise_monte_carlo_sample_filter_psd_integration_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    rows = (
        contract.get("sample_ids"),
        contract.get("accepted_sample_mask"),
        contract.get("psd_sample_ids"),
        contract.get("integration_sample_ids"),
        contract.get("frequency_hz"),
        contract.get("bin_width_hz"),
        contract.get("integration_bin_width_hz"),
    )
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in rows
    ):
        return False
    try:
        sample_ids = [int(value) for value in rows[0]]
        accepted = list(rows[1])
        psd_ids = [int(value) for value in rows[2]]
        integration_ids = [int(value) for value in rows[3]]
        frequencies = [_finite(value, "frequency_hz", positive=True) for value in rows[4]]
        widths = [_finite(value, "bin_width_hz", positive=True) for value in rows[5]]
        integration_widths = [
            _finite(value, "integration_bin_width_hz", positive=True)
            for value in rows[6]
        ]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("monte_carlo_generation_id") or "")
    digest = str(contract.get("psd_table_sha256") or "")
    expected_integration_ids = [
        sample_id for sample_id, keep in zip(sample_ids, accepted) if keep
    ]
    return (
        bool(generation)
        and contract.get("sample_filter_monte_carlo_generation_id") == generation
        and contract.get("psd_table_monte_carlo_generation_id") == generation
        and contract.get("integration_monte_carlo_generation_id") == generation
        and bool(sample_ids)
        and all(value > 0 for value in sample_ids)
        and len(set(sample_ids)) == len(sample_ids)
        and len(accepted) == len(sample_ids)
        and all(isinstance(value, bool) for value in accepted)
        and psd_ids == sample_ids
        and integration_ids == expected_integration_ids
        and bool(integration_ids)
        and len(frequencies) == len(widths) == len(integration_widths)
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and integration_widths == widths
        and contract.get("psd_sidedness") == "one-sided"
        and contract.get("integration_psd_sidedness")
        == contract.get("psd_sidedness")
        and _is_sha256(digest)
        and contract.get("integration_input_sha256") == digest
    )


def _stepped_transient_measure_row_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "stepped_transient_accepted_grid_measure_row_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    keys = (
        "step_ids",
        "parameter_names",
        "parameter_tuples",
        "measure_step_ids",
        "measure_parameter_tuples",
        "measure_names",
        "accepted_time_grid_s",
        "measure_time_grid_s",
        "measure_row_keys",
        "decoded_measure_row_keys",
    )
    rows = tuple(contract.get(key) for key in keys)
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in rows
    ):
        return False
    try:
        step_ids = [int(value) for value in rows[0]]
        names = [str(value) for value in rows[1]]
        parameter_tuples = [
            [_finite(value, "parameter_tuple") for value in row] for row in rows[2]
        ]
        measure_step_ids = [int(value) for value in rows[3]]
        measure_tuples = [
            [_finite(value, "measure_parameter_tuple") for value in row]
            for row in rows[4]
        ]
        measure_names = [str(value) for value in rows[5]]
        accepted_grid = [_finite(value, "accepted_time_grid_s") for value in rows[6]]
        measure_grid = [_finite(value, "measure_time_grid_s") for value in rows[7]]
        row_keys = [str(value) for value in rows[8]]
        decoded_keys = [str(value) for value in rows[9]]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("transient_generation_id") or "")
    grid_digest = str(contract.get("accepted_grid_sha256") or "")
    tuple_digest = str(contract.get("parameter_tuple_table_sha256") or "")
    if len(step_ids) != len(parameter_tuples):
        return False
    step_parameters = dict(zip(step_ids, parameter_tuples))
    expected_measure_tuples = [step_parameters.get(step_id) for step_id in measure_step_ids]
    expected_keys = [
        f"{step_id}:{measure_name}"
        for step_id, measure_name in zip(measure_step_ids, measure_names)
    ]
    return (
        bool(generation)
        and contract.get("accepted_grid_transient_generation_id") == generation
        and contract.get("parameter_tuple_transient_generation_id") == generation
        and contract.get("measure_row_transient_generation_id") == generation
        and bool(step_ids)
        and all(value > 0 for value in step_ids)
        and len(set(step_ids)) == len(step_ids)
        and bool(names)
        and all(names)
        and len(set(names)) == len(names)
        and all(len(row) == len(names) for row in parameter_tuples)
        and bool(measure_step_ids)
        and len(measure_step_ids) == len(measure_tuples) == len(measure_names)
        and all(name for name in measure_names)
        and measure_tuples == expected_measure_tuples
        and len(row_keys) == len(measure_step_ids)
        and row_keys == expected_keys
        and decoded_keys == row_keys
        and len(accepted_grid) >= 2
        and accepted_grid[0] >= 0.0
        and all(right > left for left, right in zip(accepted_grid, accepted_grid[1:]))
        and measure_grid == accepted_grid
        and _is_sha256(grid_digest)
        and contract.get("measure_grid_sha256") == grid_digest
        and _is_sha256(tuple_digest)
        and contract.get("measure_parameter_tuple_table_sha256") == tuple_digest
    )


def _switched_converter_cycle_measure_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "switched_converter_cycle_measure_initial_state_topology_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    keys = (
        "state_variable_ids",
        "initial_state_values",
        "solver_initial_state_values",
        "accepted_time_grid_s",
        "measure_time_grid_s",
        "cycle_windows_s",
        "measure_cycle_windows_s",
        "measure_names",
        "reported_measure_names",
        "cycle_measure_values",
        "reported_cycle_measure_values",
    )
    rows = tuple(contract.get(key) for key in keys)
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in rows
    ):
        return False
    try:
        state_ids = [str(value) for value in rows[0]]
        initial_state = [_finite(value, "initial_state") for value in rows[1]]
        solver_state = [_finite(value, "solver_initial_state") for value in rows[2]]
        accepted_grid = [_finite(value, "accepted_time_grid_s") for value in rows[3]]
        measure_grid = [_finite(value, "measure_time_grid_s") for value in rows[4]]
        cycle_windows = [
            [_finite(value, "cycle_window_s") for value in window] for window in rows[5]
        ]
        measure_windows = [
            [_finite(value, "measure_cycle_window_s") for value in window]
            for window in rows[6]
        ]
        measure_names = [str(value) for value in rows[7]]
        reported_names = [str(value) for value in rows[8]]
        measure_values = [_finite(value, "cycle_measure_value") for value in rows[9]]
        reported_values = [
            _finite(value, "reported_cycle_measure_value") for value in rows[10]
        ]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("transient_generation_id") or "")
    topology_digest = str(contract.get("switching_topology_sha256") or "")
    measure_digest = str(contract.get("cycle_measure_table_sha256") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "initial_state_transient_generation_id",
                "topology_transient_generation_id",
                "accepted_grid_transient_generation_id",
                "cycle_window_transient_generation_id",
                "measure_transient_generation_id",
            )
        )
        and bool(state_ids)
        and all(state_ids)
        and len(set(state_ids)) == len(state_ids)
        and len(initial_state) == len(state_ids)
        and solver_state == initial_state
        and _is_sha256(topology_digest)
        and contract.get("solver_switching_topology_sha256") == topology_digest
        and len(accepted_grid) >= 2
        and accepted_grid[0] >= 0.0
        and all(right > left for left, right in zip(accepted_grid, accepted_grid[1:]))
        and measure_grid == accepted_grid
        and bool(cycle_windows)
        and all(
            len(window) == 2
            and accepted_grid[0] <= window[0] < window[1] <= accepted_grid[-1]
            for window in cycle_windows
        )
        and all(
            right[0] >= left[1]
            for left, right in zip(cycle_windows, cycle_windows[1:])
        )
        and measure_windows == cycle_windows
        and bool(measure_names)
        and all(measure_names)
        and len(set(measure_names)) == len(measure_names)
        and reported_names == measure_names
        and len(measure_values) == len(measure_names)
        and reported_values == measure_values
        and _is_sha256(measure_digest)
        and contract.get("reported_cycle_measure_table_sha256") == measure_digest
    )
