"""Solver-neutral identities for an ideal two-winding transformer."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed) or (positive and parsed <= 0.0):
        requirement = "positive and finite" if positive else "finite"
        raise ValueError(f"{name} must be {requirement}")
    return parsed


def _complex_pair(value: object, name: str) -> complex:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
    ):
        raise ValueError(f"{name} must contain [real, imag]")
    return complex(_finite(value[0], name), _finite(value[1], name))


def _relative_error(actual: complex | float, expected: complex | float) -> float:
    return float(abs(actual - expected) / max(abs(expected), 1.0e-300))


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


def _is_sha256(value: str) -> bool:
    return len(value) == 64 and all(
        character in "0123456789abcdef" for character in value
    )


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
        "exactly_four_timing_stages": timing_ok,
    }
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
