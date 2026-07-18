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


def _electrothermal_generation_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "electrothermal_temperature_device_model_network_timestep_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    keys = (
        "electrical_loss_time_s",
        "thermal_time_s",
        "electrical_loss_w",
        "thermal_input_loss_w",
        "junction_temperature_c",
        "reported_junction_temperature_c",
    )
    rows = tuple(contract.get(key) for key in keys)
    if not all(
        isinstance(values, Sequence) and not isinstance(values, (str, bytes))
        for values in rows
    ):
        return False
    try:
        model_temperature = _finite(
            contract.get("device_model_temperature_c"), "device_model_temperature_c"
        )
        solver_temperature = _finite(
            contract.get("solver_device_model_temperature_c"),
            "solver_device_model_temperature_c",
        )
        loss_time = [_finite(value, "electrical_loss_time_s") for value in rows[0]]
        thermal_time = [_finite(value, "thermal_time_s") for value in rows[1]]
        loss = [_finite(value, "electrical_loss_w") for value in rows[2]]
        thermal_loss = [_finite(value, "thermal_input_loss_w") for value in rows[3]]
        temperature = [_finite(value, "junction_temperature_c") for value in rows[4]]
        reported_temperature = [
            _finite(value, "reported_junction_temperature_c") for value in rows[5]
        ]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("electrothermal_generation_id") or "")
    model_digest = str(contract.get("device_model_sha256") or "")
    network_digest = str(contract.get("thermal_network_sha256") or "")
    table_digest = str(contract.get("electrothermal_table_sha256") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "device_model_electrothermal_generation_id",
                "temperature_electrothermal_generation_id",
                "loss_trace_electrothermal_generation_id",
                "thermal_network_electrothermal_generation_id",
                "time_grid_electrothermal_generation_id",
                "result_electrothermal_generation_id",
            )
        )
        and solver_temperature == model_temperature
        and _is_sha256(model_digest)
        and contract.get("solver_device_model_sha256") == model_digest
        and _is_sha256(network_digest)
        and contract.get("solver_thermal_network_sha256") == network_digest
        and len(loss_time) >= 2
        and loss_time[0] >= 0.0
        and all(right > left for left, right in zip(loss_time, loss_time[1:]))
        and thermal_time == loss_time
        and len(loss) == len(loss_time)
        and all(value >= 0.0 for value in loss)
        and thermal_loss == loss
        and len(temperature) == len(loss_time)
        and reported_temperature == temperature
        and _is_sha256(table_digest)
        and contract.get("reported_electrothermal_table_sha256") == table_digest
    )


def _monte_carlo_subcircuit_identity_ok(positive: Mapping[str, object]) -> bool:
    c = positive.get("subcircuit_monte_carlo_seed_model_include_raw_generation_identity")
    if c is None:
        return True
    if not isinstance(c, Mapping):
        return False
    try:
        sample_ids = [int(v) for v in c.get("sample_ids", [])]
        raw_ids = [int(v) for v in c.get("raw_sample_ids", [])]
        seeds = [int(v) for v in c.get("random_seeds", [])]
        raw_seeds = [int(v) for v in c.get("raw_random_seeds", [])]
    except (TypeError, ValueError):
        return False
    g = str(c.get("monte_carlo_generation_id") or "")
    model = str(c.get("model_include_sha256") or "")
    params = str(c.get("parameter_override_sha256") or "")
    table = str(c.get("raw_sample_table_sha256") or "")
    return (
        bool(g)
        and all(c.get(k) == g for k in (
            "seed_monte_carlo_generation_id", "model_include_monte_carlo_generation_id",
            "parameter_override_monte_carlo_generation_id", "raw_sample_monte_carlo_generation_id",
            "statistic_monte_carlo_generation_id"))
        and bool(sample_ids) and len(set(sample_ids)) == len(sample_ids) and raw_ids == sample_ids
        and len(seeds) == len(sample_ids) and len(set(seeds)) == len(seeds) and raw_seeds == seeds
        and _is_sha256(model) and c.get("raw_model_include_sha256") == model
        and _is_sha256(params) and c.get("raw_parameter_override_sha256") == params
        and _is_sha256(table) and c.get("statistic_sample_table_sha256") == table
    )


def _behavioral_switch_hysteresis_identity_ok(positive: Mapping[str, object]) -> bool:
    c = positive.get("behavioral_switch_hysteresis_state_timestep_measure_generation_identity")
    if c is None:
        return True
    if not isinstance(c, Mapping):
        return False
    try:
        states = [int(v) for v in c.get("hysteresis_states", [])]
        result_states = [int(v) for v in c.get("measure_hysteresis_states", [])]
        events = [float(v) for v in c.get("event_times_s", [])]
        result_events = [float(v) for v in c.get("measure_event_times_s", [])]
        times = [float(v) for v in c.get("accepted_time_s", [])]
        result_times = [float(v) for v in c.get("measure_time_s", [])]
        window = [float(v) for v in c.get("measure_window_s", [])]
        result_window = [float(v) for v in c.get("reported_measure_window_s", [])]
    except (TypeError, ValueError):
        return False
    g = str(c.get("transient_generation_id") or "")
    digest = str(c.get("measure_table_sha256") or "")
    return (
        bool(g)
        and all(c.get(k) == g for k in (
            "hysteresis_state_transient_generation_id", "event_history_transient_generation_id",
            "accepted_timestep_transient_generation_id", "measure_window_transient_generation_id",
            "measure_result_transient_generation_id"))
        and bool(states) and all(v in {0, 1} for v in states) and result_states == states
        and len(events) == len(states) and all(math.isfinite(v) for v in events)
        and all(a < b for a, b in zip(events, events[1:])) and result_events == events
        and len(times) >= 2 and all(math.isfinite(v) for v in times)
        and all(a < b for a, b in zip(times, times[1:])) and result_times == times
        and len(window) == 2 and times[0] <= window[0] < window[1] <= times[-1]
        and result_window == window and _is_sha256(digest)
        and c.get("reported_measure_table_sha256") == digest
    )


def _hierarchical_step_identity_ok(positive: Mapping[str, object]) -> bool:
    c = positive.get(
        "hierarchical_step_parameter_scope_model_bin_temperature_sample_generation_identity"
    )
    if c is None:
        return True
    if not isinstance(c, Mapping):
        return False
    generation = str(c.get("step_generation_id") or "")
    try:
        sample_ids = [int(value) for value in c.get("step_sample_ids", [])]
        result_ids = [int(value) for value in c.get("result_step_sample_ids", [])]
        parameters = [float(value) for value in c.get("parameter_values_ohm", [])]
        result_parameters = [
            float(value) for value in c.get("result_parameter_values_ohm", [])
        ]
        temperatures = [float(value) for value in c.get("temperatures_c", [])]
        result_temperatures = [
            float(value) for value in c.get("result_temperatures_c", [])
        ]
        samples = [float(value) for value in c.get("sample_values_v", [])]
        result_samples = [float(value) for value in c.get("result_sample_values_v", [])]
    except (TypeError, ValueError):
        return False
    hierarchy = str(c.get("hierarchy_path") or "")
    scope = str(c.get("parameter_scope") or "")
    model_bin = str(c.get("model_bin") or "")
    digest = str(c.get("step_table_sha256") or "")
    return (
        bool(generation)
        and all(
            c.get(key) == generation
            for key in (
                "scope_step_generation_id",
                "model_bin_step_generation_id",
                "temperature_step_generation_id",
                "sample_row_step_generation_id",
                "result_step_generation_id",
            )
        )
        and bool(hierarchy)
        and c.get("result_hierarchy_path") == hierarchy
        and bool(scope)
        and c.get("result_parameter_scope") == scope
        and bool(model_bin)
        and c.get("result_model_bin") == model_bin
        and bool(sample_ids)
        and all(value > 0 for value in sample_ids)
        and len(set(sample_ids)) == len(sample_ids)
        and result_ids == sample_ids
        and len(parameters) == len(sample_ids)
        and all(math.isfinite(value) and value > 0.0 for value in parameters)
        and result_parameters == parameters
        and len(temperatures) == len(sample_ids)
        and all(math.isfinite(value) and value >= -273.15 for value in temperatures)
        and result_temperatures == temperatures
        and len(samples) == len(sample_ids)
        and all(math.isfinite(value) for value in samples)
        and result_samples == samples
        and _is_sha256(digest)
        and c.get("result_step_table_sha256") == digest
    )


def _ac_noise_source_identity_ok(positive: Mapping[str, object]) -> bool:
    c = positive.get(
        "ac_noise_source_normalization_node_alias_complex_axis_generation_identity"
    )
    if c is None:
        return True
    if not isinstance(c, Mapping):
        return False
    generation = str(c.get("analysis_generation_id") or "")
    aliases = c.get("node_aliases")
    result_aliases = c.get("result_node_aliases")
    transfer = c.get("transfer_function_ri")
    result_transfer = c.get("result_transfer_function_ri")
    if not all(
        isinstance(value, Sequence) and not isinstance(value, (str, bytes))
        for value in (aliases, result_aliases, transfer, result_transfer)
    ):
        return False
    try:
        frequency = [float(value) for value in c.get("frequency_grid_hz", [])]
        result_frequency = [
            float(value) for value in c.get("result_frequency_grid_hz", [])
        ]
        transfer_pairs = [[float(value) for value in row] for row in transfer]
        result_pairs = [[float(value) for value in row] for row in result_transfer]
        noise = [float(value) for value in c.get("output_noise_v_per_sqrt_hz", [])]
        result_noise = [
            float(value) for value in c.get("result_output_noise_v_per_sqrt_hz", [])
        ]
    except (TypeError, ValueError):
        return False
    alias_rows = [[str(value) for value in row] for row in aliases]
    result_alias_rows = [[str(value) for value in row] for row in result_aliases]
    source = str(c.get("source_id") or "")
    normalization = str(c.get("source_normalization") or "")
    digest = str(c.get("ac_noise_table_sha256") or "")
    return (
        bool(generation)
        and all(
            c.get(key) == generation
            for key in (
                "source_analysis_generation_id",
                "node_alias_analysis_generation_id",
                "complex_axis_analysis_generation_id",
                "frequency_grid_analysis_generation_id",
                "result_analysis_generation_id",
            )
        )
        and bool(source)
        and c.get("result_source_id") == source
        and normalization == "1_V_ac"
        and c.get("result_source_normalization") == normalization
        and bool(alias_rows)
        and all(len(row) == 2 and all(row) for row in alias_rows)
        and len({row[0] for row in alias_rows}) == len(alias_rows)
        and result_alias_rows == alias_rows
        and c.get("complex_axis_convention") == "real_imaginary"
        and c.get("result_complex_axis_convention") == "real_imaginary"
        and len(frequency) >= 3
        and all(math.isfinite(value) and value > 0.0 for value in frequency)
        and all(right > left for left, right in zip(frequency, frequency[1:]))
        and result_frequency == frequency
        and len(transfer_pairs) == len(frequency)
        and all(
            len(row) == 2 and all(math.isfinite(value) for value in row)
            for row in transfer_pairs
        )
        and result_pairs == transfer_pairs
        and len(noise) == len(frequency)
        and all(math.isfinite(value) and value >= 0.0 for value in noise)
        and result_noise == noise
        and _is_sha256(digest)
        and c.get("result_ac_noise_table_sha256") == digest
    )


def _transient_startup_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "transient_startup_initial_condition_uic_operating_point_waveform_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    initial = contract.get("initial_conditions")
    result_initial = contract.get("result_initial_conditions")
    if not all(
        isinstance(value, Sequence) and not isinstance(value, (str, bytes))
        for value in (initial, result_initial)
    ):
        return False
    try:
        initial_rows = [[str(row[0]), float(row[1])] for row in initial]
        result_initial_rows = [[str(row[0]), float(row[1])] for row in result_initial]
        time_s = [float(value) for value in contract.get("accepted_time_s", [])]
        result_time_s = [float(value) for value in contract.get("result_time_s", [])]
        traces = [str(value) for value in contract.get("waveform_trace_ids", [])]
        result_traces = [
            str(value) for value in contract.get("result_waveform_trace_ids", [])
        ]
    except (IndexError, TypeError, ValueError):
        return False
    generation = str(contract.get("transient_generation_id") or "")
    operating_point_digest = str(contract.get("operating_point_sha256") or "")
    waveform_digest = str(contract.get("waveform_table_sha256") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "startup_mode_transient_generation_id",
                "initial_condition_transient_generation_id",
                "uic_transient_generation_id",
                "operating_point_transient_generation_id",
                "accepted_time_grid_transient_generation_id",
                "waveform_transient_generation_id",
                "result_transient_generation_id",
            )
        )
        and contract.get("startup_mode") == "operating_point_then_transient"
        and contract.get("result_startup_mode") == contract.get("startup_mode")
        and contract.get("uic_enabled") is False
        and contract.get("result_uic_enabled") is False
        and bool(initial_rows)
        and all(name and math.isfinite(value) for name, value in initial_rows)
        and len({name for name, _ in initial_rows}) == len(initial_rows)
        and result_initial_rows == initial_rows
        and _is_sha256(operating_point_digest)
        and contract.get("result_operating_point_sha256") == operating_point_digest
        and len(time_s) >= 2
        and time_s[0] >= 0.0
        and all(math.isfinite(value) for value in time_s)
        and all(right > left for left, right in zip(time_s, time_s[1:]))
        and result_time_s == time_s
        and bool(traces)
        and all(traces)
        and len(set(traces)) == len(traces)
        and result_traces == traces
        and _is_sha256(waveform_digest)
        and contract.get("result_waveform_table_sha256") == waveform_digest
    )


def _stepped_monte_carlo_aggregation_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "stepped_monte_carlo_measure_aggregation_failed_row_seed_weight_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        sample_ids = [int(value) for value in contract.get("sample_ids", [])]
        result_ids = [int(value) for value in contract.get("result_sample_ids", [])]
        seeds = [int(value) for value in contract.get("random_seeds", [])]
        result_seeds = [int(value) for value in contract.get("result_random_seeds", [])]
        statuses = [str(value) for value in contract.get("measure_statuses", [])]
        accepted_ids = [
            int(value) for value in contract.get("accepted_sample_ids", [])
        ]
        result_accepted_ids = [
            int(value) for value in contract.get("result_accepted_sample_ids", [])
        ]
        failed_ids = [int(value) for value in contract.get("failed_sample_ids", [])]
        result_failed_ids = [
            int(value) for value in contract.get("result_failed_sample_ids", [])
        ]
        values = [float(value) for value in contract.get("sample_values", [])]
        weights = [float(value) for value in contract.get("sample_weights", [])]
        result_weights = [
            float(value) for value in contract.get("result_sample_weights", [])
        ]
        reported_mean = float(contract.get("reported_weighted_mean"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("monte_carlo_generation_id") or "")
    expected_accepted = [
        sample_id
        for sample_id, status in zip(sample_ids, statuses)
        if status == "passed"
    ]
    expected_failed = [
        sample_id
        for sample_id, status in zip(sample_ids, statuses)
        if status == "failed"
    ]
    total_weight = sum(weights)
    weighted_mean = (
        sum(value * weight for value, weight in zip(values, weights)) / total_weight
        if total_weight > 0.0
        else math.nan
    )
    digest = str(contract.get("sample_table_sha256") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "measure_row_monte_carlo_generation_id",
                "seed_monte_carlo_generation_id",
                "filter_monte_carlo_generation_id",
                "weight_monte_carlo_generation_id",
                "aggregation_monte_carlo_generation_id",
                "result_monte_carlo_generation_id",
            )
        )
        and bool(sample_ids)
        and all(value > 0 for value in sample_ids)
        and len(set(sample_ids)) == len(sample_ids)
        and result_ids == sample_ids
        and len(seeds) == len(sample_ids)
        and len(set(seeds)) == len(seeds)
        and result_seeds == seeds
        and len(statuses) == len(sample_ids)
        and all(status in {"passed", "failed"} for status in statuses)
        and accepted_ids == expected_accepted
        and result_accepted_ids == accepted_ids
        and failed_ids == expected_failed
        and result_failed_ids == failed_ids
        and len(values) == len(weights) == len(accepted_ids) == len(sample_ids)
        and all(math.isfinite(value) for value in values)
        and all(math.isfinite(weight) and weight > 0.0 for weight in weights)
        and result_weights == weights
        and contract.get("aggregation_rule") == "weighted_mean"
        and contract.get("result_aggregation_rule") == contract.get("aggregation_rule")
        and math.isfinite(reported_mean)
        and math.isclose(reported_mean, weighted_mean, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and _is_sha256(digest)
        and contract.get("result_sample_table_sha256") == digest
    )


def _ac_sweep_measure_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "ac_sweep_mode_frequency_grid_complex_phase_unwrap_measure_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        points_per_decade = int(contract.get("points_per_decade"))
        result_points_per_decade = int(contract.get("result_points_per_decade"))
        frequencies = [float(value) for value in contract.get("frequency_grid_hz", [])]
        result_frequencies = [
            float(value) for value in contract.get("result_frequency_grid_hz", [])
        ]
        row_ids = [str(value) for value in contract.get("measure_row_ids", [])]
        result_row_ids = [
            str(value) for value in contract.get("result_measure_row_ids", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("ac_generation_id") or "")
    digest = str(contract.get("measure_table_sha256") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "sweep_mode_ac_generation_id",
                "frequency_grid_ac_generation_id",
                "complex_basis_ac_generation_id",
                "phase_unwrap_ac_generation_id",
                "measure_row_ac_generation_id",
                "result_ac_generation_id",
            )
        )
        and contract.get("sweep_mode") == "decade"
        and contract.get("result_sweep_mode") == contract.get("sweep_mode")
        and points_per_decade > 0
        and result_points_per_decade == points_per_decade
        and len(frequencies) >= 3
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and contract.get("complex_basis") == "real_imaginary"
        and contract.get("result_complex_basis") == contract.get("complex_basis")
        and contract.get("phase_unwrap") == "continuous_radians"
        and contract.get("result_phase_unwrap") == contract.get("phase_unwrap")
        and bool(row_ids)
        and all(row_ids)
        and len(set(row_ids)) == len(row_ids)
        and result_row_ids == row_ids
        and _is_sha256(digest)
        and contract.get("result_measure_table_sha256") == digest
    )


def _electrothermal_waveform_closure_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "electrothermal_device_power_temperature_model_thermal_network_timestep_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        trace_ids = [str(value) for value in contract.get("device_power_trace_ids", [])]
        result_trace_ids = [
            str(value) for value in contract.get("result_device_power_trace_ids", [])
        ]
        timestep = float(contract.get("time_step_s"))
        result_timestep = float(contract.get("result_time_step_s"))
        time_grid = [float(value) for value in contract.get("time_grid_s", [])]
        result_time_grid = [
            float(value) for value in contract.get("result_time_grid_s", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("electrothermal_generation_id") or "")
    digests = (
        ("device_power_sha256", "result_device_power_sha256"),
        ("temperature_model_sha256", "result_temperature_model_sha256"),
        ("thermal_network_sha256", "result_thermal_network_sha256"),
        ("temperature_waveform_sha256", "result_temperature_waveform_sha256"),
    )
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "device_power_electrothermal_generation_id",
                "temperature_model_electrothermal_generation_id",
                "thermal_network_electrothermal_generation_id",
                "timestep_electrothermal_generation_id",
                "result_electrothermal_generation_id",
            )
        )
        and bool(trace_ids)
        and all(trace_ids)
        and len(set(trace_ids)) == len(trace_ids)
        and result_trace_ids == trace_ids
        and bool(str(contract.get("temperature_model_id") or ""))
        and contract.get("result_temperature_model_id")
        == contract.get("temperature_model_id")
        and bool(str(contract.get("thermal_network_id") or ""))
        and contract.get("result_thermal_network_id")
        == contract.get("thermal_network_id")
        and math.isfinite(timestep)
        and timestep > 0.0
        and result_timestep == timestep
        and len(time_grid) >= 2
        and time_grid[0] == 0.0
        and all(math.isfinite(value) for value in time_grid)
        and all(right > left for left, right in zip(time_grid, time_grid[1:]))
        and all(
            math.isclose(
                right - left,
                timestep,
                rel_tol=1.0e-12,
                abs_tol=1.0e-15,
            )
            for left, right in zip(time_grid, time_grid[1:])
        )
        and result_time_grid == time_grid
        and all(
            _is_sha256(str(contract.get(source) or ""))
            and contract.get(target) == contract.get(source)
            for source, target in digests
        )
    )


def _noise_integration_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "noise_input_output_source_normalization_psd_sidedness_integration_grid_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        frequencies = [float(value) for value in contract.get("frequency_grid_hz", [])]
        result_frequencies = [
            float(value) for value in contract.get("result_frequency_grid_hz", [])
        ]
        psd_values = [float(value) for value in contract.get("psd_values", [])]
        result_psd_values = [
            float(value) for value in contract.get("result_psd_values", [])
        ]
        integrated = float(contract.get("integrated_noise_v_rms"))
        result_integrated = float(contract.get("result_integrated_noise_v_rms"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("noise_generation_id") or "")
    digest = str(contract.get("noise_result_sha256") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "input_source_noise_generation_id",
                "output_source_noise_generation_id",
                "normalization_noise_generation_id",
                "psd_noise_generation_id",
                "integration_grid_noise_generation_id",
                "result_noise_generation_id",
            )
        )
        and bool(str(contract.get("input_source_id") or ""))
        and contract.get("result_input_source_id") == contract.get("input_source_id")
        and bool(str(contract.get("input_node") or ""))
        and contract.get("result_input_node") == contract.get("input_node")
        and bool(str(contract.get("output_node") or ""))
        and contract.get("result_output_node") == contract.get("output_node")
        and contract.get("normalization") == "input_referred_voltage_density"
        and contract.get("result_normalization") == contract.get("normalization")
        and contract.get("psd_sidedness") == "one_sided"
        and contract.get("result_psd_sidedness") == contract.get("psd_sidedness")
        and contract.get("psd_unit") == "V^2/Hz"
        and contract.get("result_psd_unit") == contract.get("psd_unit")
        and len(frequencies) >= 3
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and len(psd_values) == len(frequencies)
        and all(math.isfinite(value) and value >= 0.0 for value in psd_values)
        and result_psd_values == psd_values
        and contract.get("integration_rule") == "log_frequency_trapezoid"
        and contract.get("result_integration_rule") == contract.get("integration_rule")
        and math.isfinite(integrated)
        and integrated >= 0.0
        and result_integrated == integrated
        and _is_sha256(digest)
        and contract.get("accepted_noise_result_sha256") == digest
    )


def _switch_event_timing_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "switch_hysteresis_event_order_max_timestep_measure_window_waveform_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        high = float(contract.get("threshold_high_v"))
        result_high = float(contract.get("result_threshold_high_v"))
        low = float(contract.get("threshold_low_v"))
        result_low = float(contract.get("result_threshold_low_v"))
        max_timestep = float(contract.get("max_timestep_s"))
        result_max_timestep = float(contract.get("result_max_timestep_s"))
        window = [float(value) for value in contract.get("measure_window_s", [])]
        result_window = [
            float(value) for value in contract.get("result_measure_window_s", [])
        ]
        events = [str(value) for value in contract.get("event_order", [])]
        result_events = [str(value) for value in contract.get("result_event_order", [])]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("switch_generation_id") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "hysteresis_switch_generation_id",
                "event_order_switch_generation_id",
                "timestep_switch_generation_id",
                "measure_window_switch_generation_id",
                "waveform_switch_generation_id",
                "result_switch_generation_id",
            )
        )
        and contract.get("switch_model_id") == "voltage_hysteretic_switch"
        and contract.get("result_switch_model_id") == contract.get("switch_model_id")
        and all(math.isfinite(value) for value in (high, low))
        and high > low
        and result_high == high
        and result_low == low
        and events == ["rising_on", "falling_off", "rising_on"]
        and result_events == events
        and math.isfinite(max_timestep)
        and max_timestep > 0.0
        and result_max_timestep == max_timestep
        and len(window) == 2
        and all(math.isfinite(value) for value in window)
        and 0.0 <= window[0] < window[1]
        and max_timestep <= (window[1] - window[0]) / 10.0
        and result_window == window
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("result_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("measure_table_sha256") or ""))
        and contract.get("accepted_measure_table_sha256")
        == contract.get("measure_table_sha256")
    )


def _smps_efficiency_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "smps_efficiency_source_load_steady_cycle_energy_integration_switching_waveform_timestep_result_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        window = [float(value) for value in contract.get("steady_cycle_window_s", [])]
        result_window = [
            float(value)
            for value in contract.get("result_steady_cycle_window_s", [])
        ]
        source_energy = float(contract.get("source_energy_j"))
        result_source_energy = float(contract.get("result_source_energy_j"))
        load_energy = float(contract.get("load_energy_j"))
        result_load_energy = float(contract.get("result_load_energy_j"))
        efficiency = float(contract.get("efficiency"))
        result_efficiency = float(contract.get("result_efficiency"))
        timestep = float(contract.get("max_timestep_s"))
        result_timestep = float(contract.get("result_max_timestep_s"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("efficiency_generation_id") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "source_trace_efficiency_generation_id",
                "load_trace_efficiency_generation_id",
                "cycle_window_efficiency_generation_id",
                "integration_efficiency_generation_id",
                "waveform_efficiency_generation_id",
                "timestep_efficiency_generation_id",
                "result_efficiency_generation_id",
            )
        )
        and bool(str(contract.get("source_trace_id") or ""))
        and contract.get("result_source_trace_id") == contract.get("source_trace_id")
        and bool(str(contract.get("load_trace_id") or ""))
        and contract.get("result_load_trace_id") == contract.get("load_trace_id")
        and len(window) == 2
        and all(math.isfinite(value) for value in window)
        and 0.0 <= window[0] < window[1]
        and result_window == window
        and contract.get("energy_integration_rule")
        == "trapezoid_power_over_time"
        and contract.get("result_energy_integration_rule")
        == contract.get("energy_integration_rule")
        and all(
            math.isfinite(value)
            for value in (
                source_energy,
                result_source_energy,
                load_energy,
                result_load_energy,
                efficiency,
                result_efficiency,
                timestep,
                result_timestep,
            )
        )
        and source_energy > 0.0
        and 0.0 <= load_energy <= source_energy
        and result_source_energy == source_energy
        and result_load_energy == load_energy
        and 0.0 <= efficiency <= 1.0
        and math.isclose(
            efficiency,
            load_energy / source_energy,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and result_efficiency == efficiency
        and timestep > 0.0
        and timestep <= (window[1] - window[0]) / 20.0
        and result_timestep == timestep
        and _is_sha256(str(contract.get("switching_waveform_sha256") or ""))
        and contract.get("result_switching_waveform_sha256")
        == contract.get("switching_waveform_sha256")
        and _is_sha256(str(contract.get("efficiency_result_sha256") or ""))
        and contract.get("accepted_efficiency_result_sha256")
        == contract.get("efficiency_result_sha256")
    )


def _loop_gain_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "loop_gain_break_injection_sign_phase_unwrap_crossover_margin_frequency_grid_result_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        nodes = [str(value) for value in contract.get("loop_break_nodes", [])]
        result_nodes = [
            str(value) for value in contract.get("result_loop_break_nodes", [])
        ]
        injection_sign = int(contract.get("injection_sign"))
        result_injection_sign = int(contract.get("result_injection_sign"))
        frequencies = [float(value) for value in contract.get("frequency_grid_hz", [])]
        result_frequencies = [
            float(value) for value in contract.get("result_frequency_grid_hz", [])
        ]
        gains = [float(value) for value in contract.get("loop_gain_db", [])]
        result_gains = [float(value) for value in contract.get("result_loop_gain_db", [])]
        phases = [float(value) for value in contract.get("phase_deg", [])]
        result_phases = [float(value) for value in contract.get("result_phase_deg", [])]
        crossover = float(contract.get("gain_crossover_hz"))
        result_crossover = float(contract.get("result_gain_crossover_hz"))
        phase_margin = float(contract.get("phase_margin_deg"))
        result_phase_margin = float(contract.get("result_phase_margin_deg"))
        gain_margin = float(contract.get("gain_margin_db"))
        result_gain_margin = float(contract.get("result_gain_margin_db"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("loop_gain_generation_id") or "")
    zero_gain_indices = [
        index for index, gain in enumerate(gains) if math.isclose(gain, 0.0, abs_tol=1.0e-12)
    ]
    crossover_index_ok = (
        len(zero_gain_indices) == 1
        and frequencies[zero_gain_indices[0]] == crossover
        and math.isclose(
            180.0 + phases[zero_gain_indices[0]],
            phase_margin,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
    ) if frequencies and len(frequencies) == len(gains) == len(phases) else False
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "break_loop_gain_generation_id",
                "injection_loop_gain_generation_id",
                "phase_loop_gain_generation_id",
                "crossover_loop_gain_generation_id",
                "frequency_loop_gain_generation_id",
                "result_loop_gain_generation_id",
            )
        )
        and bool(str(contract.get("loop_break_element") or ""))
        and contract.get("result_loop_break_element")
        == contract.get("loop_break_element")
        and len(nodes) == 2
        and all(nodes)
        and nodes[0] != nodes[1]
        and result_nodes == nodes
        and injection_sign in {-1, 1}
        and result_injection_sign == injection_sign
        and contract.get("phase_unwrap_rule") == "continuous_negative_180"
        and contract.get("result_phase_unwrap_rule")
        == contract.get("phase_unwrap_rule")
        and contract.get("crossover_interpolation")
        == "log_frequency_linear_db"
        and contract.get("result_crossover_interpolation")
        == contract.get("crossover_interpolation")
        and len(frequencies) >= 3
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and len(gains) == len(frequencies)
        and len(phases) == len(frequencies)
        and all(math.isfinite(value) for value in gains + phases)
        and result_gains == gains
        and result_phases == phases
        and crossover_index_ok
        and result_crossover == crossover
        and math.isfinite(phase_margin)
        and 0.0 < phase_margin < 180.0
        and result_phase_margin == phase_margin
        and math.isfinite(gain_margin)
        and gain_margin > 0.0
        and result_gain_margin == gain_margin
        and _is_sha256(str(contract.get("loop_gain_result_sha256") or ""))
        and contract.get("accepted_loop_gain_result_sha256")
        == contract.get("loop_gain_result_sha256")
    )


def _mosfet_switching_loss_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "mosfet_switching_loss_gate_charge_overlap_deadtime_event_grid_temperature_cycle_result_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        deadtime = float(contract.get("deadtime_s"))
        result_deadtime = float(contract.get("result_deadtime_s"))
        events = [float(value) for value in contract.get("event_times_s", [])]
        result_events = [
            float(value) for value in contract.get("result_event_times_s", [])
        ]
        temperature = float(contract.get("junction_temperature_c"))
        result_temperature = float(contract.get("result_junction_temperature_c"))
        window = [float(value) for value in contract.get("cycle_window_s", [])]
        result_window = [
            float(value) for value in contract.get("result_cycle_window_s", [])
        ]
        turn_on = float(contract.get("turn_on_energy_j"))
        result_turn_on = float(contract.get("result_turn_on_energy_j"))
        turn_off = float(contract.get("turn_off_energy_j"))
        result_turn_off = float(contract.get("result_turn_off_energy_j"))
        frequency = float(contract.get("switching_frequency_hz"))
        result_frequency = float(contract.get("result_switching_frequency_hz"))
        loss = float(contract.get("switching_loss_w"))
        result_loss = float(contract.get("result_switching_loss_w"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("switching_generation_id") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "gate_charge_switching_generation_id",
                "overlap_switching_generation_id",
                "deadtime_switching_generation_id",
                "event_grid_switching_generation_id",
                "temperature_switching_generation_id",
                "cycle_switching_generation_id",
                "result_switching_generation_id",
            )
        )
        and contract.get("gate_charge_trace_id") == "Qgate(M1)"
        and contract.get("result_gate_charge_trace_id")
        == contract.get("gate_charge_trace_id")
        and contract.get("overlap_power_trace_id") == "Vds(M1)*Id(M1)"
        and contract.get("result_overlap_power_trace_id")
        == contract.get("overlap_power_trace_id")
        and math.isfinite(deadtime)
        and deadtime > 0.0
        and result_deadtime == deadtime
        and len(events) >= 4
        and all(math.isfinite(value) for value in events)
        and all(right > left for left, right in zip(events, events[1:]))
        and result_events == events
        and contract.get("event_grid_rule") == "edge-aligned-local-refinement"
        and contract.get("result_event_grid_rule") == contract.get("event_grid_rule")
        and math.isfinite(temperature)
        and temperature > -273.15
        and result_temperature == temperature
        and len(window) == 2
        and all(math.isfinite(value) for value in window)
        and 0.0 <= window[0] < window[1]
        and all(window[0] <= event <= window[1] for event in events)
        and result_window == window
        and all(
            math.isfinite(value) and value >= 0.0
            for value in (turn_on, turn_off)
        )
        and result_turn_on == turn_on
        and result_turn_off == turn_off
        and math.isfinite(frequency)
        and frequency > 0.0
        and result_frequency == frequency
        and math.isfinite(loss)
        and loss >= 0.0
        and math.isclose(
            loss, (turn_on + turn_off) * frequency, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and result_loss == loss
        and _is_sha256(str(contract.get("event_grid_sha256") or ""))
        and contract.get("result_event_grid_sha256")
        == contract.get("event_grid_sha256")
        and _is_sha256(str(contract.get("switching_waveform_sha256") or ""))
        and contract.get("result_switching_waveform_sha256")
        == contract.get("switching_waveform_sha256")
        and _is_sha256(str(contract.get("switching_loss_result_sha256") or ""))
        and contract.get("accepted_switching_loss_result_sha256")
        == contract.get("switching_loss_result_sha256")
    )


def _step_response_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "step_response_initial_final_rise_threshold_settling_band_overshoot_window_waveform_result_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        initial = float(contract.get("initial_value"))
        result_initial = float(contract.get("result_initial_value"))
        final = float(contract.get("final_value"))
        result_final = float(contract.get("result_final_value"))
        thresholds = [
            float(value) for value in contract.get("rise_threshold_fractions", [])
        ]
        result_thresholds = [
            float(value)
            for value in contract.get("result_rise_threshold_fractions", [])
        ]
        crossings = [
            float(value) for value in contract.get("rise_crossing_times_s", [])
        ]
        result_crossings = [
            float(value) for value in contract.get("result_rise_crossing_times_s", [])
        ]
        rise_time = float(contract.get("rise_time_s"))
        result_rise_time = float(contract.get("result_rise_time_s"))
        settling_band = float(contract.get("settling_band_fraction"))
        result_settling_band = float(contract.get("result_settling_band_fraction"))
        settling_time = float(contract.get("settling_time_s"))
        result_settling_time = float(contract.get("result_settling_time_s"))
        peak = float(contract.get("overshoot_peak"))
        result_peak = float(contract.get("result_overshoot_peak"))
        overshoot = float(contract.get("overshoot_fraction"))
        result_overshoot = float(contract.get("result_overshoot_fraction"))
        window = [float(value) for value in contract.get("measurement_window_s", [])]
        result_window = [
            float(value) for value in contract.get("result_measurement_window_s", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("step_generation_id") or "")
    expected_overshoot = (peak - final) / abs(final - initial)
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "initial_step_generation_id",
                "final_step_generation_id",
                "rise_step_generation_id",
                "settling_step_generation_id",
                "overshoot_step_generation_id",
                "window_step_generation_id",
                "waveform_step_generation_id",
                "result_step_generation_id",
            )
        )
        and all(math.isfinite(value) for value in (initial, final))
        and final != initial
        and result_initial == initial
        and result_final == final
        and len(thresholds) == 2
        and 0.0 < thresholds[0] < thresholds[1] < 1.0
        and result_thresholds == thresholds
        and len(crossings) == 2
        and all(math.isfinite(value) for value in crossings)
        and 0.0 <= crossings[0] < crossings[1]
        and result_crossings == crossings
        and math.isfinite(rise_time)
        and math.isclose(
            rise_time, crossings[1] - crossings[0], rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and result_rise_time == rise_time
        and math.isfinite(settling_band)
        and 0.0 < settling_band < 1.0
        and result_settling_band == settling_band
        and math.isfinite(settling_time)
        and settling_time >= crossings[1]
        and result_settling_time == settling_time
        and math.isfinite(peak)
        and result_peak == peak
        and math.isfinite(overshoot)
        and overshoot >= 0.0
        and math.isclose(
            overshoot, expected_overshoot, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and result_overshoot == overshoot
        and len(window) == 2
        and all(math.isfinite(value) for value in window)
        and 0.0 <= window[0] < window[1]
        and window[0] <= crossings[0] < crossings[1] <= settling_time <= window[1]
        and result_window == window
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("result_waveform_sha256")
        == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("step_result_sha256") or ""))
        and contract.get("accepted_step_result_sha256")
        == contract.get("step_result_sha256")
    )


def _mosfet_soa_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "mosfet_soa_vds_id_pulse_width_duty_temperature_model_waveform_result_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        vds = float(contract.get("vds_v"))
        result_vds = float(contract.get("result_vds_v"))
        drain_current = float(contract.get("id_a"))
        result_drain_current = float(contract.get("result_id_a"))
        pulse_width = float(contract.get("pulse_width_s"))
        result_pulse_width = float(contract.get("result_pulse_width_s"))
        period = float(contract.get("repetition_period_s"))
        result_period = float(contract.get("result_repetition_period_s"))
        duty = float(contract.get("duty_cycle"))
        result_duty = float(contract.get("result_duty_cycle"))
        temperature = float(contract.get("junction_temperature_c"))
        result_temperature = float(contract.get("result_junction_temperature_c"))
        current_limit = float(contract.get("soa_limit_id_a"))
        result_current_limit = float(contract.get("result_soa_limit_id_a"))
        margin = float(contract.get("soa_margin_fraction"))
        result_margin = float(contract.get("result_soa_margin_fraction"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("soa_generation_id") or "")
    expected_duty = pulse_width / period
    expected_margin = (current_limit - drain_current) / current_limit
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "voltage_soa_generation_id",
                "current_soa_generation_id",
                "pulse_soa_generation_id",
                "duty_soa_generation_id",
                "temperature_soa_generation_id",
                "model_soa_generation_id",
                "waveform_soa_generation_id",
                "result_soa_generation_id",
            )
        )
        and all(
            math.isfinite(value) and value > 0.0
            for value in (vds, drain_current, pulse_width, period, current_limit)
        )
        and pulse_width < period
        and result_vds == vds
        and result_drain_current == drain_current
        and result_pulse_width == pulse_width
        and result_period == period
        and math.isfinite(duty)
        and 0.0 < duty < 1.0
        and math.isclose(duty, expected_duty, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and result_duty == duty
        and math.isfinite(temperature)
        and temperature > -273.15
        and result_temperature == temperature
        and drain_current <= current_limit
        and result_current_limit == current_limit
        and math.isfinite(margin)
        and margin >= 0.0
        and math.isclose(
            margin, expected_margin, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and result_margin == margin
        and _is_sha256(str(contract.get("model_card_sha256") or ""))
        and contract.get("result_model_card_sha256")
        == contract.get("model_card_sha256")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("result_waveform_sha256")
        == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("soa_result_sha256") or ""))
        and contract.get("accepted_soa_result_sha256")
        == contract.get("soa_result_sha256")
    )


def _monte_carlo_yield_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "monte_carlo_yield_distribution_tolerance_seed_failure_sample_owner_result_generation_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        parameters = [str(value) for value in contract.get("parameter_order", [])]
        result_parameters = [
            str(value) for value in contract.get("result_parameter_order", [])
        ]
        families = [
            str(value) for value in contract.get("distribution_families", [])
        ]
        result_families = [
            str(value)
            for value in contract.get("result_distribution_families", [])
        ]
        nominal = [float(value) for value in contract.get("nominal_values", [])]
        result_nominal = [
            float(value) for value in contract.get("result_nominal_values", [])
        ]
        tolerances = [
            float(value) for value in contract.get("relative_tolerances", [])
        ]
        result_tolerances = [
            float(value)
            for value in contract.get("result_relative_tolerances", [])
        ]
        seeds = [int(value) for value in contract.get("seed_schedule", [])]
        result_seeds = [
            int(value) for value in contract.get("result_seed_schedule", [])
        ]
        sample_ids = [int(value) for value in contract.get("sample_ids", [])]
        result_sample_ids = [
            int(value) for value in contract.get("result_sample_ids", [])
        ]
        failed_ids = [
            int(value) for value in contract.get("failed_sample_ids", [])
        ]
        result_failed_ids = [
            int(value) for value in contract.get("result_failed_sample_ids", [])
        ]
        accepted_ids = [
            int(value) for value in contract.get("accepted_sample_ids", [])
        ]
        result_accepted_ids = [
            int(value) for value in contract.get("result_accepted_sample_ids", [])
        ]
        yield_fraction = float(contract.get("yield_fraction"))
        result_yield_fraction = float(contract.get("result_yield_fraction"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("yield_generation_id") or "")
    expected_accepted = [
        sample_id for sample_id in sample_ids if sample_id not in set(failed_ids)
    ]
    expected_yield = len(expected_accepted) / len(sample_ids) if sample_ids else -1.0
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "distribution_yield_generation_id",
                "tolerance_yield_generation_id",
                "seed_yield_generation_id",
                "criterion_yield_generation_id",
                "sample_yield_generation_id",
                "owner_yield_generation_id",
                "result_yield_generation_id",
            )
        )
        and bool(parameters)
        and all(parameters)
        and len(set(parameters)) == len(parameters)
        and result_parameters == parameters
        and len(families) == len(nominal) == len(tolerances) == len(parameters)
        and all(family in {"gaussian", "uniform"} for family in families)
        and result_families == families
        and all(math.isfinite(value) and value > 0.0 for value in nominal)
        and result_nominal == nominal
        and all(math.isfinite(value) and 0.0 < value < 1.0 for value in tolerances)
        and result_tolerances == tolerances
        and sample_ids == list(range(len(sample_ids)))
        and len(seeds) == len(sample_ids)
        and all(seed >= 0 for seed in seeds)
        and len(set(seeds)) == len(seeds)
        and result_seeds == seeds
        and result_sample_ids == sample_ids
        and len(set(failed_ids)) == len(failed_ids)
        and set(failed_ids).issubset(sample_ids)
        and result_failed_ids == failed_ids
        and accepted_ids == expected_accepted
        and result_accepted_ids == accepted_ids
        and bool(str(contract.get("failure_criterion") or ""))
        and contract.get("result_failure_criterion")
        == contract.get("failure_criterion")
        and math.isfinite(yield_fraction)
        and math.isclose(
            yield_fraction, expected_yield, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and result_yield_fraction == yield_fraction
        and _is_sha256(str(contract.get("circuit_owner_sha256") or ""))
        and contract.get("result_circuit_owner_sha256")
        == contract.get("circuit_owner_sha256")
        and _is_sha256(str(contract.get("sample_table_sha256") or ""))
        and contract.get("result_sample_table_sha256")
        == contract.get("sample_table_sha256")
        and _is_sha256(str(contract.get("yield_result_sha256") or ""))
        and contract.get("accepted_yield_result_sha256")
        == contract.get("yield_result_sha256")
    )


def _is_sha256(value: str) -> bool:
    return len(value) == 64 and all(
        character in "0123456789abcdef" for character in value
    )


def _behavioral_source_event_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "behavioral_source_event_timestep_derivative_charge_energy_initial_owner_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        event_time = float(contract.get("event_time_s"))
        result_event_time = float(contract.get("result_event_time_s"))
        time_grid = [float(value) for value in contract.get("time_grid_s", [])]
        result_time_grid = [
            float(value) for value in contract.get("result_time_grid_s", [])
        ]
        charge = float(contract.get("integrated_charge_c"))
        result_charge = float(contract.get("result_integrated_charge_c"))
        energy = float(contract.get("source_energy_j"))
        result_energy = float(contract.get("result_source_energy_j"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("behavioral_generation_id") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "event_behavioral_generation_id",
                "timestep_behavioral_generation_id",
                "derivative_behavioral_generation_id",
                "charge_behavioral_generation_id",
                "energy_behavioral_generation_id",
                "initial_behavioral_generation_id",
                "owner_behavioral_generation_id",
                "result_behavioral_generation_id",
            )
        )
        and math.isfinite(event_time)
        and event_time >= 0.0
        and result_event_time == event_time
        and len(time_grid) >= 3
        and all(math.isfinite(value) and value >= 0.0 for value in time_grid)
        and all(left < right for left, right in zip(time_grid, time_grid[1:]))
        and any(
            math.isclose(value, event_time, rel_tol=0.0, abs_tol=1.0e-15)
            for value in time_grid
        )
        and result_time_grid == time_grid
        and contract.get("derivative_convention") == "right_limit_after_event"
        and contract.get("result_derivative_convention")
        == contract.get("derivative_convention")
        and isinstance(contract.get("initial_state"), Mapping)
        and contract.get("result_initial_state") == contract.get("initial_state")
        and math.isfinite(charge)
        and result_charge == charge
        and math.isfinite(energy)
        and energy >= 0.0
        and result_energy == energy
        and _is_sha256(str(contract.get("waveform_owner_sha256") or ""))
        and contract.get("result_waveform_owner_sha256")
        == contract.get("waveform_owner_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _touchstone_network_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "touchstone_impedance_frequency_parameter_port_complex_passivity_file_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        impedance = float(contract.get("reference_impedance_ohm"))
        result_impedance = float(contract.get("result_reference_impedance_ohm"))
        ports = [int(value) for value in contract.get("port_order", [])]
        result_ports = [int(value) for value in contract.get("result_port_order", [])]
        singular_value = float(contract.get("maximum_singular_value"))
        result_singular_value = float(contract.get("result_maximum_singular_value"))
        tolerance = float(contract.get("passivity_tolerance"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("touchstone_generation_id") or "")
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "impedance_touchstone_generation_id",
                "frequency_touchstone_generation_id",
                "parameter_touchstone_generation_id",
                "port_touchstone_generation_id",
                "complex_touchstone_generation_id",
                "passivity_touchstone_generation_id",
                "file_touchstone_generation_id",
                "result_touchstone_generation_id",
            )
        )
        and math.isfinite(impedance)
        and impedance > 0.0
        and result_impedance == impedance
        and contract.get("frequency_unit") in {"Hz", "kHz", "MHz", "GHz"}
        and contract.get("result_frequency_unit") == contract.get("frequency_unit")
        and contract.get("parameter_type") == "S"
        and contract.get("result_parameter_type") == contract.get("parameter_type")
        and len(ports) >= 2
        and all(port > 0 for port in ports)
        and len(set(ports)) == len(ports)
        and result_ports == ports
        and contract.get("complex_format") in {"RI", "MA", "DB"}
        and contract.get("result_complex_format") == contract.get("complex_format")
        and math.isfinite(tolerance)
        and tolerance >= 0.0
        and math.isfinite(singular_value)
        and singular_value >= 0.0
        and singular_value <= 1.0 + tolerance
        and result_singular_value == singular_value
        and _is_sha256(str(contract.get("touchstone_file_sha256") or ""))
        and contract.get("parsed_touchstone_file_sha256")
        == contract.get("touchstone_file_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _smps_startup_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "smps_startup_softstart_uvlo_switch_cycle_timestep_energy_waveform_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        softstart_time = [float(value) for value in contract.get("softstart_time_s", [])]
        result_softstart_time = [
            float(value) for value in contract.get("result_softstart_time_s", [])
        ]
        softstart_command = [
            float(value) for value in contract.get("softstart_command", [])
        ]
        result_softstart_command = [
            float(value) for value in contract.get("result_softstart_command", [])
        ]
        uvlo_on = float(contract.get("uvlo_on_v"))
        result_uvlo_on = float(contract.get("result_uvlo_on_v"))
        uvlo_off = float(contract.get("uvlo_off_v"))
        result_uvlo_off = float(contract.get("result_uvlo_off_v"))
        first_cycle = [
            float(value) for value in contract.get("first_switching_cycle_s", [])
        ]
        result_first_cycle = [
            float(value)
            for value in contract.get("result_first_switching_cycle_s", [])
        ]
        timestep_grid = [
            float(value) for value in contract.get("aligned_timestep_grid_s", [])
        ]
        result_timestep_grid = [
            float(value)
            for value in contract.get("result_aligned_timestep_grid_s", [])
        ]
        energies = [
            float(contract.get(key))
            for key in (
                "input_energy_j",
                "output_energy_j",
                "stored_energy_j",
                "loss_energy_j",
            )
        ]
        result_energies = [
            float(contract.get(key))
            for key in (
                "result_input_energy_j",
                "result_output_energy_j",
                "result_stored_energy_j",
                "result_loss_energy_j",
            )
        ]
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("startup_generation_id") or "")
    energy_scale = max(abs(energies[0]), 1.0e-30)
    energy_residual = energies[0] - sum(energies[1:])
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "softstart_startup_generation_id",
                "uvlo_startup_generation_id",
                "switch_startup_generation_id",
                "timestep_startup_generation_id",
                "energy_startup_generation_id",
                "waveform_startup_generation_id",
                "result_startup_generation_id",
            )
        )
        and len(softstart_time) >= 2
        and len(softstart_command) == len(softstart_time)
        and all(math.isfinite(value) and value >= 0.0 for value in softstart_time)
        and all(left < right for left, right in zip(softstart_time, softstart_time[1:]))
        and all(math.isfinite(value) and 0.0 <= value <= 1.0 for value in softstart_command)
        and all(left <= right for left, right in zip(softstart_command, softstart_command[1:]))
        and result_softstart_time == softstart_time
        and result_softstart_command == softstart_command
        and math.isfinite(uvlo_on)
        and math.isfinite(uvlo_off)
        and uvlo_on > uvlo_off >= 0.0
        and result_uvlo_on == uvlo_on
        and result_uvlo_off == uvlo_off
        and len(first_cycle) == 2
        and all(math.isfinite(value) for value in first_cycle)
        and softstart_time[-1] <= first_cycle[0] < first_cycle[1]
        and result_first_cycle == first_cycle
        and len(timestep_grid) >= 3
        and all(math.isfinite(value) for value in timestep_grid)
        and all(left < right for left, right in zip(timestep_grid, timestep_grid[1:]))
        and math.isclose(timestep_grid[0], first_cycle[0], rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(timestep_grid[-1], first_cycle[1], rel_tol=0.0, abs_tol=1.0e-15)
        and result_timestep_grid == timestep_grid
        and all(math.isfinite(value) and value >= 0.0 for value in energies)
        and result_energies == energies
        and abs(energy_residual) <= 1.0e-12 * energy_scale
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("result_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _correlation_matrix_is_positive_semidefinite(
    matrix: list[list[float]],
) -> bool:
    """Check a symmetric matrix with a pivoted LDL-style factorization."""
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        return False
    tolerance = 1.0e-12 * max(
        1.0, max(abs(value) for row in matrix for value in row)
    ) * size
    lower = [[0.0] * size for _ in range(size)]
    diagonal = [0.0] * size
    for row in range(size):
        lower[row][row] = 1.0
        for column in range(row):
            residual = matrix[row][column] - sum(
                lower[row][index]
                * diagonal[index]
                * lower[column][index]
                for index in range(column)
            )
            if abs(diagonal[column]) <= tolerance:
                if abs(residual) > tolerance:
                    return False
                lower[row][column] = 0.0
            else:
                lower[row][column] = residual / diagonal[column]
        pivot = matrix[row][row] - sum(
            lower[row][index] * lower[row][index] * diagonal[index]
            for index in range(row)
        )
        if pivot < -tolerance:
            return False
        diagonal[row] = 0.0 if abs(pivot) <= tolerance else pivot
    return True


def _noise_correlation_band_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "noise_source_correlation_psd_grid_bandwidth_transfer_integration_rms_model_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        source_order = [str(value) for value in contract.get("source_order", [])]
        result_source_order = [
            str(value) for value in contract.get("result_source_order", [])
        ]
        correlation = [
            [float(value) for value in row]
            for row in contract.get("source_correlation", [])
        ]
        result_correlation = [
            [float(value) for value in row]
            for row in contract.get("result_source_correlation", [])
        ]
        frequency = [float(value) for value in contract.get("frequency_hz", [])]
        result_frequency = [
            float(value) for value in contract.get("result_frequency_hz", [])
        ]
        bandwidth = [
            float(value) for value in contract.get("integration_bandwidth_hz", [])
        ]
        result_bandwidth = [
            float(value)
            for value in contract.get("result_integration_bandwidth_hz", [])
        ]
        transfer = [
            float(value) for value in contract.get("transfer_magnitude", [])
        ]
        result_transfer = [
            float(value) for value in contract.get("result_transfer_magnitude", [])
        ]
        output_psd = [
            float(value) for value in contract.get("output_psd_v2_per_hz", [])
        ]
        result_output_psd = [
            float(value)
            for value in contract.get("result_output_psd_v2_per_hz", [])
        ]
        integrated = float(contract.get("integrated_noise_v2"))
        result_integrated = float(contract.get("result_integrated_noise_v2"))
        rms_noise = float(contract.get("rms_noise_v"))
        result_rms_noise = float(contract.get("result_rms_noise_v"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("noise_generation_id") or "")
    source_count = len(source_order)
    correlation_shape_ok = (
        source_count > 0
        and len(correlation) == source_count
        and all(len(row) == source_count for row in correlation)
    )
    correlation_ok = correlation_shape_ok and all(
        math.isfinite(correlation[row][column])
        and abs(correlation[row][column]) <= 1.0
        and math.isclose(
            correlation[row][column],
            correlation[column][row],
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and (
            row != column
            or math.isclose(correlation[row][column], 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        )
        for row in range(source_count)
        for column in range(source_count)
    ) and _correlation_matrix_is_positive_semidefinite(correlation)
    trapezoid = sum(
        0.5 * (left_psd + right_psd) * (right_frequency - left_frequency)
        for left_frequency, right_frequency, left_psd, right_psd in zip(
            frequency, frequency[1:], output_psd, output_psd[1:]
        )
    )
    integration_scale = max(abs(integrated), abs(trapezoid), 1.0e-30)
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "correlation_noise_generation_id",
                "psd_noise_generation_id",
                "grid_noise_generation_id",
                "bandwidth_noise_generation_id",
                "transfer_noise_generation_id",
                "integration_noise_generation_id",
                "model_noise_generation_id",
                "result_noise_generation_id",
            )
        )
        and all(source_order)
        and len(set(source_order)) == source_count
        and result_source_order == source_order
        and correlation_ok
        and result_correlation == correlation
        and contract.get("psd_convention") == "one_sided_v2_per_hz"
        and contract.get("result_psd_convention") == contract.get("psd_convention")
        and len(frequency) >= 2
        and all(math.isfinite(value) and value > 0.0 for value in frequency)
        and all(left < right for left, right in zip(frequency, frequency[1:]))
        and result_frequency == frequency
        and len(bandwidth) == 2
        and frequency[0] <= bandwidth[0] < bandwidth[1] <= frequency[-1]
        and math.isclose(bandwidth[0], frequency[0], rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(bandwidth[1], frequency[-1], rel_tol=0.0, abs_tol=1.0e-15)
        and result_bandwidth == bandwidth
        and len(transfer) == len(frequency)
        and all(math.isfinite(value) and value >= 0.0 for value in transfer)
        and result_transfer == transfer
        and len(output_psd) == len(frequency)
        and all(math.isfinite(value) and value >= 0.0 for value in output_psd)
        and result_output_psd == output_psd
        and math.isfinite(integrated)
        and integrated >= 0.0
        and abs(integrated - trapezoid) <= 1.0e-12 * integration_scale
        and result_integrated == integrated
        and math.isfinite(rms_noise)
        and rms_noise >= 0.0
        and math.isclose(rms_noise * rms_noise, integrated, rel_tol=1.0e-12, abs_tol=1.0e-30)
        and result_rms_noise == rms_noise
        and _is_sha256(str(contract.get("noise_model_sha256") or ""))
        and contract.get("result_noise_model_sha256") == contract.get("noise_model_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
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
