"""Internal circuit evidence checks for LTspice identities."""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from ._identity_common import (
    _finite,
    _is_sha256,
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
