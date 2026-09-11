"""Internal dynamic evidence checks for LTspice identities."""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from ._identity_common import (
    _is_sha256,
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


def _transmission_line_transient_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "transmission_line_z0_delay_reflection_arrival_polarity_causality_energy_waveform_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        z0 = float(contract.get("characteristic_impedance_ohm"))
        source_z = float(contract.get("source_impedance_ohm"))
        load_z = float(contract.get("load_impedance_ohm"))
        source_gamma = float(contract.get("source_reflection_coefficient"))
        load_gamma = float(contract.get("load_reflection_coefficient"))
        delay = float(contract.get("one_way_delay_s"))
        incident_arrival = float(contract.get("incident_arrival_s"))
        reflected_arrival = float(contract.get("reflected_source_arrival_s"))
        time = [float(value) for value in contract.get("time_s", [])]
        result_time = [float(value) for value in contract.get("result_time_s", [])]
        waveform = [
            float(value) for value in contract.get("source_observation_v", [])
        ]
        result_waveform = [
            float(value)
            for value in contract.get("result_source_observation_v", [])
        ]
        pre_arrival = float(contract.get("pre_arrival_max_abs_v"))
        incident_energy = float(contract.get("incident_energy_j"))
        reflected_energy = float(contract.get("reflected_energy_j"))
        accepted_energy = float(contract.get("accepted_energy_j"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("transmission_line_generation_id") or "")
    values = (z0, source_z, load_z, source_gamma, load_gamma, delay)
    expected_source_gamma = (source_z - z0) / (source_z + z0)
    expected_load_gamma = (load_z - z0) / (load_z + z0)
    expected_reflected_polarity = "positive" if load_gamma > 0.0 else "negative"
    causal_samples = [
        abs(value) for sample_time, value in zip(time, waveform) if sample_time < delay
    ]
    causal_max = max(causal_samples, default=0.0)
    energy_scale = max(abs(incident_energy), 1.0e-30)
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "impedance_transmission_line_generation_id",
                "delay_transmission_line_generation_id",
                "reflection_transmission_line_generation_id",
                "arrival_transmission_line_generation_id",
                "polarity_transmission_line_generation_id",
                "causality_transmission_line_generation_id",
                "energy_transmission_line_generation_id",
                "waveform_transmission_line_generation_id",
                "result_transmission_line_generation_id",
            )
        )
        and all(math.isfinite(value) for value in values)
        and z0 > 0.0
        and source_z > 0.0
        and load_z > 0.0
        and delay > 0.0
        and float(contract.get("result_characteristic_impedance_ohm")) == z0
        and float(contract.get("result_source_impedance_ohm")) == source_z
        and float(contract.get("result_load_impedance_ohm")) == load_z
        and math.isclose(source_gamma, expected_source_gamma, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(load_gamma, expected_load_gamma, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and float(contract.get("result_source_reflection_coefficient")) == source_gamma
        and float(contract.get("result_load_reflection_coefficient")) == load_gamma
        and math.isclose(incident_arrival, delay, rel_tol=0.0, abs_tol=1.0e-18)
        and math.isclose(reflected_arrival, 2.0 * delay, rel_tol=0.0, abs_tol=1.0e-18)
        and float(contract.get("result_one_way_delay_s")) == delay
        and float(contract.get("result_incident_arrival_s")) == incident_arrival
        and float(contract.get("result_reflected_source_arrival_s")) == reflected_arrival
        and contract.get("incident_pulse_polarity") == "positive"
        and contract.get("result_incident_pulse_polarity")
        == contract.get("incident_pulse_polarity")
        and contract.get("reflected_pulse_polarity") == expected_reflected_polarity
        and contract.get("result_reflected_pulse_polarity")
        == contract.get("reflected_pulse_polarity")
        and len(time) >= 4
        and len(waveform) == len(time)
        and all(math.isfinite(value) for value in (*time, *waveform))
        and all(left < right for left, right in zip(time, time[1:]))
        and result_time == time
        and result_waveform == waveform
        and math.isfinite(pre_arrival)
        and math.isclose(pre_arrival, causal_max, rel_tol=0.0, abs_tol=1.0e-15)
        and pre_arrival <= 1.0e-12
        and float(contract.get("result_pre_arrival_max_abs_v")) == pre_arrival
        and all(
            math.isfinite(value) and value >= 0.0
            for value in (incident_energy, reflected_energy, accepted_energy)
        )
        and abs(incident_energy - reflected_energy - accepted_energy)
        <= 1.0e-12 * energy_scale
        and math.isclose(
            reflected_energy / incident_energy,
            load_gamma * load_gamma,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and float(contract.get("result_incident_energy_j")) == incident_energy
        and float(contract.get("result_reflected_energy_j")) == reflected_energy
        and float(contract.get("result_accepted_energy_j")) == accepted_energy
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("result_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _sampled_loop_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "sampled_loop_sideband_injection_period_fft_bin_phase_nyquist_crossover_waveform_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        sideband = int(contract.get("sideband_order"))
        period = float(contract.get("switching_period_s"))
        switching_frequency = float(contract.get("switching_frequency_hz"))
        duration = float(contract.get("record_duration_s"))
        injection_frequency = float(contract.get("injection_frequency_hz"))
        coherent_bin = int(contract.get("coherent_fft_bin"))
        frequency = [float(value) for value in contract.get("loop_frequency_hz", [])]
        result_frequency = [
            float(value) for value in contract.get("result_loop_frequency_hz", [])
        ]
        magnitude = [float(value) for value in contract.get("loop_magnitude", [])]
        result_magnitude = [
            float(value) for value in contract.get("result_loop_magnitude", [])
        ]
        phase = [float(value) for value in contract.get("loop_phase_deg", [])]
        result_phase = [
            float(value) for value in contract.get("result_loop_phase_deg", [])
        ]
        crossover = float(contract.get("crossover_frequency_hz"))
        encirclements = int(
            contract.get("nyquist_clockwise_encirclements_minus_one")
        )
        open_rhp = int(contract.get("open_loop_right_half_plane_poles"))
        closed_rhp = int(contract.get("closed_loop_right_half_plane_poles"))
    except (TypeError, ValueError):
        return False
    generation = str(contract.get("sampled_loop_generation_id") or "")
    unity_indices = [
        index
        for index, value in enumerate(magnitude)
        if math.isclose(value, 1.0, rel_tol=1.0e-12, abs_tol=1.0e-15)
    ]
    coherent_bin_exact = injection_frequency * duration
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "sideband_sampled_loop_generation_id",
                "injection_sampled_loop_generation_id",
                "period_sampled_loop_generation_id",
                "fft_sampled_loop_generation_id",
                "phase_sampled_loop_generation_id",
                "nyquist_sampled_loop_generation_id",
                "crossover_sampled_loop_generation_id",
                "waveform_sampled_loop_generation_id",
                "result_sampled_loop_generation_id",
            )
        )
        and sideband == 0
        and contract.get("result_sideband_order") == sideband
        and contract.get("sideband_selection") == "fundamental"
        and contract.get("result_sideband_selection") == contract.get("sideband_selection")
        and contract.get("injection_point") == "control_to_duty_break"
        and contract.get("result_injection_point") == contract.get("injection_point")
        and math.isfinite(period)
        and period > 0.0
        and math.isclose(switching_frequency, 1.0 / period, rel_tol=1.0e-12)
        and float(contract.get("result_switching_period_s")) == period
        and float(contract.get("result_switching_frequency_hz")) == switching_frequency
        and math.isfinite(duration)
        and duration > period
        and float(contract.get("result_record_duration_s")) == duration
        and math.isfinite(injection_frequency)
        and 0.0 < injection_frequency < 0.5 * switching_frequency
        and float(contract.get("result_injection_frequency_hz")) == injection_frequency
        and coherent_bin > 0
        and math.isclose(coherent_bin_exact, coherent_bin, rel_tol=0.0, abs_tol=1.0e-12)
        and contract.get("result_coherent_fft_bin") == coherent_bin
        and len(frequency) >= 3
        and len(magnitude) == len(frequency) == len(phase)
        and all(math.isfinite(value) and value > 0.0 for value in frequency)
        and all(left < right for left, right in zip(frequency, frequency[1:]))
        and all(math.isfinite(value) and value > 0.0 for value in magnitude)
        and all(math.isfinite(value) and -360.0 <= value <= 360.0 for value in phase)
        and result_frequency == frequency
        and result_magnitude == magnitude
        and result_phase == phase
        and len(unity_indices) == 1
        and math.isclose(crossover, frequency[unity_indices[0]], rel_tol=1.0e-12)
        and float(contract.get("result_crossover_frequency_hz")) == crossover
        and open_rhp >= 0
        and closed_rhp >= 0
        and closed_rhp == open_rhp + encirclements
        and contract.get("result_nyquist_clockwise_encirclements_minus_one")
        == encirclements
        and contract.get("result_open_loop_right_half_plane_poles") == open_rhp
        and contract.get("result_closed_loop_right_half_plane_poles") == closed_rhp
        and closed_rhp == 0
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("result_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _periodic_steady_state_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "periodic_steady_state_charge_flux_cycle_energy_efficiency_phase_waveform_owner_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("periodic_generation_id") or "")
    try:
        start_charge = [float(item) for item in contract["cycle_start_charge_c"]]
        end_charge = [float(item) for item in contract["cycle_end_charge_c"]]
        result_start_charge = [float(item) for item in contract["result_cycle_start_charge_c"]]
        result_end_charge = [float(item) for item in contract["result_cycle_end_charge_c"]]
        start_flux = [float(item) for item in contract["cycle_start_flux_wb"]]
        end_flux = [float(item) for item in contract["cycle_end_flux_wb"]]
        result_start_flux = [float(item) for item in contract["result_cycle_start_flux_wb"]]
        result_end_flux = [float(item) for item in contract["result_cycle_end_flux_wb"]]
        tolerance = float(contract["closure_tolerance"])
        result_tolerance = float(contract["result_closure_tolerance"])
        input_energy = float(contract["input_energy_j"])
        output_energy = float(contract["output_energy_j"])
        loss_energy = float(contract["loss_energy_j"])
        efficiency = float(contract["efficiency"])
        phase_window = [float(item) for item in contract["phase_window_rad"]]
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "charge_periodic_generation_id",
                "flux_periodic_generation_id",
                "energy_periodic_generation_id",
                "efficiency_periodic_generation_id",
                "phase_periodic_generation_id",
                "waveform_periodic_generation_id",
                "owner_periodic_generation_id",
                "result_periodic_generation_id",
            )
        )
        and bool(start_charge)
        and len(end_charge) == len(start_charge)
        and all(math.isfinite(item) for item in (*start_charge, *end_charge))
        and math.isfinite(tolerance)
        and 0.0 < tolerance <= 1.0e-6
        and result_tolerance == tolerance
        and max(abs(left - right) for left, right in zip(start_charge, end_charge))
        <= tolerance
        and result_start_charge == start_charge
        and result_end_charge == end_charge
        and bool(start_flux)
        and len(end_flux) == len(start_flux)
        and all(math.isfinite(item) for item in (*start_flux, *end_flux))
        and max(abs(left - right) for left, right in zip(start_flux, end_flux))
        <= tolerance
        and result_start_flux == start_flux
        and result_end_flux == end_flux
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (input_energy, output_energy, loss_energy)
        )
        and input_energy > 0.0
        and math.isclose(
            input_energy,
            output_energy + loss_energy,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isfinite(efficiency)
        and 0.0 <= efficiency <= 1.0
        and math.isclose(efficiency, output_energy / input_energy, rel_tol=1.0e-12)
        and contract.get("result_input_energy_j") == input_energy
        and contract.get("result_output_energy_j") == output_energy
        and contract.get("result_loss_energy_j") == loss_energy
        and contract.get("result_efficiency") == efficiency
        and len(phase_window) == 2
        and math.isclose(phase_window[0], 0.0, abs_tol=1.0e-12)
        and math.isclose(phase_window[1], 2.0 * math.pi, rel_tol=1.0e-12)
        and contract.get("result_phase_window_rad") == phase_window
        and bool(str(contract.get("waveform_owner") or ""))
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _small_signal_consistency_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "bias_small_signal_jacobian_ac_transient_pole_zero_normalization_circuit_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("small_signal_generation_id") or "")
    try:
        bias = [float(item) for item in contract["bias_state"]]
        result_bias = [float(item) for item in contract["result_bias_state"]]
        jacobian = [[float(item) for item in row] for row in contract["jacobian"]]
        result_jacobian = [[float(item) for item in row] for row in contract["result_jacobian"]]
        frequencies = [float(item) for item in contract["frequency_rad_s"]]
        result_frequencies = [float(item) for item in contract["result_frequency_rad_s"]]
        transfer = [complex(float(row[0]), float(row[1])) for row in contract["ac_transfer_ri"]]
        result_transfer = [complex(float(row[0]), float(row[1])) for row in contract["result_ac_transfer_ri"]]
        times = [float(item) for item in contract["time_s"]]
        result_times = [float(item) for item in contract["result_time_s"]]
        impulse = [float(item) for item in contract["impulse_response"]]
        result_impulse = [float(item) for item in contract["result_impulse_response"]]
        step = [float(item) for item in contract["step_response"]]
        result_step = [float(item) for item in contract["result_step_response"]]
        poles = [complex(float(row[0]), float(row[1])) for row in contract["poles_ri"]]
        result_poles = [complex(float(row[0]), float(row[1])) for row in contract["result_poles_ri"]]
        zeros = [complex(float(row[0]), float(row[1])) for row in contract["zeros_ri"]]
        result_zeros = [complex(float(row[0]), float(row[1])) for row in contract["result_zeros_ri"]]
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    first_order_pole = poles[0].real if len(poles) == 1 else math.nan
    physical_response = (
        len(times) == len(impulse) == len(step) >= 2
        and math.isfinite(first_order_pole)
        and first_order_pole < 0.0
        and abs(poles[0].imag) <= 1.0e-12
        and all(
            math.isclose(
                impulse[index],
                math.exp(first_order_pole * time),
                rel_tol=1.0e-10,
                abs_tol=1.0e-12,
            )
            and math.isclose(
                step[index],
                1.0 - math.exp(first_order_pole * time),
                rel_tol=1.0e-10,
                abs_tol=1.0e-12,
            )
            for index, time in enumerate(times)
        )
    )
    expected_transfer = [1.0 / complex(-first_order_pole, omega) for omega in frequencies]
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "bias_small_signal_generation_id",
                "jacobian_small_signal_generation_id",
                "ac_small_signal_generation_id",
                "transient_small_signal_generation_id",
                "pole_zero_small_signal_generation_id",
                "normalization_small_signal_generation_id",
                "circuit_small_signal_generation_id",
                "result_small_signal_generation_id",
            )
        )
        and len(bias) == 1
        and all(math.isfinite(item) for item in bias)
        and result_bias == bias
        and jacobian == [[-1.0]]
        and result_jacobian == jacobian
        and len(frequencies) == len(transfer) >= 2
        and all(math.isfinite(item) and item >= 0.0 for item in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and all(abs(actual - expected) <= 1.0e-10 for actual, expected in zip(transfer, expected_transfer))
        and result_transfer == transfer
        and all(left < right for left, right in zip(times, times[1:]))
        and result_times == times
        and physical_response
        and result_impulse == impulse
        and result_step == step
        and result_poles == poles
        and not zeros
        and result_zeros == zeros
        and contract.get("normalization") == "monic_denominator_unit_dc_gain"
        and contract.get("result_normalization") == contract.get("normalization")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("linearized_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _switched_limit_cycle_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "switched_converter_poincare_limitcycle_floquet_event_flux_charge_cycle_energy_step_waveform_owner_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("limit_cycle_generation_id") or "")
    try:
        start_state = [float(item) for item in contract["poincare_state_start"]]
        end_state = [float(item) for item in contract["poincare_state_end"]]
        result_start_state = [float(item) for item in contract["result_poincare_state_start"]]
        result_end_state = [float(item) for item in contract["result_poincare_state_end"]]
        tolerance = float(contract["poincare_tolerance"])
        result_tolerance = float(contract["result_poincare_tolerance"])
        events = [str(item) for item in contract["event_sequence"]]
        result_events = [str(item) for item in contract["result_event_sequence"]]
        event_times = [float(item) for item in contract["event_times_s"]]
        result_event_times = [float(item) for item in contract["result_event_times_s"]]
        floquet = [
            complex(float(row[0]), float(row[1])) for row in contract["floquet_multipliers_ri"]
        ]
        result_floquet = [
            complex(float(row[0]), float(row[1]))
            for row in contract["result_floquet_multipliers_ri"]
        ]
        start_flux = [float(item) for item in contract["inductor_flux_start_wb"]]
        end_flux = [float(item) for item in contract["inductor_flux_end_wb"]]
        result_start_flux = [float(item) for item in contract["result_inductor_flux_start_wb"]]
        result_end_flux = [float(item) for item in contract["result_inductor_flux_end_wb"]]
        start_charge = [float(item) for item in contract["capacitor_charge_start_c"]]
        end_charge = [float(item) for item in contract["capacitor_charge_end_c"]]
        result_start_charge = [float(item) for item in contract["result_capacitor_charge_start_c"]]
        result_end_charge = [float(item) for item in contract["result_capacitor_charge_end_c"]]
        input_energy = float(contract["input_cycle_energy_j"])
        output_energy = float(contract["output_cycle_energy_j"])
        loss_energy = float(contract["loss_cycle_energy_j"])
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    finite_rows = all(
        math.isfinite(item)
        for row in (start_state, end_state, start_flux, end_flux, start_charge, end_charge)
        for item in row
    )
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "poincare_limit_cycle_generation_id",
                "floquet_limit_cycle_generation_id",
                "event_limit_cycle_generation_id",
                "flux_limit_cycle_generation_id",
                "charge_limit_cycle_generation_id",
                "energy_limit_cycle_generation_id",
                "step_limit_cycle_generation_id",
                "waveform_limit_cycle_generation_id",
                "owner_limit_cycle_generation_id",
                "result_limit_cycle_generation_id",
            )
        )
        and bool(start_state)
        and len(end_state) == len(start_state)
        and finite_rows
        and math.isfinite(tolerance)
        and 0.0 < tolerance <= 1.0e-6
        and result_tolerance == tolerance
        and max(abs(left - right) for left, right in zip(start_state, end_state)) <= tolerance
        and result_start_state == start_state
        and result_end_state == end_state
        and len(events) == len(event_times) >= 3
        and events[0] == events[-1] == "switch_on"
        and all(left != right for left, right in zip(events, events[1:]))
        and result_events == events
        and all(math.isfinite(item) for item in event_times)
        and math.isclose(event_times[0], 0.0, abs_tol=1.0e-15)
        and all(left < right for left, right in zip(event_times, event_times[1:]))
        and result_event_times == event_times
        and bool(floquet)
        and all(
            math.isfinite(item.real) and math.isfinite(item.imag) and abs(item) < 1.0
            for item in floquet
        )
        and result_floquet == floquet
        and bool(start_flux)
        and len(end_flux) == len(start_flux)
        and max(abs(left - right) for left, right in zip(start_flux, end_flux)) <= tolerance
        and result_start_flux == start_flux
        and result_end_flux == end_flux
        and bool(start_charge)
        and len(end_charge) == len(start_charge)
        and max(abs(left - right) for left, right in zip(start_charge, end_charge)) <= tolerance
        and result_start_charge == start_charge
        and result_end_charge == end_charge
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (input_energy, output_energy, loss_energy)
        )
        and input_energy > 0.0
        and math.isclose(
            input_energy, output_energy + loss_energy, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and contract.get("result_input_cycle_energy_j") == input_energy
        and contract.get("result_output_cycle_energy_j") == output_energy
        and contract.get("result_loss_cycle_energy_j") == loss_energy
        and bool(str(contract.get("step_owner") or ""))
        and contract.get("accepted_step_owner") == contract.get("step_owner")
        and _is_sha256(str(contract.get("step_sha256") or ""))
        and contract.get("accepted_step_sha256") == contract.get("step_sha256")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("accepted_waveform_sha256") == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _noise_referral_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "noise_input_output_referred_density_correlation_integration_gain_circuit_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("noise_generation_id") or "")
    try:
        frequency = [float(item) for item in contract["frequency_hz"]]
        result_frequency = [float(item) for item in contract["result_frequency_hz"]]
        input_density = [float(item) for item in contract["input_referred_density_v_per_sqrt_hz"]]
        result_input_density = [
            float(item) for item in contract["result_input_referred_density_v_per_sqrt_hz"]
        ]
        output_density = [float(item) for item in contract["output_referred_density_v_per_sqrt_hz"]]
        result_output_density = [
            float(item) for item in contract["result_output_referred_density_v_per_sqrt_hz"]
        ]
        gain = [float(item) for item in contract["transfer_gain_magnitude"]]
        result_gain = [float(item) for item in contract["result_transfer_gain_magnitude"]]
        correlation = [
            [float(item) for item in row] for row in contract["source_correlation_matrix"]
        ]
        result_correlation = [
            [float(item) for item in row] for row in contract["result_source_correlation_matrix"]
        ]
        bandwidth = [float(item) for item in contract["integration_bandwidth_hz"]]
        result_bandwidth = [float(item) for item in contract["result_integration_bandwidth_hz"]]
        bin_width = [float(item) for item in contract["integration_bin_width_hz"]]
        result_bin_width = [float(item) for item in contract["result_integration_bin_width_hz"]]
        rms_noise = float(contract["total_output_rms_noise_v"])
        result_rms_noise = float(contract["result_total_output_rms_noise_v"])
    except (KeyError, TypeError, ValueError):
        return False
    correlation_ok = (
        bool(correlation)
        and all(len(row) == len(correlation) for row in correlation)
        and all(math.isfinite(item) for row in correlation for item in row)
        and all(
            math.isclose(correlation[i][i], 1.0, abs_tol=1.0e-12) for i in range(len(correlation))
        )
        and all(
            math.isclose(correlation[i][j], correlation[j][i], abs_tol=1.0e-12)
            and abs(correlation[i][j]) <= 1.0
            for i in range(len(correlation))
            for j in range(len(correlation))
        )
    )
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "input_noise_generation_id",
                "output_noise_generation_id",
                "correlation_noise_generation_id",
                "integration_noise_generation_id",
                "gain_noise_generation_id",
                "circuit_noise_generation_id",
                "result_noise_generation_id",
            )
        )
        and len(frequency)
        == len(input_density)
        == len(output_density)
        == len(gain)
        == len(bin_width)
        >= 2
        and all(math.isfinite(item) and item > 0.0 for item in frequency)
        and all(left < right for left, right in zip(frequency, frequency[1:]))
        and result_frequency == frequency
        and all(
            math.isfinite(item) and item >= 0.0 for item in (*input_density, *output_density, *gain)
        )
        and result_input_density == input_density
        and result_output_density == output_density
        and result_gain == gain
        and all(
            math.isclose(output, source * scale, rel_tol=1.0e-12, abs_tol=1.0e-18)
            for source, output, scale in zip(input_density, output_density, gain)
        )
        and correlation_ok
        and result_correlation == correlation
        and contract.get("spectral_density_unit") == "V/sqrt(Hz)"
        and contract.get("result_spectral_density_unit") == contract.get("spectral_density_unit")
        and bandwidth == [frequency[0], frequency[-1]]
        and result_bandwidth == bandwidth
        and all(math.isfinite(item) and item > 0.0 for item in bin_width)
        and result_bin_width == bin_width
        and math.isfinite(rms_noise)
        and rms_noise >= 0.0
        and math.isclose(
            rms_noise * rms_noise,
            sum(value * value * width for value, width in zip(output_density, bin_width)),
            rel_tol=1.0e-12,
            abs_tol=1.0e-24,
        )
        and result_rms_noise == rms_noise
        and bool(str(contract.get("circuit_owner") or ""))
        and contract.get("accepted_circuit_owner") == contract.get("circuit_owner")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _pll_closure_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "pll_lock_phase_noise_jitter_loop_gain_waveform_circuit_owner_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("pll_generation_id") or "")
    try:
        lock_time = float(contract["lock_time_s"])
        result_lock_time = float(contract["result_lock_time_s"])
        duration = float(contract["waveform_duration_s"])
        result_duration = float(contract["result_waveform_duration_s"])
        phase_start = float(contract["phase_error_start_rad"])
        result_phase_start = float(contract["result_phase_error_start_rad"])
        phase_end = float(contract["phase_error_end_rad"])
        result_phase_end = float(contract["result_phase_error_end_rad"])
        tolerance = float(contract["phase_error_tolerance_rad"])
        result_tolerance = float(contract["result_phase_error_tolerance_rad"])
        offsets = [float(item) for item in contract["phase_noise_offset_hz"]]
        result_offsets = [
            float(item) for item in contract["result_phase_noise_offset_hz"]
        ]
        phase_noise = [float(item) for item in contract["phase_noise_dbc_per_hz"]]
        result_phase_noise = [
            float(item) for item in contract["result_phase_noise_dbc_per_hz"]
        ]
        jitter = float(contract["integrated_jitter_s"])
        result_jitter = float(contract["result_integrated_jitter_s"])
        crossover = float(contract["loop_gain_crossover_hz"])
        result_crossover = float(contract["result_loop_gain_crossover_hz"])
        crossover_gain = float(contract["loop_gain_magnitude_at_crossover"])
        result_crossover_gain = float(
            contract["result_loop_gain_magnitude_at_crossover"]
        )
        phase_margin = float(contract["phase_margin_deg"])
        result_phase_margin = float(contract["result_phase_margin_deg"])
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "lock_pll_generation_id",
                "phase_error_pll_generation_id",
                "noise_pll_generation_id",
                "jitter_pll_generation_id",
                "loop_gain_pll_generation_id",
                "waveform_pll_generation_id",
                "circuit_pll_generation_id",
                "result_pll_generation_id",
            )
        )
        and contract.get("lock_state") == "locked"
        and contract.get("result_lock_state") == contract.get("lock_state")
        and all(math.isfinite(item) for item in (lock_time, duration))
        and 0.0 < lock_time <= duration
        and result_lock_time == lock_time
        and result_duration == duration
        and all(
            math.isfinite(item)
            for item in (phase_start, phase_end, tolerance)
        )
        and 0.0 < tolerance <= 1.0e-2
        and abs(phase_end) <= tolerance
        and abs(phase_end) < abs(phase_start)
        and result_phase_start == phase_start
        and result_phase_end == phase_end
        and result_tolerance == tolerance
        and len(offsets) == len(phase_noise) >= 2
        and all(math.isfinite(item) and item > 0.0 for item in offsets)
        and all(left < right for left, right in zip(offsets, offsets[1:]))
        and all(math.isfinite(item) and item <= 0.0 for item in phase_noise)
        and all(right <= left for left, right in zip(phase_noise, phase_noise[1:]))
        and result_offsets == offsets
        and result_phase_noise == phase_noise
        and math.isfinite(jitter)
        and jitter > 0.0
        and result_jitter == jitter
        and math.isfinite(crossover)
        and crossover > 0.0
        and result_crossover == crossover
        and math.isfinite(crossover_gain)
        and math.isclose(crossover_gain, 1.0, rel_tol=5.0e-2, abs_tol=5.0e-2)
        and result_crossover_gain == crossover_gain
        and math.isfinite(phase_margin)
        and 0.0 < phase_margin < 180.0
        and result_phase_margin == phase_margin
        and bool(str(contract.get("waveform_owner") or ""))
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("accepted_waveform_sha256")
        == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _electrothermal_fixedpoint_owner_identity_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get(
        "electrothermal_thermal_impedance_power_temperature_device_fixedpoint_stability_circuit_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("electrothermal_generation_id") or "")
    try:
        thermal_impedance = float(contract["thermal_impedance_k_per_w"])
        result_thermal_impedance = float(
            contract["result_thermal_impedance_k_per_w"]
        )
        power = float(contract["dissipated_power_w"])
        result_power = float(contract["result_dissipated_power_w"])
        ambient = float(contract["ambient_temperature_c"])
        result_ambient = float(contract["result_ambient_temperature_c"])
        junction = float(contract["junction_temperature_c"])
        result_junction = float(contract["result_junction_temperature_c"])
        reference_power = float(contract["reference_power_w"])
        result_reference_power = float(contract["result_reference_power_w"])
        reference_temperature = float(contract["reference_temperature_c"])
        result_reference_temperature = float(
            contract["result_reference_temperature_c"]
        )
        coefficient = float(contract["power_temperature_coefficient_per_k"])
        result_coefficient = float(
            contract["result_power_temperature_coefficient_per_k"]
        )
        residual = float(contract["fixedpoint_residual_c"])
        result_residual = float(contract["result_fixedpoint_residual_c"])
        residual_tolerance = float(contract["fixedpoint_tolerance_c"])
        result_residual_tolerance = float(
            contract["result_fixedpoint_tolerance_c"]
        )
        stability_slope = float(contract["stability_slope"])
        result_stability_slope = float(contract["result_stability_slope"])
    except (KeyError, TypeError, ValueError):
        return False
    device_power = reference_power * (
        1.0 + coefficient * (junction - reference_temperature)
    )
    computed_residual = junction - (ambient + thermal_impedance * device_power)
    computed_slope = thermal_impedance * reference_power * coefficient
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "thermal_impedance_generation_id",
                "power_generation_id",
                "temperature_generation_id",
                "device_law_generation_id",
                "fixedpoint_generation_id",
                "stability_generation_id",
                "circuit_generation_id",
                "result_generation_id",
            )
        )
        and all(
            math.isfinite(item)
            for item in (
                thermal_impedance,
                power,
                ambient,
                junction,
                reference_power,
                reference_temperature,
                coefficient,
                residual,
                residual_tolerance,
                stability_slope,
            )
        )
        and thermal_impedance > 0.0
        and power > 0.0
        and reference_power > 0.0
        and result_thermal_impedance == thermal_impedance
        and result_power == power
        and result_ambient == ambient
        and result_junction == junction
        and math.isclose(
            junction,
            ambient + thermal_impedance * power,
            rel_tol=1.0e-12,
            abs_tol=1.0e-10,
        )
        and contract.get("device_law") == "linear_temperature_power"
        and contract.get("result_device_law") == contract.get("device_law")
        and result_reference_power == reference_power
        and result_reference_temperature == reference_temperature
        and result_coefficient == coefficient
        and math.isclose(device_power, power, rel_tol=1.0e-12, abs_tol=1.0e-10)
        and residual_tolerance > 0.0
        and abs(residual) <= residual_tolerance
        and math.isclose(residual, computed_residual, rel_tol=1.0e-10, abs_tol=1.0e-10)
        and result_residual == residual
        and result_residual_tolerance == residual_tolerance
        and math.isclose(
            stability_slope, computed_slope, rel_tol=1.0e-12, abs_tol=1.0e-12
        )
        and abs(stability_slope) < 1.0
        and result_stability_slope == stability_slope
        and bool(str(contract.get("circuit_owner") or ""))
        and contract.get("accepted_circuit_owner") == contract.get("circuit_owner")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _oscillator_closure_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "oscillator_startup_barkhausen_frequency_amplitude_limitcycle_energy_owner_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("oscillator_generation_id") or "")
    try:
        startup_time = [float(item) for item in contract["startup_time_s"]]
        result_startup_time = [
            float(item) for item in contract["result_startup_time_s"]
        ]
        startup_amplitude = [
            float(item) for item in contract["startup_amplitude_v"]
        ]
        result_startup_amplitude = [
            float(item) for item in contract["result_startup_amplitude_v"]
        ]
        barkhausen_frequency = float(contract["barkhausen_frequency_hz"])
        result_barkhausen_frequency = float(
            contract["result_barkhausen_frequency_hz"]
        )
        loop_gain = float(contract["loop_gain_magnitude"])
        result_loop_gain = float(contract["result_loop_gain_magnitude"])
        loop_phase = float(contract["loop_phase_deg"])
        result_loop_phase = float(contract["result_loop_phase_deg"])
        steady_frequency = float(contract["steady_frequency_hz"])
        result_steady_frequency = float(contract["result_steady_frequency_hz"])
        steady_amplitude = float(contract["steady_amplitude_v"])
        result_steady_amplitude = float(contract["result_steady_amplitude_v"])
        start_state = [float(item) for item in contract["limitcycle_start_state"]]
        end_state = [float(item) for item in contract["limitcycle_end_state"]]
        result_start_state = [
            float(item) for item in contract["result_limitcycle_start_state"]
        ]
        result_end_state = [
            float(item) for item in contract["result_limitcycle_end_state"]
        ]
        state_tolerance = float(contract["limitcycle_state_tolerance"])
        result_state_tolerance = float(
            contract["result_limitcycle_state_tolerance"]
        )
        source_energy = float(contract["cycle_source_energy_j"])
        result_source_energy = float(contract["result_cycle_source_energy_j"])
        dissipated_energy = float(contract["cycle_dissipated_energy_j"])
        result_dissipated_energy = float(
            contract["result_cycle_dissipated_energy_j"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    try:
        growth_rates = [
            math.log(right_amp / left_amp) / (right_time - left_time)
            for left_time, right_time, left_amp, right_amp in zip(
                startup_time,
                startup_time[1:],
                startup_amplitude,
                startup_amplitude[1:],
            )
        ]
    except (ValueError, ZeroDivisionError):
        return False
    phase_residual = (loop_phase + 180.0) % 360.0 - 180.0
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "startup_oscillator_generation_id",
                "barkhausen_oscillator_generation_id",
                "frequency_oscillator_generation_id",
                "amplitude_oscillator_generation_id",
                "limitcycle_oscillator_generation_id",
                "energy_oscillator_generation_id",
                "waveform_oscillator_generation_id",
                "circuit_oscillator_generation_id",
                "result_oscillator_generation_id",
            )
        )
        and len(startup_time) == len(startup_amplitude) >= 4
        and all(math.isfinite(item) for item in startup_time + startup_amplitude)
        and startup_time[0] >= 0.0
        and all(left < right for left, right in zip(startup_time, startup_time[1:]))
        and all(item > 0.0 for item in startup_amplitude)
        and all(
            left < right for left, right in zip(startup_amplitude, startup_amplitude[1:])
        )
        and result_startup_time == startup_time
        and result_startup_amplitude == startup_amplitude
        and all(math.isfinite(item) and item > 0.0 for item in growth_rates)
        and max(growth_rates) - min(growth_rates)
        <= 1.0e-10 * max(growth_rates)
        and all(
            math.isfinite(item)
            for item in (
                barkhausen_frequency,
                loop_gain,
                loop_phase,
                steady_frequency,
                steady_amplitude,
            )
        )
        and barkhausen_frequency > 0.0
        and steady_frequency > 0.0
        and math.isclose(loop_gain, 1.0, rel_tol=5.0e-2, abs_tol=5.0e-2)
        and abs(phase_residual) <= 5.0
        and math.isclose(
            steady_frequency, barkhausen_frequency, rel_tol=2.0e-2, abs_tol=0.0
        )
        and steady_amplitude > 0.0
        and result_barkhausen_frequency == barkhausen_frequency
        and result_loop_gain == loop_gain
        and result_loop_phase == loop_phase
        and result_steady_frequency == steady_frequency
        and result_steady_amplitude == steady_amplitude
        and len(start_state) == len(end_state) >= 2
        and result_start_state == start_state
        and result_end_state == end_state
        and math.isfinite(state_tolerance)
        and state_tolerance > 0.0
        and result_state_tolerance == state_tolerance
        and max(abs(left - right) for left, right in zip(start_state, end_state))
        <= state_tolerance
        and source_energy > 0.0
        and dissipated_energy > 0.0
        and math.isclose(
            source_energy, dissipated_energy, rel_tol=1.0e-9, abs_tol=1.0e-18
        )
        and result_source_energy == source_energy
        and result_dissipated_energy == dissipated_energy
        and bool(str(contract.get("waveform_owner") or ""))
        and contract.get("accepted_waveform_owner") == contract.get("waveform_owner")
        and _is_sha256(str(contract.get("waveform_sha256") or ""))
        and contract.get("accepted_waveform_sha256")
        == contract.get("waveform_sha256")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _conducted_emi_owner_identity_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get(
        "conducted_emi_lisn_spectrum_window_band_limit_power_owner_result_identity"
    )
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("emi_generation_id") or "")
    try:
        lisn_impedance = float(contract["lisn_impedance_ohm"])
        result_lisn_impedance = float(contract["result_lisn_impedance_ohm"])
        bandwidth = float(contract["resolution_bandwidth_hz"])
        result_bandwidth = float(contract["result_resolution_bandwidth_hz"])
        frequency = [float(item) for item in contract["frequency_hz"]]
        result_frequency = [float(item) for item in contract["result_frequency_hz"]]
        spectrum = [float(item) for item in contract["spectrum_dbuv"]]
        result_spectrum = [float(item) for item in contract["result_spectrum_dbuv"]]
        limits = [float(item) for item in contract["limit_dbuv"]]
        result_limits = [float(item) for item in contract["result_limit_dbuv"]]
        margins = [float(item) for item in contract["margin_db"]]
        result_margins = [float(item) for item in contract["result_margin_db"]]
        source_power = float(contract["source_input_power_w"])
        result_source_power = float(contract["result_source_input_power_w"])
        load_power = float(contract["load_power_w"])
        result_load_power = float(contract["result_load_power_w"])
        lisn_loss = float(contract["lisn_loss_w"])
        result_lisn_loss = float(contract["result_lisn_loss_w"])
        switching_loss = float(contract["switching_loss_w"])
        result_switching_loss = float(contract["result_switching_loss_w"])
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            contract.get(key) == generation
            for key in (
                "lisn_emi_generation_id",
                "waveform_emi_generation_id",
                "fft_emi_generation_id",
                "detector_emi_generation_id",
                "limit_emi_generation_id",
                "power_emi_generation_id",
                "circuit_emi_generation_id",
                "result_emi_generation_id",
            )
        )
        and math.isfinite(lisn_impedance)
        and math.isclose(lisn_impedance, 50.0, rel_tol=1.0e-2, abs_tol=0.0)
        and result_lisn_impedance == lisn_impedance
        and contract.get("waveform_window") in {"hann", "blackman_harris"}
        and contract.get("result_waveform_window") == contract.get("waveform_window")
        and contract.get("fft_normalization") == "rms_single_sided"
        and contract.get("result_fft_normalization")
        == contract.get("fft_normalization")
        and contract.get("detector") == "quasi_peak"
        and contract.get("result_detector") == contract.get("detector")
        and math.isfinite(bandwidth)
        and bandwidth > 0.0
        and result_bandwidth == bandwidth
        and len(frequency) == len(spectrum) == len(limits) == len(margins) >= 3
        and all(math.isfinite(item) and item > 0.0 for item in frequency)
        and all(left < right for left, right in zip(frequency, frequency[1:]))
        and frequency[0] >= 1.5e5
        and frequency[-1] <= 3.0e7
        and result_frequency == frequency
        and all(math.isfinite(item) for item in spectrum + limits + margins)
        and result_spectrum == spectrum
        and result_limits == limits
        and result_margins == margins
        and all(
            margin >= 0.0
            and math.isclose(margin, limit - value, abs_tol=1.0e-9)
            for value, limit, margin in zip(spectrum, limits, margins)
        )
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (source_power, load_power, lisn_loss, switching_loss)
        )
        and source_power > 0.0
        and math.isclose(
            source_power,
            load_power + lisn_loss + switching_loss,
            rel_tol=1.0e-9,
            abs_tol=1.0e-12,
        )
        and result_source_power == source_power
        and result_load_power == load_power
        and result_lisn_loss == lisn_loss
        and result_switching_loss == switching_loss
        and bool(str(contract.get("circuit_owner") or ""))
        and contract.get("accepted_circuit_owner") == contract.get("circuit_owner")
        and _is_sha256(str(contract.get("circuit_sha256") or ""))
        and contract.get("accepted_circuit_sha256") == contract.get("circuit_sha256")
        and _is_sha256(str(contract.get("result_sha256") or ""))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )
