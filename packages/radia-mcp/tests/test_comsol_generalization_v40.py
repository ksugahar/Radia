from __future__ import annotations

import math

from test_comsol_generalization_v36 import _summary, gate
from test_comsol_generalization_v39 import _with_v39_poroelastic_and_induction_identity


_PROMOTED_CASE_IDS = (
    "v40_public_thermoacoustic_meanflow_convected_wavenumber_flux_impedance_power_mismatch",
    "v40_public_battery_electrothermal_soc_current_heat_temperature_energy_safety_mismatch",
)
_THERMOACOUSTIC_KEY = (
    "thermoacoustic_meanflow_convected_wavenumber_flux_impedance_power_mesh_"
    "result_generation_identity"
)
_BATTERY_KEY = (
    "battery_electrothermal_soc_current_heat_temperature_energy_safety_mesh_"
    "result_generation_identity"
)


def _with_v40_thermoacoustic_and_battery_identity(summary: dict) -> dict:
    summary = _with_v39_poroelastic_and_induction_identity(summary)
    generation = "thermoacoustic-719"
    frequency, sound_speed, mach = 1000.0, 343.0, 0.1
    density, pressure, area = 1.2, 1.0, 0.1
    flow_speed = mach * sound_speed
    particle_velocity = pressure / (density * sound_speed)
    intensity = pressure * particle_velocity
    power = intensity * area
    mirrored = {
        "frequency_hz": frequency, "sound_speed_m_per_s": sound_speed,
        "mean_flow_mach": mach, "mean_flow_speed_m_per_s": flow_speed,
        "downstream_wavenumber_rad_per_m": 2.0 * math.pi * frequency / (sound_speed + flow_speed),
        "upstream_wavenumber_rad_per_m": 2.0 * math.pi * frequency / (sound_speed - flow_speed),
        "density_kg_per_m3": density, "pressure_rms_pa": pressure,
        "particle_velocity_rms_m_per_s": particle_velocity,
        "acoustic_intensity_w_per_m2": intensity, "boundary_area_m2": area,
        "boundary_impedance_pa_s_per_m": density * sound_speed,
        "boundary_flux_power_w": power, "impedance_work_w": power,
        "dissipated_power_w": power, "power_balance_residual_w": 0.0,
    }
    summary[_THERMOACOUSTIC_KEY] = {
        "thermoacoustic_generation": generation,
        **{key: generation for key in (
            "meanflow_generation", "wavenumber_generation", "flux_generation",
            "impedance_generation", "power_generation", "mesh_generation",
            "result_generation",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "comp1/mesh1:thermoacoustic-719",
        "accepted_mesh_owner": "comp1/mesh1:thermoacoustic-719",
        "mesh_sha256": "1" * 64, "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }

    generation = "battery-electrothermal-719"
    capacity, initial_soc, current, voltage = 3600.0, 0.8, 2.0, 3.7
    timestep, resistance = 10.0, 0.05
    irreversible_heat, reversible_heat = current * current * resistance * timestep, 0.5
    thermal_energy = irreversible_heat + reversible_heat
    mass, heat_capacity, initial_temperature = 0.05, 1000.0, 298.15
    mirrored = {
        "capacity_c": capacity, "initial_state_of_charge": initial_soc,
        "terminal_current_a": current, "terminal_voltage_v": voltage,
        "time_step_s": timestep,
        "final_state_of_charge": initial_soc - current * timestep / capacity,
        "internal_resistance_ohm": resistance,
        "irreversible_heat_j": irreversible_heat, "reversible_heat_j": reversible_heat,
        "thermal_energy_j": thermal_energy,
        "electrical_energy_j": voltage * current * timestep,
        "cell_mass_kg": mass, "specific_heat_j_per_kg_k": heat_capacity,
        "initial_temperature_k": initial_temperature,
        "final_temperature_k": initial_temperature + thermal_energy / (mass * heat_capacity),
        "maximum_safe_temperature_k": 333.15, "thermal_balance_residual_j": 0.0,
    }
    summary[_BATTERY_KEY] = {
        "battery_generation": generation,
        **{key: generation for key in (
            "soc_generation", "current_generation", "heat_generation",
            "temperature_generation", "energy_generation", "safety_generation",
            "mesh_generation", "result_generation",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "comp1/mesh1:battery-719",
        "accepted_mesh_owner": "comp1/mesh1:battery-719",
        "mesh_sha256": "3" * 64, "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v40_public_positive_thermoacoustic_and_battery_contracts() -> None:
    assert gate(_with_v40_thermoacoustic_and_battery_identity(_summary()))["status"] == "ok"


def test_v40_public_thermoacoustic_meanflow_convected_wavenumber_flux_impedance_power_mismatch() -> None:
    summary = _with_v40_thermoacoustic_and_battery_identity(_summary())
    identity = summary[_THERMOACOUSTIC_KEY]
    identity.update({
        "wavenumber_generation": "thermoacoustic-718",
        "result_generation": "thermoacoustic-717",
        "result_mean_flow_mach": -2.0,
        "result_downstream_wavenumber_rad_per_m": -1.0,
        "result_acoustic_intensity_w_per_m2": -1.0,
        "result_boundary_impedance_pa_s_per_m": -1.0,
        "result_boundary_flux_power_w": -1.0,
        "result_power_balance_residual_w": 9.0,
        "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["thermoacoustic_results_use_current_meanflow_wavenumber_flux_impedance_power_mesh_and_result"]


def test_v40_public_battery_electrothermal_soc_current_heat_temperature_energy_safety_mismatch() -> None:
    summary = _with_v40_thermoacoustic_and_battery_identity(_summary())
    identity = summary[_BATTERY_KEY]
    identity.update({
        "heat_generation": "battery-electrothermal-718",
        "result_generation": "battery-electrothermal-717",
        "result_final_state_of_charge": 2.0,
        "result_terminal_current_a": -2.0,
        "result_irreversible_heat_j": -1.0,
        "result_thermal_energy_j": -1.0,
        "result_final_temperature_k": 500.0,
        "result_thermal_balance_residual_j": 9.0,
        "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "b" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["battery_results_use_current_soc_current_heat_temperature_energy_safety_mesh_and_result"]


def test_v40_public_rejects_self_consistent_nonconvected_wavenumber() -> None:
    summary = _with_v40_thermoacoustic_and_battery_identity(_summary())
    identity = summary[_THERMOACOUSTIC_KEY]
    value = 2.0 * math.pi * identity["frequency_hz"] / identity["sound_speed_m_per_s"]
    identity["downstream_wavenumber_rad_per_m"] = value
    identity["result_downstream_wavenumber_rad_per_m"] = value
    assert gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_battery_heat_creation() -> None:
    summary = _with_v40_thermoacoustic_and_battery_identity(_summary())
    identity = summary[_BATTERY_KEY]
    identity["thermal_energy_j"] = identity["result_thermal_energy_j"] = 1.0
    identity["thermal_balance_residual_j"] = identity["result_thermal_balance_residual_j"] = 0.0
    assert gate(summary)["status"] == "needs_attention"
