from __future__ import annotations

import math

from test_comsol_generalization_v36 import _summary, gate
from test_comsol_generalization_v38 import (
    _with_v38_thermoviscous_and_piezoelectric_identity,
)


_PROMOTED_CASE_IDS = (
    "v39_public_poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mismatch",
    "v39_public_rotating_induction_slip_frequency_current_loss_torque_power_frame_mismatch",
)


def _with_v39_poroelastic_and_induction_identity(summary: dict) -> dict:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(summary)
    generation = "poroelastic-biot-377"
    alpha, pressure, strain, modulus = 0.8, 1.0e5, 1.0e-3, 1.0e9
    permeability, viscosity, gradient = 1.0e-12, 1.0e-3, 1.0e6
    volume, timestep = 1.0e-2, 0.1
    flux = -permeability * gradient / viscosity
    mirrored = {
        "biot_coefficient": alpha, "pore_pressure_pa": pressure,
        "volumetric_strain": strain, "biot_modulus_pa": modulus,
        "fluid_content_increment": alpha * strain + pressure / modulus,
        "permeability_m2": permeability, "dynamic_viscosity_pa_s": viscosity,
        "pressure_gradient_pa_per_m": gradient, "darcy_flux_m_per_s": flux,
        "domain_volume_m3": volume, "time_step_s": timestep,
        "interface_traction_pa": -alpha * pressure,
        "storage_energy_j": 0.5 * pressure * pressure / modulus * volume,
        "skeleton_coupling_work_j": alpha * pressure * strain * volume,
        "fluid_dissipation_j": viscosity / permeability * flux * flux * volume * timestep,
    }
    summary["poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_result_generation_identity"] = {
        "poroelastic_generation": generation,
        **{key: generation for key in (
            "biot_generation", "pressure_generation", "displacement_generation",
            "flux_generation", "storage_generation", "dissipation_generation",
            "interface_generation", "mesh_generation", "result_generation",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "interface_normal": "porous_skeleton_to_free_fluid",
        "result_interface_normal": "porous_skeleton_to_free_fluid",
        "mesh_owner": "comp1/mesh1:poroelastic-377",
        "accepted_mesh_owner": "comp1/mesh1:poroelastic-377",
        "mesh_sha256": "1" * 64, "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }

    generation = "rotating-induction-377"
    frequency, pole_pairs, rotor_speed, torque = 50.0, 2, 150.0, 20.0
    synchronous = 2.0 * math.pi * frequency / pole_pairs
    slip = (synchronous - rotor_speed) / synchronous
    airgap_power, resistance = torque * synchronous, 0.2
    current = math.sqrt((slip * airgap_power) / (3.0 * resistance))
    mirrored = {
        "supply_frequency_hz": frequency, "pole_pairs": pole_pairs,
        "synchronous_speed_rad_per_s": synchronous,
        "rotor_speed_rad_per_s": rotor_speed, "slip": slip,
        "rotor_electrical_frequency_hz": slip * frequency,
        "rotor_phase_current_a_rms": current,
        "rotor_phase_resistance_ohm": resistance,
        "rotor_copper_loss_w": 3.0 * current * current * resistance,
        "airgap_torque_nm": torque, "airgap_power_w": airgap_power,
        "mechanical_power_w": torque * rotor_speed,
    }
    summary["rotating_induction_slip_frequency_current_loss_torque_power_frame_mesh_result_generation_identity"] = {
        "induction_generation": generation,
        **{key: generation for key in (
            "slip_generation", "frequency_generation", "current_generation",
            "loss_generation", "torque_generation", "power_generation",
            "frame_generation", "mesh_generation", "result_generation",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "rotating_frame": "rotor_mechanical_frame",
        "result_rotating_frame": "rotor_mechanical_frame",
        "mesh_owner": "comp1/mesh1:induction-377",
        "accepted_mesh_owner": "comp1/mesh1:induction-377",
        "mesh_sha256": "3" * 64, "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v39_public_positive_poroelastic_and_induction_contracts() -> None:
    assert gate(_with_v39_poroelastic_and_induction_identity(_summary()))["status"] == "ok"


def test_v39_public_poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mismatch() -> None:
    summary = _with_v39_poroelastic_and_induction_identity(_summary())
    identity = summary["poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_result_generation_identity"]
    identity.update({
        "flux_generation": "poroelastic-biot-376", "result_generation": "poroelastic-biot-375",
        "result_fluid_content_increment": -1.0, "result_darcy_flux_m_per_s": 1.0e-3,
        "result_interface_traction_pa": 8.0e4, "result_storage_energy_j": -1.0,
        "result_skeleton_coupling_work_j": -1.0, "result_fluid_dissipation_j": -1.0,
        "result_interface_normal": "free_fluid_to_skeleton",
        "accepted_mesh_owner": "comp1/mesh0:old", "accepted_result_sha256": "9" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["poroelastic_results_use_current_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_and_result"]


def test_v39_public_rotating_induction_slip_frequency_current_loss_torque_power_frame_mismatch() -> None:
    summary = _with_v39_poroelastic_and_induction_identity(_summary())
    identity = summary["rotating_induction_slip_frequency_current_loss_torque_power_frame_mesh_result_generation_identity"]
    identity.update({
        "slip_generation": "rotating-induction-376", "result_generation": "rotating-induction-375",
        "result_slip": -1.0, "result_rotor_electrical_frequency_hz": -50.0,
        "result_rotor_phase_current_a_rms": -1.0, "result_rotor_copper_loss_w": -1.0,
        "result_airgap_torque_nm": -20.0, "result_mechanical_power_w": 100.0,
        "result_rotating_frame": "stator_frame", "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["rotating_induction_results_use_current_slip_frequency_current_loss_torque_power_frame_mesh_and_result"]


def test_v39_public_rejects_self_consistent_wrong_darcy_flux() -> None:
    summary = _with_v39_poroelastic_and_induction_identity(_summary())
    identity = summary["poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_result_generation_identity"]
    identity["darcy_flux_m_per_s"] = identity["result_darcy_flux_m_per_s"] = 1.0e-3
    assert gate(summary)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_induction_power_creation() -> None:
    summary = _with_v39_poroelastic_and_induction_identity(_summary())
    identity = summary["rotating_induction_slip_frequency_current_loss_torque_power_frame_mesh_result_generation_identity"]
    identity["airgap_power_w"] = identity["result_airgap_power_w"] = 1.0
    assert gate(summary)["status"] == "needs_attention"
