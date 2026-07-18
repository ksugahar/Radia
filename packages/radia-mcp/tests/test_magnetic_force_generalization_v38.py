from __future__ import annotations

import math

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v37 import _identity_v37


_PROMOTED_CASE_IDS = (
    "v38_public_harmonic_conductor_skin_proximity_impedance_current_loss_poynting_mesh_mismatch",
    "v38_public_heat_radiation_convection_emissivity_ambient_flux_temperature_energy_mismatch",
)


def _identity_v38():
    identity = _identity_v37()
    generation = "harmonic-conductor-258"
    frequency = 10_000.0
    conductivity = 5.8e7
    skin_depth = math.sqrt(
        2.0
        / (
            2.0
            * math.pi
            * frequency
            * (4.0e-7 * math.pi)
            * conductivity
        )
    )
    voltage = [1.0, 0.5]
    current = [10.0, -2.0]
    denominator = current[0] ** 2 + current[1] ** 2
    impedance = [
        (voltage[0] * current[0] + voltage[1] * current[1]) / denominator,
        (voltage[1] * current[0] - voltage[0] * current[1]) / denominator,
    ]
    loss = 0.5 * (
        voltage[0] * current[0] + voltage[1] * current[1]
    )
    identity[
        "harmonic_conductor_skin_proximity_impedance_current_voltage_loss_poynting_frequency_mesh_owner_result_identity"
    ] = {
        "harmonic_generation": generation,
        **{
            key: generation
            for key in (
                "skin_generation", "proximity_generation",
                "impedance_generation", "current_generation",
                "voltage_generation", "loss_generation",
                "poynting_generation", "frequency_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )
        },
        "frequency_hz": frequency,
        "result_frequency_hz": frequency,
        "conductivity_s_m": conductivity,
        "result_conductivity_s_m": conductivity,
        "relative_permeability": 1.0,
        "result_relative_permeability": 1.0,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "proximity_current_density_peak_a_m2": [1.0e6, 1.2e6, 0.9e6],
        "result_proximity_current_density_peak_a_m2": [1.0e6, 1.2e6, 0.9e6],
        "terminal_voltage_peak_phasor_v": voltage,
        "result_terminal_voltage_peak_phasor_v": voltage,
        "terminal_current_peak_phasor_a": current,
        "result_terminal_current_peak_phasor_a": current,
        "complex_impedance_ohm": impedance,
        "result_complex_impedance_ohm": impedance,
        "copper_loss_w": loss,
        "result_copper_loss_w": loss,
        "inward_poynting_power_w": loss,
        "result_inward_poynting_power_w": loss,
        "phasor_convention": "peak_cosine",
        "result_phasor_convention": "peak_cosine",
        "mesh_levels": [1, 2, 3],
        "result_mesh_levels": [1, 2, 3],
        "impedance_relative_changes": [0.10, 0.02, 0.004],
        "result_impedance_relative_changes": [0.10, 0.02, 0.004],
        "maximum_final_relative_change": 0.01,
        "result_maximum_final_relative_change": 0.01,
        "field_owner": "harmonic:conductor-258",
        "accepted_field_owner": "harmonic:conductor-258",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "harmonic_result_sha256": "2" * 64,
        "accepted_harmonic_result_sha256": "2" * 64,
    }

    generation = "heat-radiation-258"
    sigma = 5.670374419e-8
    coefficient = 10.0
    emissivity = 0.8
    ambient = 300.0
    boundary = 400.0
    convection = coefficient * (boundary - ambient)
    radiation = emissivity * sigma * (boundary**4 - ambient**4)
    conductive = convection + radiation
    area = 0.2 * 0.1
    loss = conductive * area
    identity[
        "heat_convection_radiation_emissivity_ambient_flux_temperature_geometry_mesh_energy_result_identity"
    ] = {
        "heat_generation": generation,
        **{
            key: generation
            for key in (
                "convection_generation", "radiation_generation",
                "emissivity_generation", "ambient_generation",
                "flux_generation", "temperature_generation",
                "geometry_generation", "mesh_generation",
                "energy_generation", "result_generation",
            )
        },
        "convection_coefficient_w_m2_k": coefficient,
        "result_convection_coefficient_w_m2_k": coefficient,
        "emissivity": emissivity,
        "result_emissivity": emissivity,
        "stefan_boltzmann_w_m2_k4": sigma,
        "result_stefan_boltzmann_w_m2_k4": sigma,
        "ambient_temperature_k": ambient,
        "result_ambient_temperature_k": ambient,
        "boundary_temperature_k": boundary,
        "result_boundary_temperature_k": boundary,
        "convection_flux_w_m2": convection,
        "result_convection_flux_w_m2": convection,
        "radiation_flux_w_m2": radiation,
        "result_radiation_flux_w_m2": radiation,
        "conductive_outward_flux_w_m2": conductive,
        "result_conductive_outward_flux_w_m2": conductive,
        "geometry_weighting": "planar_depth",
        "result_geometry_weighting": "planar_depth",
        "boundary_length_m": 0.2,
        "result_boundary_length_m": 0.2,
        "planar_depth_m": 0.1,
        "result_planar_depth_m": 0.1,
        "effective_boundary_area_m2": area,
        "result_effective_boundary_area_m2": area,
        "boundary_heat_loss_w": loss,
        "result_boundary_heat_loss_w": loss,
        "domain_heat_generation_w": loss,
        "result_domain_heat_generation_w": loss,
        "energy_balance_residual_w": 0.0,
        "result_energy_balance_residual_w": 0.0,
        "energy_tolerance_w": 1.0e-9,
        "result_energy_tolerance_w": 1.0e-9,
        "mesh_owner": "heat:mesh-258",
        "accepted_mesh_owner": "heat:mesh-258",
        "heat_result_sha256": "3" * 64,
        "accepted_heat_result_sha256": "3" * 64,
    }
    return identity


def test_v38_public_positive_harmonic_conductor_and_heat_boundary_closure():
    assert _gate(_identity_v38())["status"] == "ok"


def test_v38_public_harmonic_conductor_skin_proximity_impedance_current_loss_poynting_mesh_mismatch():
    identity = _identity_v38()
    row = identity[
        "harmonic_conductor_skin_proximity_impedance_current_voltage_loss_poynting_frequency_mesh_owner_result_identity"
    ]
    row.update(
        {
            "skin_generation": "harmonic-conductor-257",
            "result_frequency_hz": -10_000.0,
            "result_skin_depth_m": -1.0,
            "result_complex_impedance_ohm": [-1.0, 2.0],
            "result_copper_loss_w": -4.5,
            "result_inward_poynting_power_w": 45.0,
            "result_mesh_levels": [3, 2, 1],
            "accepted_field_owner": "stale:harmonic",
        }
    )
    assert _gate(identity)["status"] == "needs_attention"


def test_v38_public_heat_radiation_convection_emissivity_ambient_flux_temperature_energy_mismatch():
    identity = _identity_v38()
    row = identity[
        "heat_convection_radiation_emissivity_ambient_flux_temperature_geometry_mesh_energy_result_identity"
    ]
    row.update(
        {
            "radiation_generation": "heat-radiation-257",
            "result_emissivity": 1.5,
            "result_ambient_temperature_k": -300.0,
            "result_radiation_flux_w_m2": -1.0,
            "result_geometry_weighting": "axisymmetric",
            "result_boundary_heat_loss_w": -1.0,
            "result_energy_balance_residual_w": 11.0,
            "accepted_mesh_owner": "stale:heat",
        }
    )
    assert _gate(identity)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_but_wrong_skin_depth():
    identity = _identity_v38()
    row = identity[
        "harmonic_conductor_skin_proximity_impedance_current_voltage_loss_poynting_frequency_mesh_owner_result_identity"
    ]
    row["skin_depth_m"] *= 2.0
    row["result_skin_depth_m"] = row["skin_depth_m"]
    assert _gate(identity)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_but_wrong_radiation_flux():
    identity = _identity_v38()
    row = identity[
        "heat_convection_radiation_emissivity_ambient_flux_temperature_geometry_mesh_energy_result_identity"
    ]
    row["radiation_flux_w_m2"] *= 0.5
    row["result_radiation_flux_w_m2"] = row["radiation_flux_w_m2"]
    assert _gate(identity)["status"] == "needs_attention"
