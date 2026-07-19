from __future__ import annotations

import math

from test_comsol_generalization_v36 import _summary, gate
from test_comsol_generalization_v41 import _with_v41_piezoelectric_and_bearing_identity


_PROMOTED_CASE_IDS = (
    "v42_public_inductionheating_skin_proximity_joule_thermal_flux_temperature_energy_mismatch",
    "v42_public_species_transport_reaction_diffusion_flux_massbalance_rate_temperature_mismatch",
)
_INDUCTION_KEY = (
    "inductionheating_skin_proximity_joule_thermal_flux_temperature_energy_"
    "mesh_result_generation_identity"
)
_SPECIES_KEY = (
    "species_transport_reaction_diffusion_flux_massbalance_rate_temperature_"
    "mesh_result_generation_identity"
)


def _with_v42_induction_and_species_identity(summary: dict) -> dict:
    summary = _with_v41_piezoelectric_and_bearing_identity(summary)
    generation = "induction-heating-725"
    frequency = 10_000.0
    conductivity = 5.8e7
    permeability = 4.0e-7 * math.pi
    joule_loss, magnetic_loss = 780.0, 20.0
    input_power = joule_loss + magnetic_loss
    ambient, maximum = 293.15, 373.15
    mirrored = {
        "frequency_hz": frequency,
        "conductivity_s_per_m": conductivity,
        "relative_permeability": 1.0,
        "skin_depth_m": math.sqrt(2.0 / (2.0 * math.pi * frequency * permeability * conductivity)),
        "surface_current_density_a_per_m": [1200.0, 1500.0, 1100.0],
        "proximity_current_density_a_per_m2": [4.0e7, 6.0e7, 3.5e7],
        "joule_loss_w": joule_loss,
        "magnetic_loss_w": magnetic_loss,
        "electromagnetic_input_power_w": input_power,
        "outward_thermal_flux_w": input_power,
        "ambient_temperature_k": ambient,
        "maximum_temperature_k": maximum,
        "temperature_rise_k": maximum - ambient,
        "electromagnetic_power_balance_residual_w": 0.0,
        "thermal_power_balance_residual_w": 0.0,
    }
    summary[_INDUCTION_KEY] = {
        "induction_generation": generation,
        **{key: generation for key in (
            "skin_generation", "proximity_generation", "joule_generation",
            "thermal_generation", "temperature_generation", "energy_generation",
            "mesh_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "component/mesh:induction-725",
        "accepted_mesh_owner": "component/mesh:induction-725",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "reacting-species-725"
    temperature, gas_constant = 350.0, 8.314462618
    preexponential, activation_energy = 1.0e6, 40_000.0
    rate = preexponential * math.exp(-activation_energy / (gas_constant * temperature))
    concentration, volume = 2.0, 0.01
    consumption = rate * concentration * volume
    mirrored = {
        "diffusivity_m2_per_s": 2.0e-9,
        "temperature_k": temperature,
        "gas_constant_j_per_mol_k": gas_constant,
        "preexponential_factor_per_s": preexponential,
        "activation_energy_j_per_mol": activation_energy,
        "reaction_rate_constant_per_s": rate,
        "mean_concentration_mol_per_m3": concentration,
        "domain_volume_m3": volume,
        "integrated_species_mol": concentration * volume,
        "integrated_consumption_mol_per_s": consumption,
        "inward_boundary_flux_mol_per_s": consumption,
        "mass_balance_residual_mol_per_s": 0.0,
    }
    summary[_SPECIES_KEY] = {
        "species_generation": generation,
        **{key: generation for key in (
            "diffusion_generation", "reaction_generation", "flux_generation",
            "mass_generation", "temperature_generation", "mesh_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "component/mesh:species-725",
        "accepted_mesh_owner": "component/mesh:species-725",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v42_public_positive_induction_and_species_contracts() -> None:
    assert gate(_with_v42_induction_and_species_identity(_summary()))["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v42_public_inductionheating_skin_proximity_joule_thermal_flux_temperature_energy_mismatch() -> None:
    summary = _with_v42_induction_and_species_identity(_summary())
    summary[_INDUCTION_KEY].update({
        "skin_generation": "induction-heating-724",
        "result_generation": "induction-heating-723",
        "result_skin_depth_m": -1.0,
        "result_proximity_current_density_a_per_m2": [-1.0],
        "result_joule_loss_w": -10.0,
        "result_outward_thermal_flux_w": 10.0,
        "result_maximum_temperature_k": 250.0,
        "result_electromagnetic_power_balance_residual_w": 99.0,
        "accepted_mesh_owner": "component/mesh:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["induction_heating_results_use_current_skin_proximity_joule_thermal_temperature_energy_mesh_and_result"]


def test_v42_public_species_transport_reaction_diffusion_flux_massbalance_rate_temperature_mismatch() -> None:
    summary = _with_v42_induction_and_species_identity(_summary())
    summary[_SPECIES_KEY].update({
        "reaction_generation": "reacting-species-724",
        "result_generation": "reacting-species-723",
        "result_diffusivity_m2_per_s": -1.0,
        "result_temperature_k": 0.0,
        "result_reaction_rate_constant_per_s": -1.0,
        "result_integrated_species_mol": -1.0,
        "result_inward_boundary_flux_mol_per_s": -1.0,
        "result_mass_balance_residual_mol_per_s": 1.0,
        "accepted_mesh_owner": "component/mesh:old",
        "accepted_result_sha256": "b" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["reacting_species_results_use_current_diffusion_rate_flux_mass_temperature_mesh_and_result"]


def test_v42_public_rejects_self_consistent_wrong_skin_depth() -> None:
    summary = _with_v42_induction_and_species_identity(_summary())
    identity = summary[_INDUCTION_KEY]
    wrong = 2.0 * identity["skin_depth_m"]
    identity["skin_depth_m"] = identity["result_skin_depth_m"] = wrong
    assert gate(summary)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_species_mass_leak() -> None:
    summary = _with_v42_induction_and_species_identity(_summary())
    identity = summary[_SPECIES_KEY]
    consumption = identity["integrated_consumption_mol_per_s"]
    identity["inward_boundary_flux_mol_per_s"] = 2.0 * consumption
    identity["result_inward_boundary_flux_mol_per_s"] = 2.0 * consumption
    identity["mass_balance_residual_mol_per_s"] = consumption
    identity["result_mass_balance_residual_mol_per_s"] = consumption
    assert gate(summary)["status"] == "needs_attention"
