from copy import deepcopy

from radia_mcp.radia_ngsolve.conservation_identity_v54 import PIEZO, SPECIES, validate_public_v54_identity


CASE_IDS = {
    "v54_public_piezoelectric_energy_reciprocity_voltage_charge_work_phase_owner_mismatch",
    "v54_public_reactingflow_species_massfraction_rate_flux_time_solution_owner_mismatch",
}


def _records() -> dict[str, object]:
    piezo_generation = "comsol-public-piezo-v54"
    species_generation = "comsol-public-species-v54"
    fractions = {"H2": 0.10, "O2": 0.20, "H2O": 0.15, "N2": 0.55}
    rates = {"H2": -2.0, "O2": -1.0, "H2O": 2.0, "N2": 0.0}
    return {
        PIEZO: {
            "generation": piezo_generation,
            **{name: piezo_generation for name in ("electric_generation", "mechanical_generation", "phase_generation", "owner_generation", "result_generation")},
            "voltage_v": 10.0, "result_voltage_v": 10.0,
            "charge_c": 2.0e-6, "result_charge_c": 2.0e-6,
            "electric_work_j": 1.0e-5, "result_electric_work_j": 1.0e-5,
            "mechanical_work_j": 1.0e-5, "result_mechanical_work_j": 1.0e-5,
            "harmonic_phase_rad": 0.25, "result_harmonic_phase_rad": 0.25,
            "solution_owner": "solution:piezo-v54", "result_solution_owner": "solution:piezo-v54",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        SPECIES: {
            "generation": species_generation,
            **{name: species_generation for name in ("fraction_generation", "rate_generation", "flux_generation", "time_generation", "owner_generation", "result_generation")},
            "species_mass_fraction": fractions, "result_species_mass_fraction": fractions,
            "stoichiometric_rate_mol_m3_s": rates, "result_stoichiometric_rate_mol_m3_s": rates,
            "total_species_mass_rate_kg_s": 2.0e-5, "result_total_species_mass_rate_kg_s": 2.0e-5,
            "boundary_mass_flux_kg_s": -2.0e-5, "result_boundary_mass_flux_kg_s": -2.0e-5,
            "time_s": 0.02, "result_time_s": 0.02,
            "solution_owner": "solution:reacting-v54", "result_solution_owner": "solution:reacting-v54",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def test_v54_public_positive_replay_is_accepted() -> None:
    assert validate_public_v54_identity(_records())["status"] == "ok"


def test_v54_frozen_public_mutations_are_rejected() -> None:
    value = deepcopy(_records())
    value[PIEZO].update({"result_voltage_v": 5.0, "result_electric_work_j": 2.0e-5, "result_mechanical_work_j": 8.0e-6, "result_harmonic_phase_rad": -0.25, "result_solution_owner": "solution:stale"})
    value[SPECIES].update({"result_species_mass_fraction": {"H2": 0.8, "O2": 0.8}, "result_boundary_mass_flux_kg_s": 2.0e-5, "result_time_s": 0.01, "result_solution_owner": "solution:stale"})
    assert validate_public_v54_identity(value)["status"] == "needs_attention"


def test_v54_self_consistent_conservation_errors_are_rejected() -> None:
    value = deepcopy(_records())
    value[PIEZO]["electric_work_j"] = value[PIEZO]["result_electric_work_j"] = 2.0e-5
    value[SPECIES]["species_mass_fraction"] = value[SPECIES]["result_species_mass_fraction"] = {"H2": 0.8, "O2": 0.8}
    value[SPECIES]["boundary_mass_flux_kg_s"] = value[SPECIES]["result_boundary_mass_flux_kg_s"] = 2.0e-5
    assert validate_public_v54_identity(value)["status"] == "needs_attention"
