from copy import deepcopy

from radia_mcp.radia_ngsolve.dissipation_reaction_identity_v55 import (
    ELECTROCHEM,
    THERMO,
    validate_public_v55_identity,
)


CASE_IDS = {
    "v55_public_thermoelastic_damping_complexeigenfrequency_energy_dissipation_normalization_owner_mismatch",
    "v55_public_electrochem_current_species_stoichiometry_boundaryflux_time_owner_mismatch",
}


def _records() -> dict[str, object]:
    thermo_generation = "comsol-public-thermo-v55"
    electrochem_generation = "comsol-public-electrochem-v55"
    frequency_hz = 12500.0
    decay_hz = -2.5
    quality_factor = frequency_hz / (-2.0 * decay_hz)
    stored_energy_j = 4.0e-6
    cycle_dissipation_j = 2.0 * 3.141592653589793 * stored_energy_j / quality_factor
    rates = {"Li": -1.0e-6, "Li_plus": 1.0e-6, "electron": 1.0e-6}
    charges = {"Li": 0, "Li_plus": 1, "electron": -1}
    return {
        THERMO: {
            "generation": thermo_generation,
            **{
                name: thermo_generation
                for name in (
                    "eigenfrequency_generation",
                    "energy_generation",
                    "dissipation_generation",
                    "normalization_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "complex_eigenfrequency_hz": [frequency_hz, decay_hz],
            "result_complex_eigenfrequency_hz": [frequency_hz, decay_hz],
            "stored_energy_j": stored_energy_j,
            "result_stored_energy_j": stored_energy_j,
            "cycle_dissipation_j": cycle_dissipation_j,
            "result_cycle_dissipation_j": cycle_dissipation_j,
            "quality_factor": quality_factor,
            "result_quality_factor": quality_factor,
            "modal_normalization": "unit_total_stored_energy",
            "result_modal_normalization": "unit_total_stored_energy",
            "solution_owner": "solution:thermoelastic-v55",
            "result_solution_owner": "solution:thermoelastic-v55",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        ELECTROCHEM: {
            "generation": electrochem_generation,
            **{
                name: electrochem_generation
                for name in (
                    "current_generation",
                    "species_generation",
                    "stoichiometry_generation",
                    "flux_generation",
                    "time_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "terminal_current_a": 0.09648533212,
            "result_terminal_current_a": 0.09648533212,
            "species_rate_mol_s": rates,
            "result_species_rate_mol_s": rates,
            "species_charge_number": charges,
            "result_species_charge_number": charges,
            "boundary_species_flux_mol_s": {"Li_plus": -1.0e-6},
            "result_boundary_species_flux_mol_s": {"Li_plus": -1.0e-6},
            "time_s": 0.5,
            "result_time_s": 0.5,
            "solution_owner": "solution:electrochem-v55",
            "result_solution_owner": "solution:electrochem-v55",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
    }


def test_v55_public_positive_replay_is_accepted() -> None:
    assert validate_public_v55_identity(_records())["status"] == "ok"


def test_v55_frozen_public_mutations_are_rejected() -> None:
    value = deepcopy(_records())
    value[THERMO].update(
        {
            "result_complex_eigenfrequency_hz": [12500.0, 2.5],
            "result_stored_energy_j": 8.0e-6,
            "result_cycle_dissipation_j": 0.0,
            "result_modal_normalization": "peak_displacement",
            "result_solution_owner": "solution:stale",
        }
    )
    value[ELECTROCHEM].update(
        {
            "result_terminal_current_a": 1.0,
            "result_species_rate_mol_s": {"Li_plus": -1.0},
            "result_species_charge_number": {"Li_plus": -1},
            "result_boundary_species_flux_mol_s": {"Li_plus": 1.0},
            "result_time_s": 0.25,
            "result_solution_owner": "solution:stale",
        }
    )
    assert validate_public_v55_identity(value)["status"] == "needs_attention"


def test_v55_self_consistent_nonphysical_damping_is_rejected() -> None:
    value = deepcopy(_records())
    value[THERMO]["complex_eigenfrequency_hz"] = value[THERMO][
        "result_complex_eigenfrequency_hz"
    ] = [12500.0, 2.5]
    value[THERMO]["cycle_dissipation_j"] = value[THERMO][
        "result_cycle_dissipation_j"
    ] = 0.0
    assert validate_public_v55_identity(value)["status"] == "needs_attention"


def test_v55_self_consistent_charge_or_flux_imbalance_is_rejected() -> None:
    value = deepcopy(_records())
    bad_rates = {"Li": -1.0e-6, "Li_plus": 1.0e-6, "electron": 0.5e-6}
    value[ELECTROCHEM]["species_rate_mol_s"] = value[ELECTROCHEM][
        "result_species_rate_mol_s"
    ] = bad_rates
    value[ELECTROCHEM]["boundary_species_flux_mol_s"] = value[ELECTROCHEM][
        "result_boundary_species_flux_mol_s"
    ] = {"Li_plus": 1.0e-6}
    assert validate_public_v55_identity(value)["status"] == "needs_attention"


def test_v55_numeric_sha256_values_are_rejected() -> None:
    value = _records()
    numeric_digest = int("9" * 64)
    for row in value.values():
        row["result_sha256"] = numeric_digest
        row["accepted_result_sha256"] = numeric_digest
    assert validate_public_v55_identity(value)["status"] == "needs_attention"
