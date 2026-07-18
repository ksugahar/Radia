from __future__ import annotations

from test_comsol_generalization_v33 import (
    _summary,
    _with_v33_contact_and_dae_identity,
    gate,
)


_PROMOTED_CASE_IDS = (
    "v34_public_nonlinear_arclength_tangent_branch_turning_point_residual_owner_mismatch",
    "v34_public_electrochemical_species_flux_charge_mass_reaction_boundary_energy_mismatch",
)


def _with_v34_arclength_and_electrochemical_identity(summary: dict) -> dict:
    summary = _with_v33_contact_and_dae_identity(summary)
    generation = "arclength-211"
    summary[
        "nonlinear_arclength_tangent_branch_turning_residual_mesh_result_generation_identity"
    ] = {
        "continuation_generation": generation,
        **{
            key: generation
            for key in (
                "arclength_generation",
                "tangent_generation",
                "branch_generation",
                "turning_generation",
                "residual_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "previous_augmented_state": [1.0, 0.0, 1.0],
        "result_previous_augmented_state": [1.0, 0.0, 1.0],
        "predictor_tangent": [0.6, 0.0, 0.8],
        "result_predictor_tangent": [0.6, 0.0, 0.8],
        "arclength_step": 0.05,
        "result_arclength_step": 0.05,
        "predictor_augmented_state": [1.03, 0.0, 1.04],
        "result_predictor_augmented_state": [1.03, 0.0, 1.04],
        "corrected_augmented_state": [1.03, 0.0, 1.04],
        "result_corrected_augmented_state": [1.03, 0.0, 1.04],
        "branch_id": "upper_branch",
        "result_branch_id": "upper_branch",
        "turning_point_side": "pre_turn_positive_parameter_tangent",
        "result_turning_point_side": "pre_turn_positive_parameter_tangent",
        "corrected_residual_norm": 1.0e-10,
        "result_corrected_residual_norm": 1.0e-10,
        "residual_tolerance": 1.0e-8,
        "result_residual_tolerance": 1.0e-8,
        "continuation_mesh_sha256": "1" * 64,
        "result_continuation_mesh_sha256": "1" * 64,
        "continuation_result_sha256": "2" * 64,
        "accepted_continuation_result_sha256": "2" * 64,
    }
    generation = "electrochemical-211"
    summary[
        "electrochemical_species_flux_charge_mass_reaction_energy_time_mesh_result_generation_identity"
    ] = {
        "electrochemical_generation": generation,
        **{
            key: generation
            for key in (
                "species_generation",
                "flux_generation",
                "charge_generation",
                "mass_generation",
                "reaction_generation",
                "energy_generation",
                "time_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "species_order": ["A_plus", "B_minus", "C_neutral"],
        "result_species_order": ["A_plus", "B_minus", "C_neutral"],
        "charge_numbers": [1, -1, 0],
        "result_charge_numbers": [1, -1, 0],
        "molar_mass_basis": [1.0, 1.0, 1.0],
        "result_molar_mass_basis": [1.0, 1.0, 1.0],
        "reaction_stoichiometry": [-1.0, -1.0, 2.0],
        "result_reaction_stoichiometry": [-1.0, -1.0, 2.0],
        "reaction_extent_mol": 0.5,
        "result_reaction_extent_mol": 0.5,
        "initial_inventory_mol": [1.0, 1.0, 0.0],
        "result_initial_inventory_mol": [1.0, 1.0, 0.0],
        "final_inventory_mol": [0.5, 0.5, 1.0],
        "result_final_inventory_mol": [0.5, 0.5, 1.0],
        "integrated_boundary_flux_mol": [0.0, 0.0, 0.0],
        "result_integrated_boundary_flux_mol": [0.0, 0.0, 0.0],
        "integrated_electric_current_c": 0.0,
        "result_integrated_electric_current_c": 0.0,
        "initial_free_energy_j": 2.0,
        "result_initial_free_energy_j": 2.0,
        "final_free_energy_j": 1.8,
        "result_final_free_energy_j": 1.8,
        "dissipated_free_energy_j": 0.2,
        "result_dissipated_free_energy_j": 0.2,
        "time_s": [0.0, 1.0],
        "result_time_s": [0.0, 1.0],
        "electrochemical_mesh_sha256": "3" * 64,
        "result_electrochemical_mesh_sha256": "3" * 64,
        "electrochemical_result_sha256": "4" * 64,
        "accepted_electrochemical_result_sha256": "4" * 64,
    }
    return summary


def test_v34_public_positive_arclength_and_electrochemical_contracts() -> None:
    result = gate(_with_v34_arclength_and_electrochemical_identity(_summary()))
    assert result["status"] == "ok"


def test_v34_public_nonlinear_arclength_tangent_branch_turning_point_residual_owner_mismatch() -> None:
    summary = _with_v34_arclength_and_electrochemical_identity(_summary())
    summary[
        "nonlinear_arclength_tangent_branch_turning_residual_mesh_result_generation_identity"
    ].update(
        {
            "arclength_generation": "arclength-210",
            "branch_generation": "arclength-209",
            "result_generation": "arclength-208",
            "result_previous_augmented_state": [0.0, 1.0, 1.0],
            "result_predictor_tangent": [-0.8, 0.0, 0.6],
            "result_arclength_step": -0.05,
            "result_predictor_augmented_state": [1.2, 0.0, 0.9],
            "result_corrected_augmented_state": [1.5, 0.0, 0.5],
            "result_branch_id": "lower_branch",
            "result_turning_point_side": "post_turn_negative_parameter_tangent",
            "result_corrected_residual_norm": 1.0e-2,
            "result_residual_tolerance": 1.0e-12,
            "result_continuation_mesh_sha256": "9" * 64,
            "accepted_continuation_result_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_continuation_uses_current_arclength_tangent_branch_turning_residual_mesh_and_result"
    ]


def test_v34_public_electrochemical_species_flux_charge_mass_reaction_boundary_energy_mismatch() -> None:
    summary = _with_v34_arclength_and_electrochemical_identity(_summary())
    summary[
        "electrochemical_species_flux_charge_mass_reaction_energy_time_mesh_result_generation_identity"
    ].update(
        {
            "species_generation": "electrochemical-210",
            "reaction_generation": "electrochemical-209",
            "result_generation": "electrochemical-208",
            "result_species_order": ["C_neutral", "B_minus", "A_plus"],
            "result_charge_numbers": [0, -1, 1],
            "result_molar_mass_basis": [2.0, 1.0, 1.0],
            "result_reaction_stoichiometry": [1.0, -1.0, 0.0],
            "result_reaction_extent_mol": 0.8,
            "result_initial_inventory_mol": [0.0, 1.0, 1.0],
            "result_final_inventory_mol": [2.0, 0.1, 0.1],
            "result_integrated_boundary_flux_mol": [1.0, 0.0, 0.0],
            "result_integrated_electric_current_c": 96485.0,
            "result_initial_free_energy_j": 1.0,
            "result_final_free_energy_j": 2.0,
            "result_dissipated_free_energy_j": -1.0,
            "result_time_s": [1.0, 0.0],
            "result_electrochemical_mesh_sha256": "b" * 64,
            "accepted_electrochemical_result_sha256": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "electrochemical_results_use_current_species_flux_charge_mass_reaction_energy_time_mesh_and_result"
    ]


def test_v34_public_rejects_self_consistent_nonunit_arclength_tangent() -> None:
    summary = _with_v34_arclength_and_electrochemical_identity(_summary())
    identity = summary[
        "nonlinear_arclength_tangent_branch_turning_residual_mesh_result_generation_identity"
    ]
    identity["predictor_tangent"] = [0.6, 0.0, 0.6]
    identity["result_predictor_tangent"] = [0.6, 0.0, 0.6]
    assert gate(summary)["status"] == "needs_attention"


def test_v34_public_rejects_self_consistent_nonstoichiometric_inventory() -> None:
    summary = _with_v34_arclength_and_electrochemical_identity(_summary())
    identity = summary[
        "electrochemical_species_flux_charge_mass_reaction_energy_time_mesh_result_generation_identity"
    ]
    identity["final_inventory_mol"] = [0.4, 0.5, 1.0]
    identity["result_final_inventory_mol"] = [0.4, 0.5, 1.0]
    assert gate(summary)["status"] == "needs_attention"
