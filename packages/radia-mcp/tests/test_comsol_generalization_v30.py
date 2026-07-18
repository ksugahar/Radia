from __future__ import annotations

from test_comsol_generalization_v29 import (
    _summary,
    _with_v29_sliding_and_radiation_identity,
    gate,
)


_PROMOTED_CASE_IDS = (
    "v30_public_joule_heat_source_mapping_resistivity_temperature_time_average_energy_balance_mismatch",
    "v30_public_nonlinear_eigenmode_mac_branch_normalization_parameter_continuation_mismatch",
)


def _with_v30_joule_and_eigenmode_identity(summary: dict) -> dict:
    summary = _with_v29_sliding_and_radiation_identity(summary)
    generation = "joule-heat-closure-171"
    summary[
        "joule_heat_source_current_density_resistivity_temperature_frame_time_average_energy_mesh_result_generation_identity"
    ] = {
        "joule_generation": generation,
        "mapping_joule_generation": generation,
        "resistivity_joule_generation": generation,
        "temperature_joule_generation": generation,
        "frame_joule_generation": generation,
        "averaging_joule_generation": generation,
        "energy_joule_generation": generation,
        "mesh_joule_generation": generation,
        "result_joule_generation": generation,
        "current_density_field_id": "ec.J-current-171",
        "result_current_density_field_id": "ec.J-current-171",
        "temperature_field_id": "ht.T-current-171",
        "result_temperature_field_id": "ht.T-current-171",
        "resistivity_model_id": "rho(T)-copper-171",
        "result_resistivity_model_id": "rho(T)-copper-171",
        "source_frame": "material_spatial",
        "result_source_frame": "material_spatial",
        "averaging_window_s": [0.02, 0.04],
        "result_averaging_window_s": [0.02, 0.04],
        "time_average_method": "trapezoidal_period_average",
        "result_time_average_method": "trapezoidal_period_average",
        "electric_loss_w": 12.5,
        "result_electric_loss_w": 12.5,
        "heat_source_integral_w": 12.5,
        "result_heat_source_integral_w": 12.5,
        "energy_balance_relative_tolerance": 1.0e-9,
        "coupled_mesh_sha256": "1" * 64,
        "result_coupled_mesh_sha256": "1" * 64,
        "joule_heat_result_sha256": "2" * 64,
        "accepted_joule_heat_result_sha256": "2" * 64,
    }
    generation = "nonlinear-eigenmode-171"
    summary[
        "nonlinear_eigenmode_continuation_parameter_normalization_phase_mac_branch_eigenvalue_mesh_result_generation_identity"
    ] = {
        "eigenmode_generation": generation,
        "continuation_eigenmode_generation": generation,
        "normalization_eigenmode_generation": generation,
        "phase_eigenmode_generation": generation,
        "mac_eigenmode_generation": generation,
        "branch_eigenmode_generation": generation,
        "eigenvalue_eigenmode_generation": generation,
        "mesh_eigenmode_generation": generation,
        "result_eigenmode_generation": generation,
        "continuation_parameter_name": "prestress_scale",
        "result_continuation_parameter_name": "prestress_scale",
        "continuation_parameter_values": [0.0, 0.5, 1.0],
        "result_continuation_parameter_values": [0.0, 0.5, 1.0],
        "mode_normalization": "unit_mass",
        "result_mode_normalization": "unit_mass",
        "phase_anchor_dof": "tip-z",
        "result_phase_anchor_dof": "tip-z",
        "mac_reference_branch_ids": [1, 2],
        "result_mac_reference_branch_ids": [1, 2],
        "mode_branch_ids": [[1, 2], [1, 2], [1, 2]],
        "result_mode_branch_ids": [[1, 2], [1, 2], [1, 2]],
        "eigenvalues_ri": [
            [[100.0, 0.0], [150.0, 0.0]],
            [[102.0, 0.0], [148.0, 0.0]],
            [[105.0, 0.0], [145.0, 0.0]],
        ],
        "result_eigenvalues_ri": [
            [[100.0, 0.0], [150.0, 0.0]],
            [[102.0, 0.0], [148.0, 0.0]],
            [[105.0, 0.0], [145.0, 0.0]],
        ],
        "mac_assignment_sha256": "3" * 64,
        "result_mac_assignment_sha256": "3" * 64,
        "eigenmode_mesh_sha256": "4" * 64,
        "result_eigenmode_mesh_sha256": "4" * 64,
        "eigenmode_result_sha256": "5" * 64,
        "accepted_eigenmode_result_sha256": "5" * 64,
    }
    return summary


def test_v30_public_positive_joule_heat_and_nonlinear_eigenmode() -> None:
    assert gate(_with_v30_joule_and_eigenmode_identity(_summary()))["status"] == "ok"


def test_v30_public_joule_heat_source_mapping_resistivity_temperature_time_average_energy_balance_mismatch() -> None:
    summary = _with_v30_joule_and_eigenmode_identity(_summary())
    summary[
        "joule_heat_source_current_density_resistivity_temperature_frame_time_average_energy_mesh_result_generation_identity"
    ].update(
        {
            "mapping_joule_generation": "joule-heat-closure-170",
            "averaging_joule_generation": "joule-heat-closure-169",
            "result_joule_generation": "joule-heat-closure-168",
            "result_current_density_field_id": "ec.J-old",
            "result_temperature_field_id": "ht.T-old",
            "result_resistivity_model_id": "rho-constant",
            "result_source_frame": "global_spatial",
            "result_averaging_window_s": [0.0, 0.01],
            "result_time_average_method": "final_sample",
            "result_electric_loss_w": 10.0,
            "result_heat_source_integral_w": 20.0,
            "result_coupled_mesh_sha256": "9" * 64,
            "accepted_joule_heat_result_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "joule_heat_uses_current_mapping_resistivity_temperature_frame_average_energy_mesh_and_result"
    ]


def test_v30_public_nonlinear_eigenmode_mac_branch_normalization_parameter_continuation_mismatch() -> None:
    summary = _with_v30_joule_and_eigenmode_identity(_summary())
    summary[
        "nonlinear_eigenmode_continuation_parameter_normalization_phase_mac_branch_eigenvalue_mesh_result_generation_identity"
    ].update(
        {
            "continuation_eigenmode_generation": "nonlinear-eigenmode-170",
            "mac_eigenmode_generation": "nonlinear-eigenmode-169",
            "result_eigenmode_generation": "nonlinear-eigenmode-168",
            "result_continuation_parameter_name": "temperature",
            "result_continuation_parameter_values": [1.0, 0.5, 0.0],
            "result_mode_normalization": "unit_max",
            "result_phase_anchor_dof": "base-x",
            "result_mac_reference_branch_ids": [2, 1],
            "result_mode_branch_ids": [[2, 1], [1, 2]],
            "result_eigenvalues_ri": [[[100.0, 0.0], [151.0, 0.0]]],
            "result_mac_assignment_sha256": "b" * 64,
            "result_eigenmode_mesh_sha256": "c" * 64,
            "accepted_eigenmode_result_sha256": "d" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_eigenmodes_use_current_continuation_normalization_phase_mac_branch_eigenvalues_mesh_and_result"
    ]
