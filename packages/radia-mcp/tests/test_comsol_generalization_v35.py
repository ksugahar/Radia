from __future__ import annotations

from test_comsol_generalization_v34 import (
    _summary,
    _with_v34_arclength_and_electrochemical_identity,
    gate,
)


_PROMOTED_CASE_IDS = (
    "v35_public_multirate_electromechanical_event_interpolation_work_power_timegrid_mismatch",
    "v35_public_adjoint_sensitivity_objective_chainrule_constraint_fd_mesh_owner_mismatch",
)


def _with_v35_multirate_and_adjoint_identity(summary: dict) -> dict:
    summary = _with_v34_arclength_and_electrochemical_identity(summary)
    generation = "multirate-coupling-221"
    summary[
        "multirate_electromechanical_event_interpolation_work_power_timegrid_frame_mesh_result_generation_identity"
    ] = {
        "coupling_generation": generation,
        **{
            key: generation
            for key in (
                "electrical_generation",
                "mechanical_generation",
                "event_generation",
                "timegrid_generation",
                "power_generation",
                "work_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "electrical_time_s": [0.0, 0.0005, 0.001],
        "result_electrical_time_s": [0.0, 0.0005, 0.001],
        "mechanical_time_s": [0.0, 0.00025, 0.0005, 0.00075, 0.001],
        "result_mechanical_time_s": [0.0, 0.00025, 0.0005, 0.00075, 0.001],
        "event_time_s": 0.0005,
        "result_event_time_s": 0.0005,
        "event_interpolation_side": "right_continuous_after_event",
        "result_event_interpolation_side": "right_continuous_after_event",
        "substep_owner": "coupler:electrical2_mechanical4",
        "result_substep_owner": "coupler:electrical2_mechanical4",
        "coordinate_frame": "stationary_xyz",
        "result_coordinate_frame": "stationary_xyz",
        "electrical_input_energy_j": 0.012,
        "result_electrical_input_energy_j": 0.012,
        "mechanical_output_work_j": 0.009,
        "result_mechanical_output_work_j": 0.009,
        "dissipated_energy_j": 0.003,
        "result_dissipated_energy_j": 0.003,
        "energy_balance_tolerance_j": 1.0e-10,
        "result_energy_balance_tolerance_j": 1.0e-10,
        "coupling_mesh_sha256": "1" * 64,
        "result_coupling_mesh_sha256": "1" * 64,
        "coupling_result_sha256": "2" * 64,
        "accepted_coupling_result_sha256": "2" * 64,
    }
    generation = "adjoint-sensitivity-221"
    summary[
        "adjoint_objective_design_chainrule_constraint_fd_mesh_solution_gradient_generation_identity"
    ] = {
        "sensitivity_generation": generation,
        **{
            key: generation
            for key in (
                "objective_generation",
                "design_generation",
                "chainrule_generation",
                "constraint_generation",
                "fd_generation",
                "mesh_generation",
                "solution_generation",
                "result_generation",
            )
        },
        "objective_tag": "torque_ripple_rms",
        "result_objective_tag": "torque_ripple_rms",
        "design_variable": "magnet_arc_rad",
        "result_design_variable": "magnet_arc_rad",
        "design_scale": 0.1,
        "result_design_scale": 0.1,
        "active_constraint": "magnet_volume_constant",
        "result_active_constraint": "magnet_volume_constant",
        "adjoint_gradient": 2.5,
        "chainrule_gradient": 2.5,
        "finite_difference_gradient": 2.500001,
        "gradient_tolerance": 1.0e-4,
        "result_gradient_tolerance": 1.0e-4,
        "fd_perturbation": 1.0e-5,
        "result_fd_perturbation": 1.0e-5,
        "sensitivity_mesh_sha256": "3" * 64,
        "result_sensitivity_mesh_sha256": "3" * 64,
        "primal_solution_sha256": "4" * 64,
        "result_primal_solution_sha256": "4" * 64,
        "gradient_result_sha256": "5" * 64,
        "accepted_gradient_result_sha256": "5" * 64,
    }
    return summary


def test_v35_public_positive_multirate_and_adjoint_contracts() -> None:
    assert gate(_with_v35_multirate_and_adjoint_identity(_summary()))["status"] == "ok"


def test_v35_public_multirate_electromechanical_event_interpolation_work_power_timegrid_mismatch() -> None:
    summary = _with_v35_multirate_and_adjoint_identity(_summary())
    summary[
        "multirate_electromechanical_event_interpolation_work_power_timegrid_frame_mesh_result_generation_identity"
    ].update(
        {
            "event_generation": "multirate-coupling-220",
            "timegrid_generation": "multirate-coupling-219",
            "result_generation": "multirate-coupling-218",
            "result_electrical_time_s": [0.0, 0.0004, 0.001],
            "result_mechanical_time_s": [0.0, 0.0005, 0.001],
            "result_event_time_s": 0.0006,
            "result_event_interpolation_side": "left_continuous_before_event",
            "result_substep_owner": "mechanical_only",
            "result_coordinate_frame": "rotor_cylindrical",
            "result_electrical_input_energy_j": 0.010,
            "result_mechanical_output_work_j": 0.011,
            "result_dissipated_energy_j": -0.001,
            "result_energy_balance_tolerance_j": 1.0e-15,
            "result_coupling_mesh_sha256": "a" * 64,
            "accepted_coupling_result_sha256": "b" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "multirate_electromechanical_results_use_current_event_timegrids_work_power_frame_mesh_and_result"
    ]


def test_v35_public_adjoint_sensitivity_objective_chainrule_constraint_fd_mesh_owner_mismatch() -> None:
    summary = _with_v35_multirate_and_adjoint_identity(_summary())
    summary[
        "adjoint_objective_design_chainrule_constraint_fd_mesh_solution_gradient_generation_identity"
    ].update(
        {
            "objective_generation": "adjoint-sensitivity-220",
            "fd_generation": "adjoint-sensitivity-219",
            "result_generation": "adjoint-sensitivity-218",
            "result_objective_tag": "mean_torque",
            "result_design_variable": "magnet_arc_deg",
            "result_design_scale": 10.0,
            "result_active_constraint": "none",
            "chainrule_gradient": -2.5,
            "finite_difference_gradient": 0.25,
            "result_gradient_tolerance": 1.0e-12,
            "result_fd_perturbation": 0.1,
            "result_sensitivity_mesh_sha256": "c" * 64,
            "result_primal_solution_sha256": "d" * 64,
            "accepted_gradient_result_sha256": "e" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "adjoint_sensitivities_use_current_objective_design_chainrule_constraint_fd_mesh_solution_and_result"
    ]


def test_v35_public_rejects_self_consistent_multirate_energy_creation() -> None:
    summary = _with_v35_multirate_and_adjoint_identity(_summary())
    identity = summary[
        "multirate_electromechanical_event_interpolation_work_power_timegrid_frame_mesh_result_generation_identity"
    ]
    identity["mechanical_output_work_j"] = 0.013
    identity["result_mechanical_output_work_j"] = 0.013
    assert gate(summary)["status"] == "needs_attention"


def test_v35_public_rejects_self_consistent_adjoint_fd_disagreement() -> None:
    summary = _with_v35_multirate_and_adjoint_identity(_summary())
    identity = summary[
        "adjoint_objective_design_chainrule_constraint_fd_mesh_solution_gradient_generation_identity"
    ]
    identity["finite_difference_gradient"] = 3.0
    assert gate(summary)["status"] == "needs_attention"
