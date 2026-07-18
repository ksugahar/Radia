from __future__ import annotations

from test_comsol_generalization_v32 import (
    _summary,
    _with_v32_nonlinear_and_eigenmode_identity,
    gate,
)


_PROMOTED_CASE_IDS = (
    "v33_public_contact_complementarity_gap_pressure_active_set_friction_dissipation_mismatch",
    "v33_public_field_circuit_dae_charge_current_event_energy_constraint_mismatch",
)


def _with_v33_contact_and_dae_identity(summary: dict) -> dict:
    summary = _with_v32_nonlinear_and_eigenmode_identity(summary)
    generation = "contact-complementarity-201"
    summary[
        "contact_gap_pressure_active_set_friction_dissipation_normal_mesh_result_generation_identity"
    ] = {
        "contact_generation": generation,
        **{
            key: generation
            for key in (
                "gap_generation",
                "pressure_generation",
                "active_set_generation",
                "friction_generation",
                "dissipation_generation",
                "normal_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "contact_pair": "contact1",
        "result_contact_pair": "contact1",
        "active_contact_ids": [3, 7],
        "result_active_contact_ids": [3, 7],
        "normal_gap_m": [0.0, 0.0],
        "result_normal_gap_m": [0.0, 0.0],
        "normal_pressure_pa": [1.0e6, 2.0e6],
        "result_normal_pressure_pa": [1.0e6, 2.0e6],
        "tangential_slip_m": [0.0, 1.0e-4],
        "result_tangential_slip_m": [0.0, 1.0e-4],
        "friction_traction_pa": [0.0, 2.0e5],
        "result_friction_traction_pa": [0.0, 2.0e5],
        "contact_area_m2": [1.0e-4, 1.0e-4],
        "result_contact_area_m2": [1.0e-4, 1.0e-4],
        "friction_coefficient": 0.2,
        "result_friction_coefficient": 0.2,
        "friction_dissipation_j": 0.002,
        "result_friction_dissipation_j": 0.002,
        "normal_orientation": "outward_slave_to_master",
        "result_normal_orientation": "outward_slave_to_master",
        "contact_mesh_sha256": "1" * 64,
        "result_contact_mesh_sha256": "1" * 64,
        "contact_result_sha256": "2" * 64,
        "accepted_contact_result_sha256": "2" * 64,
    }
    generation = "field-circuit-dae-201"
    summary[
        "field_circuit_dae_charge_current_event_energy_time_dataset_result_generation_identity"
    ] = {
        "dae_generation": generation,
        **{
            key: generation
            for key in (
                "charge_generation",
                "current_generation",
                "event_generation",
                "energy_generation",
                "time_generation",
                "dataset_generation",
                "result_generation",
            )
        },
        "time_s": [0.0, 0.5e-3, 1.0e-3],
        "result_time_s": [0.0, 0.5e-3, 1.0e-3],
        "switch_event_time_s": 0.5e-3,
        "result_switch_event_time_s": 0.5e-3,
        "event_side": "right_limit_after_event",
        "result_event_side": "right_limit_after_event",
        "charge_c": [0.0, 1.0e-6, 1.5e-6],
        "result_charge_c": [0.0, 1.0e-6, 1.5e-6],
        "integrated_current_c": [0.0, 1.0e-6, 1.5e-6],
        "result_integrated_current_c": [0.0, 1.0e-6, 1.5e-6],
        "algebraic_residual_c": [0.0, 1.0e-14, 0.0],
        "accepted_algebraic_residual_c": [0.0, 1.0e-14, 0.0],
        "algebraic_tolerance_c": 1.0e-12,
        "current_sign_convention": "positive_into_field_device",
        "result_current_sign_convention": "positive_into_field_device",
        "stored_energy_before_j": 0.002,
        "result_stored_energy_before_j": 0.002,
        "stored_energy_after_j": 0.0018,
        "result_stored_energy_after_j": 0.0018,
        "switch_dissipation_j": 0.0002,
        "result_switch_dissipation_j": 0.0002,
        "dataset_owner": "dset1/sol2",
        "result_dataset_owner": "dset1/sol2",
        "dae_dataset_sha256": "3" * 64,
        "result_dae_dataset_sha256": "3" * 64,
        "dae_result_sha256": "4" * 64,
        "accepted_dae_result_sha256": "4" * 64,
    }
    return summary


def test_v33_public_positive_contact_and_dae_contracts() -> None:
    assert gate(_with_v33_contact_and_dae_identity(_summary()))["status"] == "ok"


def test_v33_public_contact_complementarity_gap_pressure_active_set_friction_dissipation_mismatch() -> None:
    summary = _with_v33_contact_and_dae_identity(_summary())
    summary[
        "contact_gap_pressure_active_set_friction_dissipation_normal_mesh_result_generation_identity"
    ].update(
        {
            "gap_generation": "contact-complementarity-200",
            "active_set_generation": "contact-complementarity-199",
            "result_generation": "contact-complementarity-198",
            "result_contact_pair": "contact_old",
            "result_active_contact_ids": [7, 9],
            "result_normal_gap_m": [1.0e-3, -2.0e-4],
            "result_normal_pressure_pa": [-1.0e6, 2.0e6],
            "result_tangential_slip_m": [1.0e-3, 0.0],
            "result_friction_traction_pa": [4.0e5, 0.0],
            "result_contact_area_m2": [2.0e-4, 1.0e-4],
            "result_friction_coefficient": 0.1,
            "result_friction_dissipation_j": -0.01,
            "result_normal_orientation": "outward_master_to_slave",
            "result_contact_mesh_sha256": "9" * 64,
            "accepted_contact_result_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "contact_results_satisfy_current_complementarity_active_set_friction_dissipation_normal_mesh_and_result"
    ]


def test_v33_public_field_circuit_dae_charge_current_event_energy_constraint_mismatch() -> None:
    summary = _with_v33_contact_and_dae_identity(_summary())
    summary[
        "field_circuit_dae_charge_current_event_energy_time_dataset_result_generation_identity"
    ].update(
        {
            "charge_generation": "field-circuit-dae-200",
            "event_generation": "field-circuit-dae-199",
            "result_generation": "field-circuit-dae-198",
            "result_time_s": [0.0, 0.6e-3, 1.0e-3],
            "result_switch_event_time_s": 0.7e-3,
            "result_event_side": "left_limit_before_event",
            "result_charge_c": [0.0, -1.0e-6, 0.5e-6],
            "result_integrated_current_c": [0.0, 0.4e-6, 1.5e-6],
            "accepted_algebraic_residual_c": [0.0, 1.0e-5, 0.0],
            "result_current_sign_convention": "positive_out_of_field_device",
            "result_stored_energy_after_j": 0.0022,
            "result_switch_dissipation_j": -0.0002,
            "result_dataset_owner": "dset_old/sol1",
            "result_dae_dataset_sha256": "b" * 64,
            "accepted_dae_result_sha256": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "field_circuit_dae_results_use_current_charge_current_event_energy_time_dataset_and_result"
    ]
