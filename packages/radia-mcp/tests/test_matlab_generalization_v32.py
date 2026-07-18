from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v31 import _summary_v31


_PROMOTED_CASE_IDS = (
    "v32_public_adaptive_cq_timestep_contour_rebuild_interpolation_history_error_restart_mismatch",
    "v32_public_fembem_modal_transient_mass_damping_initial_condition_energy_balance_mismatch",
)


def _summary_v32():
    payload = deepcopy(_summary_v31())
    generation = "adaptive-cq-361"
    payload[
        "adaptive_cq_timestep_contour_history_interpolation_error_restart_operator_mesh_result_identity"
    ] = {
        "cq_generation": generation,
        **{
            key: generation
            for key in (
                "timestep_cq_generation",
                "contour_cq_generation",
                "history_cq_generation",
                "error_cq_generation",
                "restart_cq_generation",
                "operator_cq_generation",
                "mesh_cq_generation",
                "result_cq_generation",
            )
        },
        "timestep_schedule_s": [1.0e-5, 1.0e-5, 5.0e-6, 5.0e-6],
        "result_timestep_schedule_s": [1.0e-5, 1.0e-5, 5.0e-6, 5.0e-6],
        "contour_rebuild_indices": [0, 2],
        "result_contour_rebuild_indices": [0, 2],
        "history_interpolation": "barycentric_causal",
        "result_history_interpolation": "barycentric_causal",
        "local_error_estimates": [1.0e-4, 8.0e-5, 2.0e-5, 1.0e-5],
        "result_local_error_estimates": [1.0e-4, 8.0e-5, 2.0e-5, 1.0e-5],
        "local_error_tolerance": 1.0e-3,
        "result_local_error_tolerance": 1.0e-3,
        "restart_index": 2,
        "result_restart_index": 2,
        "operator_owner_sha256": "1" * 64,
        "result_operator_owner_sha256": "1" * 64,
        "history_sha256": "2" * 64,
        "loaded_history_sha256": "2" * 64,
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    generation = "modal-fembem-361"
    payload[
        "fembem_modal_transient_mass_damping_initial_projection_truncation_energy_mesh_history_result_identity"
    ] = {
        "modal_generation": generation,
        **{
            key: generation
            for key in (
                "mass_modal_generation",
                "damping_modal_generation",
                "initial_modal_generation",
                "truncation_modal_generation",
                "energy_modal_generation",
                "mesh_modal_generation",
                "history_modal_generation",
                "result_modal_generation",
            )
        },
        "mass_normalization": "M_orthonormal",
        "result_mass_normalization": "M_orthonormal",
        "damping_model": "rayleigh",
        "result_damping_model": "rayleigh",
        "rayleigh_coefficients": [0.01, 1.0e-5],
        "result_rayleigh_coefficients": [0.01, 1.0e-5],
        "initial_displacement_projection": [1.0, 0.2, 0.0],
        "result_initial_displacement_projection": [1.0, 0.2, 0.0],
        "initial_velocity_projection": [0.0, 0.0, 0.0],
        "result_initial_velocity_projection": [0.0, 0.0, 0.0],
        "modal_count": 3,
        "result_modal_count": 3,
        "truncation_frequency_hz": 1200.0,
        "result_truncation_frequency_hz": 1200.0,
        "initial_energy_j": 0.5,
        "result_initial_energy_j": 0.5,
        "radiated_energy_j": 0.2,
        "result_radiated_energy_j": 0.2,
        "dissipated_energy_j": 0.1,
        "result_dissipated_energy_j": 0.1,
        "final_energy_j": 0.2,
        "result_final_energy_j": 0.2,
        "mesh_sha256": "5" * 64,
        "result_mesh_sha256": "5" * 64,
        "time_history_owner": "fembem/case-361/modal-transient",
        "accepted_time_history_owner": "fembem/case-361/modal-transient",
        "result_sha256": "6" * 64,
        "accepted_result_sha256": "6" * 64,
    }
    return payload


def test_v32_public_positive_adaptive_cq_and_modal_fembem_transient() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v32())["status"] == "ok"


def test_v32_public_adaptive_cq_timestep_contour_rebuild_interpolation_history_error_restart_mismatch() -> None:
    payload = _summary_v32()
    identity = payload[
        "adaptive_cq_timestep_contour_history_interpolation_error_restart_operator_mesh_result_identity"
    ]
    identity.update(
        {
            "timestep_cq_generation": "adaptive-cq-360",
            "operator_cq_generation": "adaptive-cq-359",
            "result_cq_generation": "adaptive-cq-358",
            "result_timestep_schedule_s": [1.0e-5, 5.0e-6, 1.0e-5],
            "result_contour_rebuild_indices": [0],
            "result_history_interpolation": "linear_noncausal",
            "result_local_error_estimates": [0.2, 0.1, 0.05],
            "result_local_error_tolerance": 1.0e-4,
            "result_restart_index": 1,
            "result_operator_owner_sha256": "b" * 64,
            "loaded_history_sha256": "c" * 64,
            "result_mesh_sha256": "d" * 64,
            "accepted_result_sha256": "e" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "adaptive_cq_uses_current_timesteps_contour_rebuild_history_interpolation_error_restart_operator_mesh_and_result"
    ]


def test_v32_public_fembem_modal_transient_mass_damping_initial_condition_energy_balance_mismatch() -> None:
    payload = _summary_v32()
    identity = payload[
        "fembem_modal_transient_mass_damping_initial_projection_truncation_energy_mesh_history_result_identity"
    ]
    identity.update(
        {
            "mass_modal_generation": "modal-fembem-360",
            "history_modal_generation": "modal-fembem-359",
            "result_modal_generation": "modal-fembem-358",
            "result_mass_normalization": "euclidean",
            "result_damping_model": "none",
            "result_rayleigh_coefficients": [0.0, 0.0],
            "result_initial_displacement_projection": [0.0, 1.0],
            "result_initial_velocity_projection": [1.0, 0.0],
            "result_modal_count": 2,
            "result_truncation_frequency_hz": 600.0,
            "result_initial_energy_j": 0.4,
            "result_radiated_energy_j": 0.6,
            "result_dissipated_energy_j": -0.1,
            "result_final_energy_j": 0.3,
            "result_mesh_sha256": "f" * 64,
            "accepted_time_history_owner": "fembem/old-history",
            "accepted_result_sha256": "0" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "modal_fembem_transient_uses_current_mass_damping_initial_projection_truncation_energy_mesh_history_and_result"
    ]
