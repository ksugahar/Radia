from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v32 import _summary_v32


_PROMOTED_CASE_IDS = (
    "v33_public_calderon_projector_v_k_kt_w_mass_duality_normal_quadrature_mismatch",
    "v33_public_cq_symbol_contour_conjugate_symmetry_causal_ifft_parseval_passivity_mismatch",
)


def _summary_v33():
    payload = deepcopy(_summary_v32())
    generation = "calderon-p1-381"
    payload[
        "calderon_projector_p1_v_k_kt_w_mass_duality_normal_quadrature_mesh_owner_result_identity"
    ] = {
        "calderon_generation": generation,
        **{
            key: generation
            for key in (
                "space_calderon_generation",
                "operator_calderon_generation",
                "mass_calderon_generation",
                "normal_calderon_generation",
                "quadrature_calderon_generation",
                "mesh_calderon_generation",
                "projector_calderon_generation",
                "owner_calderon_generation",
                "result_calderon_generation",
            )
        },
        "trial_space": "P1",
        "result_trial_space": "P1",
        "test_space": "P1",
        "result_test_space": "P1",
        "projector_convention": "interior_calderon_outward",
        "result_projector_convention": "interior_calderon_outward",
        "block_order": ["dirichlet", "neumann"],
        "result_block_order": ["dirichlet", "neumann"],
        "v_sign": -1,
        "result_v_sign": -1,
        "k_sign": 1,
        "result_k_sign": 1,
        "kt_sign": -1,
        "result_kt_sign": -1,
        "w_sign": -1,
        "result_w_sign": -1,
        "mass_duality_residual": 2.0e-11,
        "result_mass_duality_residual": 2.0e-11,
        "mass_duality_tolerance": 1.0e-8,
        "result_mass_duality_tolerance": 1.0e-8,
        "normal_orientation": "outward",
        "result_normal_orientation": "outward",
        "singular_quadrature": "duffy_principal_value_p1",
        "result_singular_quadrature": "duffy_principal_value_p1",
        "projector_residual": 5.0e-10,
        "result_projector_residual": 5.0e-10,
        "projector_tolerance": 1.0e-8,
        "result_projector_tolerance": 1.0e-8,
        "operator_sha256": "1" * 64,
        "result_operator_sha256": "1" * 64,
        "mass_sha256": "2" * 64,
        "result_mass_sha256": "2" * 64,
        "boundary_mesh_sha256": "3" * 64,
        "result_boundary_mesh_sha256": "3" * 64,
        "result_owner": "fembem/calderon/case-381",
        "accepted_result_owner": "fembem/calderon/case-381",
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    generation = "cq-physical-381"
    transfer_samples = [
        [1.0, 0.0],
        [0.8, 0.2],
        [0.5, 0.1],
        [0.5, -0.1],
        [0.8, -0.2],
    ]
    payload[
        "cq_symbol_contour_transfer_conjugate_causal_ifft_parseval_passivity_timestep_operator_result_identity"
    ] = {
        "cq_generation": generation,
        **{
            key: generation
            for key in (
                "symbol_cq_generation",
                "contour_cq_generation",
                "transfer_cq_generation",
                "symmetry_cq_generation",
                "causality_cq_generation",
                "parseval_cq_generation",
                "passivity_cq_generation",
                "timestep_cq_generation",
                "operator_cq_generation",
                "result_cq_generation",
            )
        },
        "multistep_symbol": "BDF2",
        "result_multistep_symbol": "BDF2",
        "symbol_coefficients": [1.5, -2.0, 0.5],
        "result_symbol_coefficients": [1.5, -2.0, 0.5],
        "contour_radius": 0.92,
        "result_contour_radius": 0.92,
        "transfer_samples_ri": transfer_samples,
        "result_transfer_samples_ri": [list(row) for row in transfer_samples],
        "time_response": [0.0, 0.2, 0.1, 0.05],
        "result_time_response": [0.0, 0.2, 0.1, 0.05],
        "negative_time_energy": 0.0,
        "result_negative_time_energy": 0.0,
        "time_domain_work": 0.25,
        "result_time_domain_work": 0.25,
        "frequency_domain_work": 0.25,
        "result_frequency_domain_work": 0.25,
        "parseval_tolerance": 1.0e-10,
        "result_parseval_tolerance": 1.0e-10,
        "minimum_real_transfer": 0.5,
        "result_minimum_real_transfer": 0.5,
        "passivity_sign": "nonnegative_real_transfer",
        "result_passivity_sign": "nonnegative_real_transfer",
        "timestep_s": 1.0e-5,
        "result_timestep_s": 1.0e-5,
        "operator_family": "p1_calderon_bem",
        "result_operator_family": "p1_calderon_bem",
        "operator_sha256": "5" * 64,
        "result_operator_sha256": "5" * 64,
        "result_owner": "fembem/cq/case-381",
        "accepted_result_owner": "fembem/cq/case-381",
        "result_sha256": "6" * 64,
        "accepted_result_sha256": "6" * 64,
    }
    return payload


def test_v33_public_positive_calderon_projector_and_cq_closure() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v33())["status"] == "ok"


def test_v33_public_calderon_projector_v_k_kt_w_mass_duality_normal_quadrature_mismatch() -> None:
    payload = _summary_v33()
    identity = payload[
        "calderon_projector_p1_v_k_kt_w_mass_duality_normal_quadrature_mesh_owner_result_identity"
    ]
    identity.update(
        {
            "operator_calderon_generation": "calderon-p1-380",
            "mesh_calderon_generation": "calderon-p1-379",
            "result_calderon_generation": "calderon-p1-378",
            "result_trial_space": "P0",
            "result_test_space": "P0",
            "result_projector_convention": "exterior_stale_normal",
            "result_block_order": ["neumann", "dirichlet"],
            "result_v_sign": 1,
            "result_k_sign": -1,
            "result_kt_sign": 1,
            "result_w_sign": 1,
            "result_mass_duality_residual": 0.2,
            "result_mass_duality_tolerance": 1.0e-4,
            "result_normal_orientation": "inward",
            "result_singular_quadrature": "centroid_regular",
            "result_projector_residual": 0.4,
            "result_projector_tolerance": 1.0e-4,
            "result_operator_sha256": "c" * 64,
            "result_mass_sha256": "d" * 64,
            "result_boundary_mesh_sha256": "e" * 64,
            "accepted_result_owner": "fembem/calderon/old",
            "accepted_result_sha256": "f" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "calderon_projector_uses_current_p1_spaces_v_k_kt_w_mass_duality_normals_quadrature_mesh_owner_and_result"
    ]


def test_v33_public_cq_symbol_contour_conjugate_symmetry_causal_ifft_parseval_passivity_mismatch() -> None:
    payload = _summary_v33()
    identity = payload[
        "cq_symbol_contour_transfer_conjugate_causal_ifft_parseval_passivity_timestep_operator_result_identity"
    ]
    identity.update(
        {
            "symbol_cq_generation": "cq-physical-380",
            "transfer_cq_generation": "cq-physical-379",
            "result_cq_generation": "cq-physical-378",
            "result_multistep_symbol": "BDF1",
            "result_symbol_coefficients": [1.0, -1.0],
            "result_contour_radius": 1.1,
            "result_transfer_samples_ri": [[1.0, 0.1], [0.8, 0.3]],
            "result_time_response": [0.4, -0.2],
            "result_negative_time_energy": 0.1,
            "result_time_domain_work": 0.1,
            "result_frequency_domain_work": 0.7,
            "result_parseval_tolerance": 1.0e-4,
            "result_minimum_real_transfer": -0.3,
            "result_passivity_sign": "negative_real_transfer",
            "result_timestep_s": -1.0e-5,
            "result_operator_family": "stale_operator",
            "result_operator_sha256": "0" * 64,
            "accepted_result_owner": "fembem/cq/old",
            "accepted_result_sha256": "1" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_uses_current_symbol_contour_conjugate_transfer_causal_ifft_parseval_passivity_timestep_operator_and_result"
    ]
