from __future__ import annotations

import math
from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v37 import _summary_v37


_PROMOTED_CASE_IDS = (
    "v38_public_simp_topology_filter_projection_volume_compliance_adjoint_fd_kkt_mismatch",
    "v38_public_fembem_model_reduction_projection_stability_passivity_moment_error_owner_mismatch",
)


def _summary_v38():
    payload = deepcopy(_summary_v37())
    generation = "simp-topology-614"
    density = [0.4, 0.5, 0.6, 0.5]
    filter_matrix = [
        [0.75, 0.25, 0.0, 0.0],
        [0.25, 0.5, 0.25, 0.0],
        [0.0, 0.25, 0.5, 0.25],
        [0.0, 0.0, 0.25, 0.75],
    ]
    filtered = [
        sum(weight * value for weight, value in zip(row, density))
        for row in filter_matrix
    ]
    beta = [1.0, 2.0, 4.0, 8.0]
    eta = 0.5
    denominator = math.tanh(beta[-1] * eta) + math.tanh(beta[-1] * (1.0 - eta))
    projected = [
        (math.tanh(beta[-1] * eta) + math.tanh(beta[-1] * (value - eta)))
        / denominator
        for value in filtered
    ]
    mirrored = {
        "design_density": density,
        "density_filter_matrix": filter_matrix,
        "filtered_density": filtered,
        "projection_beta_continuation": beta,
        "projection_eta": eta,
        "projected_density": projected,
        "volume_fraction": sum(projected) / len(projected),
        "volume_fraction_limit": 0.5,
        "compliance": 12.0,
        "adjoint_compliance_gradient": [-0.25] * 4,
        "finite_difference_compliance_gradient": [-0.2500001] * 4,
        "volume_gradient": [0.25] * 4,
        "volume_lagrange_multiplier": 1.0,
        "kkt_stationarity_residual": 0.0,
        "gradient_relative_tolerance": 1.0e-5,
        "topology_mesh_sha256": "1" * 64,
    }
    payload["simp_topology_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_result_identity"] = {
        "topology_generation": generation,
        **{key: generation for key in (
            "density_generation", "filter_generation", "projection_generation",
            "volume_generation", "compliance_generation", "adjoint_generation",
            "fd_generation", "kkt_generation", "mesh_generation",
            "owner_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "topology_owner": "optimization/topology-614",
        "accepted_topology_owner": "optimization/topology-614",
        "topology_result_sha256": "2" * 64,
        "accepted_topology_result_sha256": "2" * 64,
    }

    generation = "fembem-reduction-614"
    mirrored = {
        "full_order": 4,
        "reduced_order": 2,
        "trial_projection_basis": [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
        "test_projection_basis": [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
        "biorthogonality_gram": [[1.0, 0.0], [0.0, 1.0]],
        "full_model_poles": [[-1.0, 0.0], [-2.0, 0.0], [-3.0, 0.0], [-4.0, 0.0]],
        "reduced_model_poles": [[-1.1, 0.0], [-2.1, 0.0]],
        "minimum_passivity_eigenvalue": 0.1,
        "matched_moments_full": [[1.0, 0.0], [0.5, 0.0]],
        "matched_moments_reduced": [[1.0, 0.0], [0.5, 0.0]],
        "frequency_hz": [100.0, 200.0, 400.0],
        "frequency_response_relative_error": [0.01, 0.005, 0.002],
        "maximum_frequency_response_relative_error": 0.02,
        "reduction_mesh_sha256": "3" * 64,
    }
    payload["fembem_model_reduction_projection_order_stability_passivity_moment_frequency_error_full_mesh_owner_result_identity"] = {
        "reduction_generation": generation,
        **{key: generation for key in (
            "projection_generation", "order_generation", "stability_generation",
            "passivity_generation", "moment_generation", "frequency_generation",
            "error_generation", "full_model_generation", "mesh_generation",
            "owner_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "full_model_owner": "fembem/full-614",
        "accepted_full_model_owner": "fembem/full-614",
        "reduction_owner": "fembem/reduced-614",
        "accepted_reduction_owner": "fembem/reduced-614",
        "reduction_result_sha256": "4" * 64,
        "accepted_reduction_result_sha256": "4" * 64,
    }
    return payload


def test_v38_public_positive_simp_and_fembem_model_reduction_closure():
    assert regularized_trace_inverse_path_gate(_summary_v38())["status"] == "ok"


def test_v38_public_simp_topology_filter_projection_volume_compliance_adjoint_fd_kkt_mismatch():
    payload = _summary_v38()
    identity = payload["simp_topology_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_result_identity"]
    identity.update({
        "filter_generation": "simp-topology-613",
        "kkt_generation": "simp-topology-612",
        "result_generation": "simp-topology-611",
        "result_density_filter_matrix": [[2.0, -1.0]],
        "result_filtered_density": [2.0],
        "result_projection_beta_continuation": [8.0, 1.0],
        "result_projected_density": [-1.0, 2.0],
        "result_volume_fraction": 2.0,
        "result_compliance": -12.0,
        "result_adjoint_compliance_gradient": [9.0],
        "result_finite_difference_compliance_gradient": [-9.0],
        "result_kkt_stationarity_residual": 9.0,
        "accepted_topology_owner": "optimization/old",
        "accepted_topology_result_sha256": "a" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["simp_topology_uses_current_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_and_result"]


def test_v38_public_fembem_model_reduction_projection_stability_passivity_moment_error_owner_mismatch():
    payload = _summary_v38()
    identity = payload["fembem_model_reduction_projection_order_stability_passivity_moment_frequency_error_full_mesh_owner_result_identity"]
    identity.update({
        "projection_generation": "fembem-reduction-613",
        "passivity_generation": "fembem-reduction-612",
        "result_generation": "fembem-reduction-611",
        "result_reduced_order": 5,
        "result_trial_projection_basis": [[1.0]],
        "result_test_projection_basis": [[-1.0]],
        "result_biorthogonality_gram": [[-1.0]],
        "result_reduced_model_poles": [[1.0, 0.0]],
        "result_minimum_passivity_eigenvalue": -1.0,
        "result_matched_moments_reduced": [[9.0, 9.0]],
        "result_frequency_response_relative_error": [1.0],
        "accepted_full_model_owner": "fembem/old",
        "accepted_reduction_owner": "fembem/old",
        "accepted_reduction_result_sha256": "b" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fembem_model_reduction_uses_current_projection_order_stability_passivity_moments_frequency_error_full_model_mesh_owners_and_result"]


def test_v38_public_rejects_self_consistent_nonstochastic_density_filter():
    payload = _summary_v38()
    identity = payload["simp_topology_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_result_identity"]
    identity["density_filter_matrix"] = [[1.5, -0.5, 0.0, 0.0]] * 4
    identity["result_density_filter_matrix"] = deepcopy(identity["density_filter_matrix"])
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_unstable_reduced_pole():
    payload = _summary_v38()
    identity = payload["fembem_model_reduction_projection_order_stability_passivity_moment_frequency_error_full_mesh_owner_result_identity"]
    identity["reduced_model_poles"] = [[1.0, 0.0], [-2.1, 0.0]]
    identity["result_reduced_model_poles"] = deepcopy(identity["reduced_model_poles"])
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"
