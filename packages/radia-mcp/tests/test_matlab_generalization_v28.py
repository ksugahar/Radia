from __future__ import annotations

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v27 import _summary_v27


_PROMOTED_CASE_IDS = (
    "v28_public_hmatrix_aca_cluster_permutation_admissibility_rank_tolerance_kernel_mesh_result_mismatch",
    "v28_public_calderon_cq_operator_v_k_trace_normal_frequency_grid_inverse_transform_mismatch",
)


def _summary_v28():
    summary = _summary_v27()
    generation = "hmatrix-aca-321"
    summary[
        "hmatrix_aca_cluster_permutation_admissibility_rank_tolerance_kernel_mesh_result_generation_identity"
    ] = {
        "hmatrix_generation": generation,
        "cluster_hmatrix_generation": generation,
        "permutation_hmatrix_generation": generation,
        "admissibility_hmatrix_generation": generation,
        "rank_hmatrix_generation": generation,
        "tolerance_hmatrix_generation": generation,
        "kernel_hmatrix_generation": generation,
        "mesh_hmatrix_generation": generation,
        "result_hmatrix_generation": generation,
        "cluster_permutation": [3, 1, 4, 2],
        "result_cluster_permutation": [3, 1, 4, 2],
        "admissibility_rule": "eta-weak",
        "result_admissibility_rule": "eta-weak",
        "admissibility_eta": 2.0,
        "result_admissibility_eta": 2.0,
        "aca_rank": 8,
        "result_aca_rank": 8,
        "relative_tolerance": 1.0e-6,
        "result_relative_tolerance": 1.0e-6,
        "kernel": "helmholtz-single-layer-p1",
        "result_kernel": "helmholtz-single-layer-p1",
        "cluster_tree_sha256": "1" * 64,
        "loaded_cluster_tree_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "result_sha256": "3" * 64,
        "accepted_result_sha256": "3" * 64,
    }
    generation = "calderon-cq-321"
    summary[
        "calderon_cq_operator_v_k_trace_normal_frequency_grid_inverse_transform_mesh_result_generation_identity"
    ] = {
        "calderon_generation": generation,
        "v_calderon_generation": generation,
        "k_calderon_generation": generation,
        "trace_calderon_generation": generation,
        "normal_calderon_generation": generation,
        "frequency_calderon_generation": generation,
        "inverse_calderon_generation": generation,
        "mesh_calderon_generation": generation,
        "result_calderon_generation": generation,
        "v_operator_sha256": "4" * 64,
        "result_v_operator_sha256": "4" * 64,
        "k_operator_sha256": "5" * 64,
        "result_k_operator_sha256": "5" * 64,
        "trace_basis": "p1-nodal-boundary-trace",
        "result_trace_basis": "p1-nodal-boundary-trace",
        "trace_shape": [48, 120],
        "result_trace_shape": [48, 120],
        "boundary_normal": "volume-outward",
        "result_boundary_normal": "volume-outward",
        "laplace_frequency_ri": [[10.0, 0.0], [10.0, 20.0], [10.0, 40.0]],
        "result_laplace_frequency_ri": [[10.0, 0.0], [10.0, 20.0], [10.0, 40.0]],
        "inverse_transform": "bdf2-cq-ifft-real",
        "result_inverse_transform": "bdf2-cq-ifft-real",
        "boundary_mesh_sha256": "6" * 64,
        "result_boundary_mesh_sha256": "6" * 64,
        "result_sha256": "7" * 64,
        "accepted_result_sha256": "7" * 64,
    }
    return summary


def test_v28_public_positive_hmatrix_and_calderon_cq_identities() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v28())["status"] == "ok"


def test_v28_public_rejects_hmatrix_aca_identity_mismatch() -> None:
    summary = _summary_v28()
    identity = summary[
        "hmatrix_aca_cluster_permutation_admissibility_rank_tolerance_kernel_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "permutation_hmatrix_generation": "hmatrix-aca-320",
            "mesh_hmatrix_generation": "hmatrix-aca-319",
            "result_cluster_permutation": [1, 2, 3, 4],
            "result_admissibility_rule": "strong",
            "result_admissibility_eta": 0.5,
            "result_aca_rank": 3,
            "result_relative_tolerance": 1.0e-2,
            "result_kernel": "laplace-p0",
            "loaded_cluster_tree_sha256": "d" * 64,
            "result_mesh_sha256": "e" * 64,
            "accepted_result_sha256": "f" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hmatrix_aca_uses_current_clusters_permutation_admissibility_rank_tolerance_kernel_mesh_and_result"
    ]


def test_v28_public_rejects_calderon_cq_identity_mismatch() -> None:
    summary = _summary_v28()
    identity = summary[
        "calderon_cq_operator_v_k_trace_normal_frequency_grid_inverse_transform_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "v_calderon_generation": "calderon-cq-320",
            "normal_calderon_generation": "calderon-cq-319",
            "result_v_operator_sha256": "0" * 64,
            "result_k_operator_sha256": "1" * 64,
            "result_trace_basis": "p0-cell",
            "result_trace_shape": [47, 120],
            "result_boundary_normal": "volume-inward",
            "result_laplace_frequency_ri": [[0.0, 10.0]],
            "result_inverse_transform": "direct-real-ifft",
            "result_boundary_mesh_sha256": "2" * 64,
            "accepted_result_sha256": "3" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "calderon_cq_uses_current_v_k_trace_normals_frequency_inverse_mesh_and_result"
    ]
