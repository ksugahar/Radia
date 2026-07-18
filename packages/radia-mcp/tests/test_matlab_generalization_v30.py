from __future__ import annotations

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v29 import _summary_v29


_PROMOTED_CASE_IDS = (
    "v30_public_hmatrix_recompression_svd_tolerance_norm_rank_permutation_operator_error_mismatch",
    "v30_public_cq_restart_block_history_startup_weights_time_index_sample_count_digest_mismatch",
)


def _summary_v30():
    summary = _summary_v29()
    generation = "hmatrix-recompress-341"
    summary[
        "hmatrix_recompression_svd_tolerance_norm_rank_permutation_operator_mesh_result_identity"
    ] = {
        "hmatrix_generation": generation,
        "svd_hmatrix_generation": generation,
        "tolerance_hmatrix_generation": generation,
        "rank_hmatrix_generation": generation,
        "permutation_hmatrix_generation": generation,
        "operator_hmatrix_generation": generation,
        "mesh_hmatrix_generation": generation,
        "result_hmatrix_generation": generation,
        "svd_basis": "euclidean-orthonormal",
        "result_svd_basis": "euclidean-orthonormal",
        "tolerance": 1.0e-6,
        "result_tolerance": 1.0e-6,
        "tolerance_norm": "spectral-relative",
        "result_tolerance_norm": "spectral-relative",
        "block_ranks_before": [12, 10, 8],
        "result_block_ranks_before": [12, 10, 8],
        "block_ranks_after": [6, 5, 4],
        "result_block_ranks_after": [6, 5, 4],
        "row_permutation": [2, 0, 1],
        "result_row_permutation": [2, 0, 1],
        "column_permutation": [1, 2, 0],
        "result_column_permutation": [1, 2, 0],
        "operator_relative_error": 5.0e-7,
        "result_operator_relative_error": 5.0e-7,
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "cq-restart-341"
    summary[
        "cq_restart_block_history_startup_weight_time_index_sample_contour_operator_result_identity"
    ] = {
        "cq_generation": generation,
        "block_cq_generation": generation,
        "history_cq_generation": generation,
        "startup_cq_generation": generation,
        "time_cq_generation": generation,
        "sample_cq_generation": generation,
        "owner_cq_generation": generation,
        "result_cq_generation": generation,
        "block_size": 16,
        "result_block_size": 16,
        "completed_block_ids": [0, 1, 2],
        "result_completed_block_ids": [0, 1, 2],
        "history_sample_count": 48,
        "result_history_sample_count": 48,
        "startup_weights_ri": [[1.0, 0.0], [0.5, -0.1]],
        "result_startup_weights_ri": [[1.0, 0.0], [0.5, -0.1]],
        "restart_time_index": 48,
        "result_restart_time_index": 48,
        "total_sample_count": 128,
        "result_total_sample_count": 128,
        "contour_owner_sha256": "3" * 64,
        "result_contour_owner_sha256": "3" * 64,
        "operator_owner_sha256": "4" * 64,
        "result_operator_owner_sha256": "4" * 64,
        "history_sha256": "5" * 64,
        "loaded_history_sha256": "5" * 64,
        "result_sha256": "6" * 64,
        "accepted_result_sha256": "6" * 64,
    }
    return summary


def test_v30_public_positive_hmatrix_recompression_and_cq_restart() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v30())["status"] == "ok"


def test_v30_public_hmatrix_recompression_svd_tolerance_norm_rank_permutation_operator_error_mismatch() -> None:
    summary = _summary_v30()
    identity = summary[
        "hmatrix_recompression_svd_tolerance_norm_rank_permutation_operator_mesh_result_identity"
    ]
    identity.update(
        {
            "svd_hmatrix_generation": "hmatrix-recompress-340",
            "rank_hmatrix_generation": "hmatrix-recompress-339",
            "result_svd_basis": "mass-weighted",
            "result_tolerance": 1.0e-2,
            "result_tolerance_norm": "frobenius-absolute",
            "result_block_ranks_before": [8, 10, 12],
            "result_block_ranks_after": [9, 11, 13],
            "result_row_permutation": [0, 1, 2],
            "result_column_permutation": [0, 1, 2],
            "result_operator_relative_error": 0.1,
            "result_mesh_sha256": "b" * 64,
            "accepted_result_sha256": "c" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hmatrix_recompression_uses_current_svd_tolerance_norm_ranks_permutations_operator_mesh_and_result"
    ]


def test_v30_public_cq_restart_block_history_startup_weights_time_index_sample_count_digest_mismatch() -> None:
    summary = _summary_v30()
    identity = summary[
        "cq_restart_block_history_startup_weight_time_index_sample_contour_operator_result_identity"
    ]
    identity.update(
        {
            "block_cq_generation": "cq-restart-340",
            "history_cq_generation": "cq-restart-339",
            "result_block_size": 8,
            "result_completed_block_ids": [0, 2],
            "result_history_sample_count": 40,
            "result_startup_weights_ri": [[0.0, 0.0]],
            "result_restart_time_index": 47,
            "result_total_sample_count": 64,
            "result_contour_owner_sha256": "d" * 64,
            "result_operator_owner_sha256": "e" * 64,
            "loaded_history_sha256": "f" * 64,
            "accepted_result_sha256": "0" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_block_restart_uses_current_blocks_history_startup_weights_time_samples_owners_and_result"
    ]
