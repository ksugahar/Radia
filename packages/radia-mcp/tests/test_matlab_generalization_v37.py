from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v36 import _summary_v36


_PROMOTED_CASE_IDS = (
    "v37_public_hmatrix_dense_reference_error_memory_rank_complexity_mesh_owner_mismatch",
    "v37_public_multifrequency_fembem_adjoint_gradient_quadrature_trace_fd_owner_mismatch",
)


def _summary_v37():
    payload = deepcopy(_summary_v36())
    generation = "hmatrix-benchmark-513"
    mirrored = {
        "boundary_unknown_count": [200, 400, 800],
        "dense_reference_relative_error": [8.0e-4, 6.0e-4, 4.0e-4],
        "relative_tolerance": 1.0e-3,
        "maximum_block_rank": [8, 12, 16],
        "dense_memory_bytes": [640000, 2560000, 10240000],
        "hmatrix_memory_bytes": [80000, 190000, 450000],
        "memory_complexity_exponent": 1.25,
        "rank_complexity_exponent": 0.5,
        "boundary_mesh_sha256": "1" * 64,
    }
    payload["hmatrix_dense_reference_error_tolerance_rank_memory_complexity_mesh_operator_benchmark_result_identity"] = {
        "hmatrix_generation": generation,
        **{key: generation for key in (
            "dense_generation", "tolerance_generation", "rank_generation", "memory_generation",
            "complexity_generation", "mesh_generation", "operator_generation", "benchmark_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "operator_owner": "hmatrix/operator-513", "accepted_operator_owner": "hmatrix/operator-513",
        "benchmark_owner": "hmatrix/benchmark-513", "accepted_benchmark_owner": "hmatrix/benchmark-513",
        "hmatrix_result_sha256": "2" * 64, "accepted_hmatrix_result_sha256": "2" * 64,
    }
    generation = "multifrequency-adjoint-513"
    mirrored = {
        "frequency_hz": [100.0, 200.0, 400.0], "frequency_weights": [0.2, 0.3, 0.5],
        "objective_complex": [[1.0, 0.1], [2.0, 0.2], [3.0, 0.3]],
        "weighted_objective_complex": [2.3, 0.23], "quadrature_order": [4, 4, 6],
        "trace_node_map": [1, 2, 3], "frequency_gradient": [2.0, 4.0, 6.0],
        "accumulated_gradient": 4.6, "finite_difference_gradient": 4.600001,
        "gradient_relative_tolerance": 1.0e-5, "fembem_mesh_sha256": "3" * 64,
    }
    payload["multifrequency_fembem_adjoint_weight_objective_quadrature_trace_gradient_fd_mesh_owner_result_identity"] = {
        "adjoint_generation": generation,
        **{key: generation for key in (
            "frequency_generation", "weight_generation", "objective_generation", "quadrature_generation",
            "trace_generation", "gradient_generation", "fd_generation", "mesh_generation", "owner_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "adjoint_owner": "fembem/adjoint-513", "accepted_adjoint_owner": "fembem/adjoint-513",
        "adjoint_result_sha256": "4" * 64, "accepted_adjoint_result_sha256": "4" * 64,
    }
    return payload


def test_v37_public_positive_hmatrix_and_multifrequency_adjoint_closure():
    assert regularized_trace_inverse_path_gate(_summary_v37())["status"] == "ok"


def test_v37_public_hmatrix_dense_reference_error_memory_rank_complexity_mesh_owner_mismatch():
    payload = _summary_v37()
    identity = payload["hmatrix_dense_reference_error_tolerance_rank_memory_complexity_mesh_operator_benchmark_result_identity"]
    identity.update({
        "dense_generation": "hmatrix-benchmark-512", "memory_generation": "hmatrix-benchmark-511",
        "result_generation": "hmatrix-benchmark-510", "result_dense_reference_relative_error": [9.0],
        "result_relative_tolerance": -1.0, "result_maximum_block_rank": [999],
        "result_dense_memory_bytes": [1], "result_hmatrix_memory_bytes": [999999999],
        "result_memory_complexity_exponent": 3.0, "result_rank_complexity_exponent": 2.0,
        "result_boundary_mesh_sha256": "8" * 64, "accepted_operator_owner": "hmatrix/old",
        "accepted_benchmark_owner": "benchmark/old", "accepted_hmatrix_result_sha256": "9" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["hmatrix_benchmarks_use_current_dense_error_tolerance_rank_memory_complexity_mesh_owners_and_result"]


def test_v37_public_multifrequency_fembem_adjoint_gradient_quadrature_trace_fd_owner_mismatch():
    payload = _summary_v37()
    identity = payload["multifrequency_fembem_adjoint_weight_objective_quadrature_trace_gradient_fd_mesh_owner_result_identity"]
    identity.update({
        "weight_generation": "multifrequency-adjoint-512", "trace_generation": "multifrequency-adjoint-511",
        "result_generation": "multifrequency-adjoint-510", "result_frequency_weights": [2.0, -1.0],
        "result_objective_complex": [[99.0, 99.0]], "result_weighted_objective_complex": [99.0, 99.0],
        "result_quadrature_order": [0], "result_trace_node_map": [3, 2, 1],
        "result_frequency_gradient": [99.0], "result_accumulated_gradient": -99.0,
        "result_finite_difference_gradient": 99.0, "result_fembem_mesh_sha256": "a" * 64,
        "accepted_adjoint_owner": "fembem/old", "accepted_adjoint_result_sha256": "b" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["multifrequency_fembem_adjoints_use_current_weights_objective_quadrature_trace_gradient_fd_mesh_owner_and_result"]


def test_v37_public_rejects_self_consistent_dense_storage_as_hmatrix():
    payload = _summary_v37()
    identity = payload["hmatrix_dense_reference_error_tolerance_rank_memory_complexity_mesh_operator_benchmark_result_identity"]
    identity["hmatrix_memory_bytes"] = list(identity["dense_memory_bytes"])
    identity["result_hmatrix_memory_bytes"] = list(identity["dense_memory_bytes"])
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_wrong_gradient_accumulation():
    payload = _summary_v37()
    identity = payload["multifrequency_fembem_adjoint_weight_objective_quadrature_trace_gradient_fd_mesh_owner_result_identity"]
    identity["accumulated_gradient"] = 99.0
    identity["result_accumulated_gradient"] = 99.0
    identity["finite_difference_gradient"] = 99.0
    identity["result_finite_difference_gradient"] = 99.0
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"
