from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v40 import _summary_v40


_PROMOTED_CASE_IDS = (
    "v41_public_cq_acoustic_causality_passivity_timestep_ztransform_energy_history_mismatch",
    "v41_public_hmatrix_admissibility_cluster_rank_tolerance_matvec_error_memory_mismatch",
)
_CQ_KEY = (
    "cq_acoustic_causality_passivity_timestep_ztransform_energy_history_mesh_"
    "owner_result_identity"
)
_HMATRIX_KEY = (
    "hmatrix_admissibility_cluster_rank_tolerance_matvec_error_memory_mesh_"
    "owner_result_identity"
)


def _summary_v41() -> dict:
    payload = deepcopy(_summary_v40())
    generation = "cq-acoustic-731"
    mirrored = {
        "multistep_method": "bdf2",
        "time_step_s": 2.5e-4,
        "z_transform_radius": 0.94,
        "multistep_symbol_samples": [[0.0, 0.0], [0.5, 0.2], [1.5, 0.0]],
        "laplace_frequency_samples_rad_s": [[20.0, 0.0], [35.0, 80.0], [60.0, 0.0]],
        "excitation_history": [0.0, 1.0, 0.4, 0.1, 0.0],
        "pressure_history": [0.0, 0.2, 0.3, 0.15, 0.04],
        "causal_prefix_length": 1,
        "minimum_passivity_real_part": 0.015,
        "boundary_work_j": 0.012,
        "radiated_energy_j": 0.010,
        "dissipated_energy_j": 0.002,
        "energy_balance_residual_j": 0.0,
        "energy_balance_tolerance_j": 1.0e-8,
        "boundary_mesh_sha256": "1" * 64,
    }
    payload[_CQ_KEY] = {
        "cq_generation": generation,
        **{key: generation for key in (
            "multistep_generation", "timestep_generation", "ztransform_generation",
            "frequency_generation", "history_generation", "passivity_generation",
            "energy_generation", "mesh_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "cq_owner": "acoustic/cq-731",
        "accepted_cq_owner": "acoustic/cq-731",
        "cq_result_sha256": "2" * 64,
        "accepted_cq_result_sha256": "2" * 64,
    }

    generation = "hmatrix-731"
    mirrored = {
        "cluster_leaf_size": 32,
        "cluster_permutation": [3, 1, 4, 2],
        "admissibility_eta": 2.0,
        "block_partition": [["low_rank", "dense"], ["dense", "low_rank"]],
        "numerical_ranks": [[4, 0], [0, 3]],
        "compression_relative_tolerance": 1.0e-5,
        "measured_matvec_relative_error": 4.0e-6,
        "dense_memory_bytes": 131072,
        "compressed_memory_bytes": 32768,
        "boundary_mesh_sha256": "3" * 64,
    }
    payload[_HMATRIX_KEY] = {
        "hmatrix_generation": generation,
        **{key: generation for key in (
            "cluster_generation", "admissibility_generation", "partition_generation",
            "rank_generation", "tolerance_generation", "matvec_generation",
            "memory_generation", "mesh_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "hmatrix_owner": "acoustic/hmatrix-731",
        "accepted_hmatrix_owner": "acoustic/hmatrix-731",
        "hmatrix_result_sha256": "4" * 64,
        "accepted_hmatrix_result_sha256": "4" * 64,
    }
    return payload


def test_v41_public_positive_cq_and_hmatrix_closure() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v41())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v41_public_cq_acoustic_causality_passivity_timestep_ztransform_energy_history_mismatch() -> None:
    payload = _summary_v41()
    payload[_CQ_KEY].update({
        "timestep_generation": "cq-acoustic-730",
        "result_time_step_s": -1.0,
        "result_z_transform_radius": 1.1,
        "result_pressure_history": [1.0, 9.0],
        "result_minimum_passivity_real_part": -1.0,
        "result_energy_balance_residual_j": 1.0,
        "accepted_cq_owner": "acoustic/old",
        "accepted_cq_result_sha256": "a" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_acoustic_history_uses_current_causality_passivity_timestep_ztransform_energy_mesh_owner_and_result"
    ]


def test_v41_public_hmatrix_admissibility_cluster_rank_tolerance_matvec_error_memory_mismatch() -> None:
    payload = _summary_v41()
    payload[_HMATRIX_KEY].update({
        "cluster_generation": "hmatrix-730",
        "result_cluster_permutation": [1, 1, 9],
        "result_admissibility_eta": -1.0,
        "result_numerical_ranks": [[-1]],
        "result_measured_matvec_relative_error": 1.0,
        "result_compressed_memory_bytes": 262144,
        "accepted_hmatrix_owner": "acoustic/old",
        "accepted_hmatrix_result_sha256": "b" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hmatrix_compression_uses_current_clusters_admissibility_ranks_tolerance_matvec_memory_mesh_owner_and_result"
    ]


def test_v41_public_rejects_self_consistent_noncausal_pressure_history() -> None:
    payload = _summary_v41()
    payload[_CQ_KEY]["pressure_history"] = [1.0, 0.2, 0.3, 0.15, 0.04]
    payload[_CQ_KEY]["result_pressure_history"] = [1.0, 0.2, 0.3, 0.15, 0.04]
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_hmatrix_error_or_memory_regression() -> None:
    payload = _summary_v41()
    payload[_HMATRIX_KEY]["measured_matvec_relative_error"] = 2.0e-5
    payload[_HMATRIX_KEY]["result_measured_matvec_relative_error"] = 2.0e-5
    payload[_HMATRIX_KEY]["compressed_memory_bytes"] = 131072
    payload[_HMATRIX_KEY]["result_compressed_memory_bytes"] = 131072
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"
