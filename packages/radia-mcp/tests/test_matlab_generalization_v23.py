from __future__ import annotations

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_regularized_trace_inverse_gate import _summary, _with_v22_identity


def _summary_v23():
    summary = _with_v22_identity(_summary())
    summary["parallel_pool_worker_path_device_rng_code_generation_identity"] = {
        "pool_generation": "parallel-pool-51",
        "worker_path_pool_generation": "parallel-pool-51",
        "device_pool_generation": "parallel-pool-51",
        "rng_pool_generation": "parallel-pool-51",
        "code_pool_generation": "parallel-pool-51",
        "result_pool_generation": "parallel-pool-51",
        "worker_ids": [1, 2, 3, 4],
        "result_worker_ids": [1, 2, 3, 4],
        "worker_code_paths": ["toolbox/a"] * 4,
        "result_worker_code_paths": ["toolbox/a"] * 4,
        "device_assignments": ["cpu:0", "cpu:1", "cpu:2", "cpu:3"],
        "result_device_assignments": ["cpu:0", "cpu:1", "cpu:2", "cpu:3"],
        "random_stream_seeds": [101, 202, 303, 404],
        "result_random_stream_seeds": [101, 202, 303, 404],
        "worker_code_sha256": "1" * 64,
        "result_worker_code_sha256": "1" * 64,
        "parallel_result_sha256": "2" * 64,
        "assembled_parallel_result_sha256": "2" * 64,
    }
    summary["autodiff_tape_variable_order_mesh_objective_generation_identity"] = {
        "tape_generation": "autodiff-tape-51",
        "variable_order_tape_generation": "autodiff-tape-51",
        "mesh_tape_generation": "autodiff-tape-51",
        "objective_scaling_tape_generation": "autodiff-tape-51",
        "primal_solve_tape_generation": "autodiff-tape-51",
        "gradient_result_tape_generation": "autodiff-tape-51",
        "variable_ids": ["radius", "thickness", "impedance"],
        "gradient_variable_ids": ["radius", "thickness", "impedance"],
        "mesh_sha256": "3" * 64,
        "gradient_mesh_sha256": "3" * 64,
        "objective_id": "radiated_power",
        "gradient_objective_id": "radiated_power",
        "objective_scale": 0.001,
        "gradient_objective_scale": 0.001,
        "primal_state_sha256": "4" * 64,
        "gradient_primal_state_sha256": "4" * 64,
        "gradient_table_sha256": "5" * 64,
        "reported_gradient_table_sha256": "5" * 64,
    }
    return summary


def test_v23_public_positive_parallel_pool_and_autodiff_identity():
    assert regularized_trace_inverse_path_gate(_summary_v23())["status"] == "ok"


def test_v23_public_parallel_pool_worker_path_device_rng_code_generation_mismatch():
    summary = _summary_v23()
    summary["parallel_pool_worker_path_device_rng_code_generation_identity"].update(
        {
            "worker_path_pool_generation": "parallel-pool-50",
            "device_pool_generation": "parallel-pool-49",
            "rng_pool_generation": "parallel-pool-48",
            "code_pool_generation": "parallel-pool-47",
            "result_pool_generation": "parallel-pool-46",
            "result_worker_ids": [1, 2, 4],
            "result_worker_code_paths": ["toolbox/a", "toolbox/old", "toolbox/a"],
            "result_device_assignments": ["cpu:0", "gpu:0", "cpu:3"],
            "result_random_stream_seeds": [101, 999, 404],
            "result_worker_code_sha256": "c" * 64,
            "assembled_parallel_result_sha256": "d" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "parallel_results_use_current_worker_paths_devices_rng_and_code"
    ]


def test_v23_public_autodiff_tape_variable_order_mesh_objective_generation_mismatch():
    summary = _summary_v23()
    summary["autodiff_tape_variable_order_mesh_objective_generation_identity"].update(
        {
            "variable_order_tape_generation": "autodiff-tape-50",
            "mesh_tape_generation": "autodiff-tape-49",
            "objective_scaling_tape_generation": "autodiff-tape-48",
            "primal_solve_tape_generation": "autodiff-tape-47",
            "gradient_result_tape_generation": "autodiff-tape-46",
            "gradient_variable_ids": ["thickness", "radius", "impedance"],
            "gradient_mesh_sha256": "e" * 64,
            "gradient_objective_id": "mass",
            "gradient_objective_scale": 1.0,
            "gradient_primal_state_sha256": "f" * 64,
            "reported_gradient_table_sha256": "0" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "autodiff_gradients_use_current_tape_variables_mesh_objective_and_primal"
    ]
