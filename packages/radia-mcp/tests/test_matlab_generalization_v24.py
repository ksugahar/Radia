from __future__ import annotations

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v23 import _summary_v23


_PROMOTED_CASE_IDS = (
    "v24_public_fembem_trace_normal_interface_node_order_unit_generation_mismatch",
    "v24_public_cq_contour_weight_startup_causality_window_result_generation_mismatch",
)


def _summary_v24():
    summary = _summary_v23()
    summary["fembem_trace_normal_interface_node_order_unit_generation_identity"] = {
        "coupling_generation": "fembem-101",
        "trace_coupling_generation": "fembem-101",
        "normal_coupling_generation": "fembem-101",
        "node_order_coupling_generation": "fembem-101",
        "unit_coupling_generation": "fembem-101",
        "operator_coupling_generation": "fembem-101",
        "result_coupling_generation": "fembem-101",
        "trace_orientation": "volume_to_boundary",
        "result_trace_orientation": "volume_to_boundary",
        "outward_normal_convention": "exterior_from_volume",
        "result_outward_normal_convention": "exterior_from_volume",
        "interface_node_ids": [1, 2, 3, 4],
        "result_interface_node_ids": [1, 2, 3, 4],
        "boundary_triangles": [[1, 2, 3], [1, 4, 2]],
        "result_boundary_triangles": [[1, 2, 3], [1, 4, 2]],
        "outward_normals": [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
        "result_outward_normals": [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
        "physical_units": {"pressure": "Pa", "normal_velocity": "m/s"},
        "result_physical_units": {"pressure": "Pa", "normal_velocity": "m/s"},
        "interface_mesh_sha256": "1" * 64,
        "result_interface_mesh_sha256": "1" * 64,
        "coupled_operator_sha256": "2" * 64,
        "result_coupled_operator_sha256": "2" * 64,
    }
    summary[
        "cq_contour_weight_startup_causality_window_result_generation_identity"
    ] = {
        "cq_generation": "cq-time-101",
        "contour_cq_generation": "cq-time-101",
        "weight_cq_generation": "cq-time-101",
        "startup_cq_generation": "cq-time-101",
        "causality_window_cq_generation": "cq-time-101",
        "time_grid_cq_generation": "cq-time-101",
        "result_cq_generation": "cq-time-101",
        "method": "BDF2",
        "result_method": "BDF2",
        "contour_points_ri": [[0.8, 0.0], [0.0, 0.8], [-0.8, 0.0], [0.0, -0.8]],
        "result_contour_points_ri": [[0.8, 0.0], [0.0, 0.8], [-0.8, 0.0], [0.0, -0.8]],
        "cq_weights_ri": [[1.5, 0.0], [-2.0, 0.0], [0.5, 0.0], [0.0, 0.0]],
        "result_cq_weights_ri": [[1.5, 0.0], [-2.0, 0.0], [0.5, 0.0], [0.0, 0.0]],
        "startup_weights_ri": [[1.0, 0.0], [-1.0, 0.0]],
        "result_startup_weights_ri": [[1.0, 0.0], [-1.0, 0.0]],
        "time_samples_s": [0.0, 0.001, 0.002, 0.003],
        "result_time_samples_s": [0.0, 0.001, 0.002, 0.003],
        "causality_window_s": [0.0, 0.003],
        "result_causality_window_s": [0.0, 0.003],
        "prehistory_norm": 0.0,
        "result_prehistory_norm": 0.0,
        "cq_result_sha256": "3" * 64,
        "reported_cq_result_sha256": "3" * 64,
    }
    return summary


def test_v24_public_positive_fembem_interface_and_cq_time_identity() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v24())["status"] == "ok"


def test_v24_public_fembem_trace_normal_node_order_unit_mismatch() -> None:
    summary = _summary_v24()
    summary["fembem_trace_normal_interface_node_order_unit_generation_identity"].update(
        {
            "trace_coupling_generation": "fembem-100",
            "normal_coupling_generation": "fembem-99",
            "node_order_coupling_generation": "fembem-98",
            "unit_coupling_generation": "fembem-97",
            "operator_coupling_generation": "fembem-96",
            "result_trace_orientation": "boundary_to_volume",
            "result_outward_normal_convention": "inward_to_volume",
            "result_interface_node_ids": [1, 3, 2, 4],
            "result_boundary_triangles": [[1, 3, 2], [1, 2, 4]],
            "result_outward_normals": [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
            "result_physical_units": {"pressure": "kPa", "normal_velocity": "mm/s"},
            "result_interface_mesh_sha256": "a" * 64,
            "result_coupled_operator_sha256": "b" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_trace_uses_current_normals_nodes_units_mesh_and_operator"
    ]


def test_v24_public_cq_contour_weight_startup_causality_window_mismatch() -> None:
    summary = _summary_v24()
    summary[
        "cq_contour_weight_startup_causality_window_result_generation_identity"
    ].update(
        {
            "contour_cq_generation": "cq-time-100",
            "weight_cq_generation": "cq-time-99",
            "startup_cq_generation": "cq-time-98",
            "causality_window_cq_generation": "cq-time-97",
            "time_grid_cq_generation": "cq-time-96",
            "result_contour_points_ri": [[1.2, 0.0], [0.0, 0.7], [-0.7, 0.0], [0.0, -0.7]],
            "result_cq_weights_ri": [[1.0, 0.0], [-1.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            "result_startup_weights_ri": [[0.0, 0.0], [0.0, 0.0]],
            "result_time_samples_s": [0.001, 0.002, 0.003, 0.004],
            "result_causality_window_s": [-0.001, 0.003],
            "result_prehistory_norm": 0.2,
            "reported_cq_result_sha256": "c" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_time_history_uses_current_contour_weights_startup_and_causality_window"
    ]
