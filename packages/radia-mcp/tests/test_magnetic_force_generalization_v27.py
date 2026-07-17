from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v26 import _summary_v26


def _summary_v27():
    summary = _summary_v26()
    identity = summary["artifact_identity"]
    normals = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
    ]
    identity["bem_panel_orientation_magnetization_frame_self_term_energy_force_generation_identity"] = {
        "panel_generation": "bem-panel-311",
        "orientation_panel_generation": "bem-panel-311",
        "magnetization_panel_generation": "bem-panel-311",
        "self_term_panel_generation": "bem-panel-311",
        "energy_panel_generation": "bem-panel-311",
        "force_panel_generation": "bem-panel-311",
        "mesh_panel_generation": "bem-panel-311",
        "result_panel_generation": "bem-panel-311",
        "panel_ids": [11, 12, 13, 14],
        "result_panel_ids": [11, 12, 13, 14],
        "outward_unit_normals": normals,
        "result_outward_unit_normals": normals,
        "panel_area_m2": [0.5, 0.5, 0.5, 0.5],
        "result_panel_area_m2": [0.5, 0.5, 0.5, 0.5],
        "magnetization_a_m": [[0.0, 0.0, 800000.0]] * 4,
        "result_magnetization_a_m": [[0.0, 0.0, 800000.0]] * 4,
        "magnetization_frame": "global-cartesian",
        "result_magnetization_frame": "global-cartesian",
        "singular_self_term": "analytic-solid-angle",
        "result_singular_self_term": "analytic-solid-angle",
        "displacement_m": [-0.001, 0.0, 0.001],
        "result_displacement_m": [-0.001, 0.0, 0.001],
        "magnetic_energy_j": [0.005, 0.0, 0.005],
        "result_magnetic_energy_j": [0.005, 0.0, 0.005],
        "negative_energy_derivative_force_n": [10.0, 0.0, -10.0],
        "result_force_n": [10.0, 0.0, -10.0],
        "panel_mesh_sha256": "1" * 64,
        "result_panel_mesh_sha256": "1" * 64,
        "force_result_sha256": "2" * 64,
        "accepted_force_result_sha256": "2" * 64,
    }
    identity["motor_reduced_basis_snapshot_operating_point_interpolation_torque_residual_generation_identity"] = {
        "reduced_generation": "motor-rom-311",
        "basis_reduced_generation": "motor-rom-311",
        "snapshot_reduced_generation": "motor-rom-311",
        "operating_point_reduced_generation": "motor-rom-311",
        "weight_reduced_generation": "motor-rom-311",
        "torque_reduced_generation": "motor-rom-311",
        "residual_reduced_generation": "motor-rom-311",
        "result_reduced_generation": "motor-rom-311",
        "basis_dimension": 3,
        "result_basis_dimension": 3,
        "snapshot_ids": ["snap-a", "snap-b", "snap-c"],
        "result_snapshot_ids": ["snap-a", "snap-b", "snap-c"],
        "snapshot_operating_points": [[1000.0, 10.0, 0.0], [2000.0, 20.0, 30.0], [3000.0, 30.0, 60.0]],
        "result_snapshot_operating_points": [[1000.0, 10.0, 0.0], [2000.0, 20.0, 30.0], [3000.0, 30.0, 60.0]],
        "query_operating_point": [2200.0, 22.0, 36.0],
        "result_query_operating_point": [2200.0, 22.0, 36.0],
        "interpolation_weights": [0.2, 0.6, 0.2],
        "result_interpolation_weights": [0.2, 0.6, 0.2],
        "snapshot_torque_nm": [1.0, 2.0, 3.0],
        "result_snapshot_torque_nm": [1.0, 2.0, 3.0],
        "reduced_torque_nm": 2.0,
        "result_reduced_torque_nm": 2.0,
        "relative_residual": 1.0e-6,
        "accepted_relative_residual": 1.0e-4,
        "result_relative_residual": 1.0e-6,
        "basis_sha256": "3" * 64,
        "loaded_basis_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v27_public_positive_bem_panel_and_motor_reduced_basis_identity():
    assert magnetic_force_method_profile_gate(_summary_v27())["status"] == "ok"


def test_v27_public_demag_bem_panel_orientation_magnetization_frame_self_term_energy_force_mismatch():
    summary = _summary_v27()
    summary["artifact_identity"][
        "bem_panel_orientation_magnetization_frame_self_term_energy_force_generation_identity"
    ].update({
        "orientation_panel_generation": "bem-panel-310",
        "self_term_panel_generation": "bem-panel-309",
        "mesh_panel_generation": "bem-panel-308",
        "result_outward_unit_normals": [[-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        "result_panel_area_m2": [0.5, -0.5],
        "result_magnetization_a_m": [[800000.0, 0.0, 0.0]],
        "result_magnetization_frame": "panel-local",
        "result_singular_self_term": "omitted",
        "result_magnetic_energy_j": [0.0, 0.005, 0.0],
        "result_force_n": [-10.0, 0.0, 10.0],
        "result_panel_mesh_sha256": "8" * 64,
        "accepted_force_result_sha256": "9" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_demag_force_uses_current_panel_orientation_magnetization_self_term_energy_mesh_and_result"
    ]


def test_v27_public_motor_reduced_basis_snapshot_operating_point_interpolation_torque_residual_mismatch():
    summary = _summary_v27()
    summary["artifact_identity"][
        "motor_reduced_basis_snapshot_operating_point_interpolation_torque_residual_generation_identity"
    ].update({
        "basis_reduced_generation": "motor-rom-310",
        "snapshot_reduced_generation": "motor-rom-309",
        "weight_reduced_generation": "motor-rom-308",
        "result_basis_dimension": 2,
        "result_snapshot_ids": ["snap-a", "snap-old"],
        "result_snapshot_operating_points": [[1000.0, 10.0, 0.0]],
        "result_query_operating_point": [4000.0, 50.0, 90.0],
        "result_interpolation_weights": [0.5, 0.5, 0.5],
        "result_snapshot_torque_nm": [1.0, 2.0],
        "result_reduced_torque_nm": 4.0,
        "result_relative_residual": 0.2,
        "loaded_basis_sha256": "a" * 64,
        "accepted_result_sha256": "b" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_reduced_torque_uses_current_basis_snapshots_operating_point_weights_residual_and_result"
    ]
