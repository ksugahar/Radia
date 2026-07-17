from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v24 import _summary_v24


_PROMOTED_CASE_IDS = (
    "v25_public_bem_demag_surface_orientation_magnetization_volume_material_frame_mismatch",
    "v25_public_linear_motor_thrust_ripple_period_position_phase_frame_generation_mismatch",
)


def _summary_v25():
    summary = _summary_v24()
    identity = summary["artifact_identity"]
    normals = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ]
    identity[
        "bem_demag_surface_orientation_magnetization_volume_material_frame_generation_identity"
    ] = {
        "solve_generation": "bem-demag-201",
        "surface_mesh_solve_generation": "bem-demag-201",
        "surface_orientation_solve_generation": "bem-demag-201",
        "magnetization_solve_generation": "bem-demag-201",
        "body_volume_solve_generation": "bem-demag-201",
        "material_region_solve_generation": "bem-demag-201",
        "coordinate_frame_solve_generation": "bem-demag-201",
        "result_solve_generation": "bem-demag-201",
        "surface_ids": [101, 102, 103, 104, 105, 106],
        "result_surface_ids": [101, 102, 103, 104, 105, 106],
        "surface_orientation": "outward_from_magnet",
        "result_surface_orientation": "outward_from_magnet",
        "outward_normals": normals,
        "result_outward_normals": normals,
        "magnetization_vector_a_per_m": [0.0, 0.0, 900000.0],
        "result_magnetization_vector_a_per_m": [0.0, 0.0, 900000.0],
        "body_volume_m3": 1.0e-6,
        "result_body_volume_m3": 1.0e-6,
        "material_region_id": 7,
        "result_material_region_id": 7,
        "coordinate_frame": "global_xyz",
        "result_coordinate_frame": "global_xyz",
        "surface_mesh_sha256": "1" * 64,
        "result_surface_mesh_sha256": "1" * 64,
        "demag_result_sha256": "2" * 64,
        "reported_demag_result_sha256": "2" * 64,
    }
    positions = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
    phases = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0, 360.0]
    thrust = [100.0, 102.0, 101.0, 99.0, 98.0, 100.0, 100.0]
    identity[
        "linear_motor_thrust_ripple_period_position_phase_frame_generation_identity"
    ] = {
        "sweep_generation": "linear-thrust-201",
        "period_sweep_generation": "linear-thrust-201",
        "position_sweep_generation": "linear-thrust-201",
        "phase_sweep_generation": "linear-thrust-201",
        "force_frame_sweep_generation": "linear-thrust-201",
        "sample_order_sweep_generation": "linear-thrust-201",
        "result_sweep_generation": "linear-thrust-201",
        "mechanical_period_m": 0.03,
        "result_mechanical_period_m": 0.03,
        "mover_positions_m": positions,
        "result_mover_positions_m": positions,
        "excitation_phase_deg": phases,
        "result_excitation_phase_deg": phases,
        "sample_order": list(range(7)),
        "result_sample_order": list(range(7)),
        "force_coordinate_frame": "global_x",
        "result_force_coordinate_frame": "global_x",
        "thrust_samples_n": thrust,
        "result_thrust_samples_n": thrust,
        "thrust_ripple_peak_to_peak_n": 4.0,
        "reported_thrust_ripple_peak_to_peak_n": 4.0,
        "thrust_table_sha256": "3" * 64,
        "result_thrust_table_sha256": "3" * 64,
    }
    return summary


def test_v25_public_positive_bem_demag_and_linear_motor_identity() -> None:
    assert magnetic_force_method_profile_gate(_summary_v25())["status"] == "ok"


def test_v25_public_bem_demag_surface_material_frame_mismatch() -> None:
    summary = _summary_v25()
    identity = summary["artifact_identity"][
        "bem_demag_surface_orientation_magnetization_volume_material_frame_generation_identity"
    ]
    identity.update(
        {
            "surface_orientation_solve_generation": "bem-demag-200",
            "magnetization_solve_generation": "bem-demag-199",
            "body_volume_solve_generation": "bem-demag-198",
            "material_region_solve_generation": "bem-demag-197",
            "coordinate_frame_solve_generation": "bem-demag-196",
            "result_surface_orientation": "inward_to_magnet",
            "result_outward_normals": [[-1.0, 0.0, 0.0]],
            "result_magnetization_vector_a_per_m": [0.0, 900000.0, 0.0],
            "result_body_volume_m3": 2.0e-6,
            "result_material_region_id": 8,
            "result_coordinate_frame": "local_xyz",
            "result_surface_mesh_sha256": "7" * 64,
            "reported_demag_result_sha256": "8" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_demag_surface_shares_orientation_magnetization_volume_material_and_frame"
    ]


def test_v25_public_linear_motor_thrust_ripple_generation_mismatch() -> None:
    summary = _summary_v25()
    identity = summary["artifact_identity"][
        "linear_motor_thrust_ripple_period_position_phase_frame_generation_identity"
    ]
    identity.update(
        {
            "period_sweep_generation": "linear-thrust-200",
            "position_sweep_generation": "linear-thrust-199",
            "phase_sweep_generation": "linear-thrust-198",
            "force_frame_sweep_generation": "linear-thrust-197",
            "sample_order_sweep_generation": "linear-thrust-196",
            "result_mechanical_period_m": 0.04,
            "result_mover_positions_m": [0.03, 0.0, 0.01],
            "result_excitation_phase_deg": [0.0, 120.0, 60.0],
            "result_sample_order": [6, 0, 2],
            "result_force_coordinate_frame": "mover_local_x",
            "result_thrust_samples_n": [100.0, 98.0, 102.0],
            "reported_thrust_ripple_peak_to_peak_n": 2.0,
            "result_thrust_table_sha256": "9" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "linear_motor_thrust_ripple_shares_period_position_phase_frame_and_order"
    ]
