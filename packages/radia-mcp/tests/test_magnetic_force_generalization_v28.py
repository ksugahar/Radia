from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v27 import _summary_v27


def _summary_v28():
    summary = _summary_v27()
    identity = summary["artifact_identity"]
    identity[
        "maglev_force_stiffness_displacement_step_coordinate_mesh_solution_derivative_generation_identity"
    ] = {
        "stiffness_generation": "maglev-stiffness-321",
        "displacement_stiffness_generation": "maglev-stiffness-321",
        "coordinate_stiffness_generation": "maglev-stiffness-321",
        "geometry_stiffness_generation": "maglev-stiffness-321",
        "mesh_stiffness_generation": "maglev-stiffness-321",
        "force_stiffness_generation": "maglev-stiffness-321",
        "derivative_stiffness_generation": "maglev-stiffness-321",
        "solution_stiffness_generation": "maglev-stiffness-321",
        "result_stiffness_generation": "maglev-stiffness-321",
        "displacement_m": [-0.001, 0.0, 0.001],
        "result_displacement_m": [-0.001, 0.0, 0.001],
        "displacement_step_m": 0.001,
        "result_displacement_step_m": 0.001,
        "coordinate_direction": "global-z-positive",
        "result_coordinate_direction": "global-z-positive",
        "force_n": [12.0, 10.0, 8.0],
        "result_force_n": [12.0, 10.0, 8.0],
        "derivative_convention": "stiffness-equals-negative-force-derivative",
        "result_derivative_convention": "stiffness-equals-negative-force-derivative",
        "stiffness_n_m": 2000.0,
        "result_stiffness_n_m": 2000.0,
        "geometry_sha256": "1" * 64,
        "result_geometry_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "solution_sha256": "3" * 64,
        "accepted_solution_sha256": "3" * 64,
    }
    identity[
        "motor_winding_harmonic_current_phase_rotor_angle_coenergy_torque_result_generation_identity"
    ] = {
        "torque_generation": "motor-coenergy-321",
        "winding_torque_generation": "motor-coenergy-321",
        "harmonic_torque_generation": "motor-coenergy-321",
        "current_torque_generation": "motor-coenergy-321",
        "phase_torque_generation": "motor-coenergy-321",
        "angle_torque_generation": "motor-coenergy-321",
        "coenergy_torque_generation": "motor-coenergy-321",
        "mesh_torque_generation": "motor-coenergy-321",
        "result_torque_generation": "motor-coenergy-321",
        "phase_order": ["U", "V", "W"],
        "result_phase_order": ["U", "V", "W"],
        "harmonic_orders": [1, 5, 7],
        "result_harmonic_orders": [1, 5, 7],
        "phase_current_harmonic_a": [
            [100.0, -50.0, -50.0],
            [8.0, -4.0, -4.0],
            [5.0, -2.5, -2.5],
        ],
        "result_phase_current_harmonic_a": [
            [100.0, -50.0, -50.0],
            [8.0, -4.0, -4.0],
            [5.0, -2.5, -2.5],
        ],
        "current_phase_deg": [0.0, -120.0, 120.0],
        "result_current_phase_deg": [0.0, -120.0, 120.0],
        "rotor_mechanical_angle_deg": [0.0, 1.0, 2.0],
        "result_rotor_mechanical_angle_deg": [0.0, 1.0, 2.0],
        "coenergy_j": [0.100, 0.102, 0.104],
        "result_coenergy_j": [0.100, 0.102, 0.104],
        "torque_convention": "positive-coenergy-angle-derivative",
        "result_torque_convention": "positive-coenergy-angle-derivative",
        "torque_nm": [0.11459155902616465, 0.11459155902616465],
        "result_torque_nm": [0.11459155902616465, 0.11459155902616465],
        "mesh_sha256": "4" * 64,
        "result_mesh_sha256": "4" * 64,
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    return summary


def test_v28_public_positive_maglev_stiffness_and_motor_coenergy_identities():
    assert magnetic_force_method_profile_gate(_summary_v28())["status"] == "ok"


def test_v28_public_maglev_force_stiffness_displacement_step_coordinate_mesh_solution_derivative_mismatch():
    summary = _summary_v28()
    identity = summary["artifact_identity"][
        "maglev_force_stiffness_displacement_step_coordinate_mesh_solution_derivative_generation_identity"
    ]
    identity.update(
        {
            "displacement_stiffness_generation": "maglev-stiffness-320",
            "mesh_stiffness_generation": "maglev-stiffness-319",
            "result_displacement_m": [0.0, 0.002, 0.004],
            "result_displacement_step_m": 0.002,
            "result_coordinate_direction": "local-x-negative",
            "result_force_n": [8.0, 10.0, 12.0],
            "result_derivative_convention": "stiffness-equals-force-derivative",
            "result_stiffness_n_m": -1000.0,
            "result_geometry_sha256": "a" * 64,
            "result_mesh_sha256": "b" * 64,
            "accepted_solution_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "maglev_stiffness_uses_current_displacement_coordinate_force_geometry_mesh_and_solution"
    ]


def test_v28_public_motor_winding_harmonic_current_phase_rotor_angle_coenergy_torque_result_mismatch():
    summary = _summary_v28()
    identity = summary["artifact_identity"][
        "motor_winding_harmonic_current_phase_rotor_angle_coenergy_torque_result_generation_identity"
    ]
    identity.update(
        {
            "winding_torque_generation": "motor-coenergy-320",
            "angle_torque_generation": "motor-coenergy-319",
            "result_phase_order": ["U", "W", "V"],
            "result_harmonic_orders": [1, 3, 5],
            "result_phase_current_harmonic_a": [[100.0, -40.0, -50.0]],
            "result_current_phase_deg": [0.0, 120.0, -120.0],
            "result_rotor_mechanical_angle_deg": [0.0, 2.0, 4.0],
            "result_coenergy_j": [0.100, 0.101, 0.099],
            "result_torque_convention": "negative-coenergy-angle-derivative",
            "result_torque_nm": [-0.03, 0.06],
            "result_mesh_sha256": "d" * 64,
            "accepted_result_sha256": "e" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_torque_uses_current_winding_harmonics_currents_phases_angles_coenergy_mesh_and_result"
    ]
