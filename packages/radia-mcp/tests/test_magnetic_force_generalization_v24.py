from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v23 import _summary_v23


_PROMOTED_CASE_IDS = (
    "v24_public_maglev_force_stiffness_equilibrium_energy_finite_difference_mismatch",
    "v24_public_motor_dual_lane_geometry_excitation_force_frame_harmonic_alignment_mismatch",
)


def _summary_v24():
    summary = _summary_v23()
    identity = summary["artifact_identity"]
    identity[
        "maglev_force_stiffness_equilibrium_energy_finite_difference_generation_identity"
    ] = {
        "maglev_generation": "maglev-101",
        "equilibrium_maglev_generation": "maglev-101",
        "displacement_maglev_generation": "maglev-101",
        "energy_maglev_generation": "maglev-101",
        "force_maglev_generation": "maglev-101",
        "stiffness_maglev_generation": "maglev-101",
        "coordinate_frame_maglev_generation": "maglev-101",
        "result_maglev_generation": "maglev-101",
        "equilibrium_displacement_m": 0.0,
        "result_equilibrium_displacement_m": 0.0,
        "equilibrium_force_n": 0.0,
        "result_equilibrium_force_n": 0.0,
        "displacement_samples_m": [-0.001, 0.0, 0.001],
        "energy_samples_j": [0.005, 0.0, 0.005],
        "force_samples_n": [10.0, 0.0, -10.0],
        "energy_finite_difference_force_n": [10.0, 0.0, -10.0],
        "stiffness_n_m": 10000.0,
        "reported_stiffness_n_m": 10000.0,
        "force_energy_sign_convention": "force=-dW/dx",
        "result_force_energy_sign_convention": "force=-dW/dx",
        "coordinate_frame": "global_z",
        "result_coordinate_frame": "global_z",
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "maglev_result_sha256": "2" * 64,
        "reported_maglev_result_sha256": "2" * 64,
    }
    identity[
        "motor_dual_lane_geometry_excitation_force_frame_harmonic_alignment_generation_identity"
    ] = {
        "comparison_generation": "dual-motor-101",
        "geometry_comparison_generation": "dual-motor-101",
        "excitation_comparison_generation": "dual-motor-101",
        "force_frame_comparison_generation": "dual-motor-101",
        "harmonic_comparison_generation": "dual-motor-101",
        "rotor_angle_comparison_generation": "dual-motor-101",
        "result_comparison_generation": "dual-motor-101",
        "lane_ids": ["ngsolve-age", "hdiv-mmm-hcurl-eddy-bubble"],
        "result_lane_ids": ["ngsolve-age", "hdiv-mmm-hcurl-eddy-bubble"],
        "geometry_revision_sha256": ["3" * 64, "3" * 64],
        "result_geometry_revision_sha256": ["3" * 64, "3" * 64],
        "excitation_table_sha256": ["4" * 64, "4" * 64],
        "result_excitation_table_sha256": ["4" * 64, "4" * 64],
        "force_coordinate_frames": ["rotor_dq", "rotor_dq"],
        "result_force_coordinate_frames": ["rotor_dq", "rotor_dq"],
        "harmonic_bins": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "result_harmonic_bins": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "rotor_angles_deg": [[0.0, 5.0, 10.0], [0.0, 5.0, 10.0]],
        "result_rotor_angles_deg": [[0.0, 5.0, 10.0], [0.0, 5.0, 10.0]],
        "force_harmonics_n": [
            [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
            [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
        ],
        "result_force_harmonics_n": [
            [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
            [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
        ],
        "comparison_result_sha256": "5" * 64,
        "reported_comparison_result_sha256": "5" * 64,
    }
    return summary


def test_v24_public_positive_maglev_energy_and_motor_dual_lane_identity() -> None:
    assert magnetic_force_method_profile_gate(_summary_v24())["status"] == "ok"


def test_v24_public_maglev_force_stiffness_energy_generation_mismatch() -> None:
    summary = _summary_v24()
    summary["artifact_identity"][
        "maglev_force_stiffness_equilibrium_energy_finite_difference_generation_identity"
    ].update(
        {
            "equilibrium_maglev_generation": "maglev-100",
            "displacement_maglev_generation": "maglev-99",
            "energy_maglev_generation": "maglev-98",
            "force_maglev_generation": "maglev-97",
            "stiffness_maglev_generation": "maglev-96",
            "coordinate_frame_maglev_generation": "maglev-95",
            "result_equilibrium_displacement_m": 0.0002,
            "result_equilibrium_force_n": 2.0,
            "energy_finite_difference_force_n": [-10.0, 0.0, 10.0],
            "reported_stiffness_n_m": -10000.0,
            "result_force_energy_sign_convention": "force=+dW/dx",
            "result_coordinate_frame": "local_z",
            "result_mesh_sha256": "8" * 64,
            "reported_maglev_result_sha256": "9" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "maglev_force_stiffness_energy_and_equilibrium_share_current_state"
    ]


def test_v24_public_motor_dual_lane_alignment_generation_mismatch() -> None:
    summary = _summary_v24()
    summary["artifact_identity"][
        "motor_dual_lane_geometry_excitation_force_frame_harmonic_alignment_generation_identity"
    ].update(
        {
            "geometry_comparison_generation": "dual-motor-100",
            "excitation_comparison_generation": "dual-motor-99",
            "force_frame_comparison_generation": "dual-motor-98",
            "harmonic_comparison_generation": "dual-motor-97",
            "rotor_angle_comparison_generation": "dual-motor-96",
            "result_geometry_revision_sha256": ["3" * 64, "a" * 64],
            "result_excitation_table_sha256": ["4" * 64, "b" * 64],
            "result_force_coordinate_frames": ["rotor_dq", "global_xyz"],
            "result_harmonic_bins": [[0, 1, 2, 3], [0, 2, 1, 4]],
            "result_rotor_angles_deg": [[0.0, 5.0, 10.0], [10.0, 5.0, 0.0]],
            "result_force_harmonics_n": [
                [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
                [[50.0, 0.0], [1.0, 0.0], [4.0, 2.0], [0.8, 0.2]],
            ],
            "reported_comparison_result_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_dual_lanes_share_geometry_excitation_frame_harmonics_and_angles"
    ]
