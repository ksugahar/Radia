from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_method_profile_gate import _summary_v22


def _summary_v23():
    summary = _summary_v22()
    identity = summary["artifact_identity"]
    identity["bem_panel_normal_material_region_demag_force_generation_identity"] = {
        "solve_generation": "bem-force-51",
        "panel_mesh_solve_generation": "bem-force-51",
        "outward_normal_solve_generation": "bem-force-51",
        "material_region_solve_generation": "bem-force-51",
        "demag_result_solve_generation": "bem-force-51",
        "force_result_solve_generation": "bem-force-51",
        "panel_ids": [101, 102, 103],
        "result_panel_ids": [101, 102, 103],
        "outward_normals": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "result_outward_normals": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "material_region_ids": [1, 1, 2],
        "result_material_region_ids": [1, 1, 2],
        "demag_field_a_per_m": [-125000.0, -118000.0, -62000.0],
        "result_demag_field_a_per_m": [-125000.0, -118000.0, -62000.0],
        "force_vectors_n": [[12.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 3.0]],
        "result_force_vectors_n": [[12.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 3.0]],
        "force_coordinate_frame": "global_xyz",
        "result_force_coordinate_frame": "global_xyz",
        "panel_force_table_sha256": "1" * 64,
        "result_panel_force_table_sha256": "1" * 64,
    }
    identity["motor_harmonic_rotor_angle_current_phase_force_frame_generation_identity"] = {
        "sweep_generation": "motor-harmonic-51",
        "rotor_angle_sweep_generation": "motor-harmonic-51",
        "current_phase_sweep_generation": "motor-harmonic-51",
        "harmonic_bin_sweep_generation": "motor-harmonic-51",
        "force_frame_sweep_generation": "motor-harmonic-51",
        "force_result_sweep_generation": "motor-harmonic-51",
        "rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "result_rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "current_phase_deg": [0.0, -120.0, 120.0],
        "result_current_phase_deg": [0.0, -120.0, 120.0],
        "harmonic_bins": [0, 1, 2, 3],
        "result_harmonic_bins": [0, 1, 2, 3],
        "force_coordinate_frame": "rotor_dq",
        "result_force_coordinate_frame": "rotor_dq",
        "force_harmonics_n": [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
        "result_force_harmonics_n": [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
        "harmonic_force_table_sha256": "2" * 64,
        "result_harmonic_force_table_sha256": "2" * 64,
    }
    return summary


def test_v23_public_positive_bem_panel_and_motor_harmonic_identity() -> None:
    assert magnetic_force_method_profile_gate(_summary_v23())["status"] == "ok"


def test_v23_public_bem_panel_normal_material_region_demag_force_generation_mismatch() -> None:
    summary = _summary_v23()
    summary["artifact_identity"][
        "bem_panel_normal_material_region_demag_force_generation_identity"
    ].update(
        {
            "panel_mesh_solve_generation": "bem-force-50",
            "outward_normal_solve_generation": "bem-force-49",
            "material_region_solve_generation": "bem-force-48",
            "demag_result_solve_generation": "bem-force-47",
            "force_result_solve_generation": "bem-force-46",
            "result_panel_ids": [101, 103, 104],
            "result_outward_normals": [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
            "result_material_region_ids": [1, 2, 3],
            "result_demag_field_a_per_m": [-90000.0, -62000.0, -30000.0],
            "result_force_vectors_n": [[-8.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
            "result_force_coordinate_frame": "local_xyz",
            "result_panel_force_table_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_demag_force_uses_current_panels_normals_materials_and_results"
    ]


def test_v23_public_motor_harmonic_rotor_angle_current_phase_force_frame_generation_mismatch() -> None:
    summary = _summary_v23()
    summary["artifact_identity"][
        "motor_harmonic_rotor_angle_current_phase_force_frame_generation_identity"
    ].update(
        {
            "rotor_angle_sweep_generation": "motor-harmonic-50",
            "current_phase_sweep_generation": "motor-harmonic-49",
            "harmonic_bin_sweep_generation": "motor-harmonic-48",
            "force_frame_sweep_generation": "motor-harmonic-47",
            "force_result_sweep_generation": "motor-harmonic-46",
            "result_rotor_angles_deg": [15.0, 10.0, 5.0, 0.0],
            "result_current_phase_deg": [0.0, 120.0, -120.0],
            "result_harmonic_bins": [0, 2, 1, 4],
            "result_force_coordinate_frame": "global_xyz",
            "result_force_harmonics_n": [[50.0, 0.0], [1.0, 0.0], [4.0, 2.0], [0.8, 0.2]],
            "result_harmonic_force_table_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_force_harmonics_use_current_angles_phases_bins_and_frame"
    ]
