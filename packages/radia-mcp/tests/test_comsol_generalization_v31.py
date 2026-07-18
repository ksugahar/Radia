from __future__ import annotations

from test_comsol_generalization_v30 import (
    _summary,
    _with_v30_joule_and_eigenmode_identity,
    gate,
)


_PROMOTED_CASE_IDS = (
    "v31_public_frequency_time_reconstruction_hermitian_window_group_delay_parseval_mismatch",
    "v31_public_rotating_force_virtual_work_stress_tensor_phase_frame_torque_balance_mismatch",
)


def _with_v31_transform_and_force_identity(summary: dict) -> dict:
    summary = _with_v30_joule_and_eigenmode_identity(summary)
    generation = "frequency-time-closure-181"
    summary[
        "frequency_time_hermitian_spacing_window_group_delay_parseval_mesh_result_generation_identity"
    ] = {
        "transform_generation": generation,
        "spectrum_transform_generation": generation,
        "window_transform_generation": generation,
        "delay_transform_generation": generation,
        "energy_transform_generation": generation,
        "mesh_transform_generation": generation,
        "result_transform_generation": generation,
        "frequencies_hz": [0.0, 100.0, 200.0],
        "result_frequencies_hz": [0.0, 100.0, 200.0],
        "frequency_spacing_hz": 100.0,
        "result_frequency_spacing_hz": 100.0,
        "spectrum_ri": [[1.0, 0.0], [0.5, -0.25], [0.2, 0.0]],
        "result_spectrum_ri": [[1.0, 0.0], [0.5, -0.25], [0.2, 0.0]],
        "hermitian_completion": "conjugate_negative_frequencies",
        "result_hermitian_completion": "conjugate_negative_frequencies",
        "window_name": "hann_periodic",
        "result_window_name": "hann_periodic",
        "window_coherent_gain": 0.5,
        "result_window_coherent_gain": 0.5,
        "group_delay_s": 0.001,
        "result_group_delay_s": 0.001,
        "time_origin_s": 0.0,
        "result_time_origin_s": 0.0,
        "frequency_energy": 2.0,
        "time_energy": 2.0,
        "parseval_relative_tolerance": 1.0e-9,
        "transform_mesh_sha256": "1" * 64,
        "result_transform_mesh_sha256": "1" * 64,
        "time_trace_sha256": "2" * 64,
        "accepted_time_trace_sha256": "2" * 64,
    }
    generation = "rotating-force-balance-181"
    summary[
        "rotating_force_virtual_work_stress_phase_frame_lever_angle_power_mesh_result_generation_identity"
    ] = {
        "force_generation": generation,
        "virtual_work_force_generation": generation,
        "stress_force_generation": generation,
        "phase_force_generation": generation,
        "frame_force_generation": generation,
        "angle_force_generation": generation,
        "power_force_generation": generation,
        "mesh_force_generation": generation,
        "result_force_generation": generation,
        "phasor_convention": "exp_positive_jwt_rms",
        "result_phasor_convention": "exp_positive_jwt_rms",
        "coordinate_frame": "rotor_material",
        "result_coordinate_frame": "rotor_material",
        "lever_arm_m": [0.05, 0.0, 0.0],
        "result_lever_arm_m": [0.05, 0.0, 0.0],
        "mechanical_angles_rad": [0.0, 0.1, 0.2],
        "result_mechanical_angles_rad": [0.0, 0.1, 0.2],
        "virtual_work_torque_nm": [1.0, 1.2, 1.1],
        "stress_tensor_torque_nm": [1.0, 1.2, 1.1],
        "torque_relative_tolerance": 1.0e-9,
        "mechanical_power_w": 110.0,
        "airgap_power_w": 110.0,
        "power_relative_tolerance": 1.0e-9,
        "force_mesh_sha256": "3" * 64,
        "result_force_mesh_sha256": "3" * 64,
        "force_result_sha256": "4" * 64,
        "accepted_force_result_sha256": "4" * 64,
    }
    return summary


def test_v31_public_positive_transform_and_force_contracts() -> None:
    result = gate(_with_v31_transform_and_force_identity(_summary()))
    assert result["status"] == "ok"


def test_v31_public_frequency_time_reconstruction_hermitian_window_group_delay_parseval_mismatch() -> None:
    summary = _with_v31_transform_and_force_identity(_summary())
    summary[
        "frequency_time_hermitian_spacing_window_group_delay_parseval_mesh_result_generation_identity"
    ].update(
        {
            "spectrum_transform_generation": "frequency-time-closure-180",
            "delay_transform_generation": "frequency-time-closure-179",
            "result_transform_generation": "frequency-time-closure-178",
            "result_frequencies_hz": [0.0, 90.0, 200.0],
            "result_frequency_spacing_hz": 90.0,
            "result_spectrum_ri": [[1.0, 0.0], [0.5, 0.25]],
            "result_hermitian_completion": "copy_without_conjugation",
            "result_window_name": "rectangular",
            "result_window_coherent_gain": 1.0,
            "result_group_delay_s": 0.0,
            "result_time_origin_s": -0.001,
            "time_energy": 3.0,
            "result_transform_mesh_sha256": "9" * 64,
            "accepted_time_trace_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "frequency_time_reconstruction_uses_current_hermitian_spacing_window_delay_parseval_mesh_and_result"
    ]


def test_v31_public_rotating_force_virtual_work_stress_tensor_phase_frame_torque_balance_mismatch() -> None:
    summary = _with_v31_transform_and_force_identity(_summary())
    summary[
        "rotating_force_virtual_work_stress_phase_frame_lever_angle_power_mesh_result_generation_identity"
    ].update(
        {
            "virtual_work_force_generation": "rotating-force-balance-180",
            "phase_force_generation": "rotating-force-balance-179",
            "result_force_generation": "rotating-force-balance-178",
            "result_phasor_convention": "exp_negative_jwt_peak",
            "result_coordinate_frame": "global_spatial",
            "result_lever_arm_m": [0.0, 0.05, 0.0],
            "result_mechanical_angles_rad": [0.0, 0.2, 0.1],
            "stress_tensor_torque_nm": [-1.0, -1.2, -1.1],
            "airgap_power_w": 95.0,
            "result_force_mesh_sha256": "b" * 64,
            "accepted_force_result_sha256": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "rotating_force_uses_current_virtual_work_stress_phase_frame_angle_power_mesh_and_result"
    ]
