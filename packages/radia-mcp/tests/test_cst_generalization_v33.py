from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v32 import _summary_v32


_PROMOTED_CASE_IDS = (
    "v33_public_eigenmode_q_energy_conductor_dielectric_radiation_inverse_sum_branch_mismatch",
    "v33_public_tdr_reference_plane_velocity_time_zero_impedance_reflection_causality_energy_mismatch",
)


def _summary_v33():
    summary = _summary_v32()
    for index, row in enumerate(summary["runs"]):
        generation = f"eigenmode-q-{381 + index}"
        q_conductor = 10000.0
        q_dielectric = 20000.0
        q_radiation = 50000.0
        q_total = 1.0 / (
            1.0 / q_conductor + 1.0 / q_dielectric + 1.0 / q_radiation
        )
        row[
            "eigenmode_frequency_branch_energy_conductor_dielectric_radiation_q_mesh_owner_result_identity"
        ] = {
            "eigenmode_generation": generation,
            **{
                key: generation
                for key in (
                    "frequency_generation",
                    "branch_generation",
                    "energy_generation",
                    "conductor_q_generation",
                    "dielectric_q_generation",
                    "radiation_q_generation",
                    "inverse_sum_generation",
                    "mesh_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "mode_id": "mode-1",
            "result_mode_id": "mode-1",
            "mode_branch": "fundamental",
            "result_mode_branch": "fundamental",
            "frequency_hz": 10.0e9,
            "result_frequency_hz": 10.0e9,
            "electric_energy_j": 0.5,
            "result_electric_energy_j": 0.5,
            "magnetic_energy_j": 0.5,
            "result_magnetic_energy_j": 0.5,
            "stored_energy_j": 1.0,
            "result_stored_energy_j": 1.0,
            "q_conductor": q_conductor,
            "result_q_conductor": q_conductor,
            "q_dielectric": q_dielectric,
            "result_q_dielectric": q_dielectric,
            "q_radiation": q_radiation,
            "result_q_radiation": q_radiation,
            "q_total": q_total,
            "result_q_total": q_total,
            "mesh_sha256": "1" * 64,
            "result_mesh_sha256": "1" * 64,
            "mode_owner": "cavity/case-381/mode-1",
            "accepted_mode_owner": "cavity/case-381/mode-1",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        }
        generation = f"tdr-reference-{381 + index}"
        times = [0.0, 1.0e-9, 2.0e-9, 3.0e-9, 4.0e-9]
        waveform = [0.0, 0.0, 0.0, 0.1, 0.05]
        row[
            "tdr_reference_plane_velocity_time_zero_impedance_arrival_window_causality_energy_owner_result_identity"
        ] = {
            "tdr_generation": generation,
            **{
                key: generation
                for key in (
                    "reference_generation",
                    "velocity_generation",
                    "time_zero_generation",
                    "impedance_generation",
                    "arrival_generation",
                    "window_generation",
                    "causality_generation",
                    "energy_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "reference_plane_m": 0.0,
            "result_reference_plane_m": 0.0,
            "propagation_velocity_m_per_s": 2.0e8,
            "result_propagation_velocity_m_per_s": 2.0e8,
            "time_zero_s": 1.0e-9,
            "result_time_zero_s": 1.0e-9,
            "characteristic_impedance_ohm": 50.0,
            "result_characteristic_impedance_ohm": 50.0,
            "reflection_distance_m": 0.2,
            "result_reflection_distance_m": 0.2,
            "reflection_arrival_s": 3.0e-9,
            "result_reflection_arrival_s": 3.0e-9,
            "time_window_s": [0.0, 4.0e-9],
            "result_time_window_s": [0.0, 4.0e-9],
            "time_samples_s": times,
            "result_time_samples_s": list(times),
            "reflection_waveform": waveform,
            "result_reflection_waveform": list(waveform),
            "pre_arrival_max_abs": 0.0,
            "result_pre_arrival_max_abs": 0.0,
            "incident_energy_j": 1.0,
            "result_incident_energy_j": 1.0,
            "reflected_energy_j": 0.1,
            "result_reflected_energy_j": 0.1,
            "accepted_energy_j": 0.9,
            "result_accepted_energy_j": 0.9,
            "waveform_owner": "tdr/case-381/port-1",
            "accepted_waveform_owner": "tdr/case-381/port-1",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        }
    return summary


def test_v33_public_positive_eigenmode_q_and_tdr_closure():
    assert nonlinear_inductance_sweep_gate(_summary_v33())["status"] == "ok"


def test_v33_public_eigenmode_q_energy_conductor_dielectric_radiation_inverse_sum_branch_mismatch():
    summary = _summary_v33()
    identity = summary["runs"][0][
        "eigenmode_frequency_branch_energy_conductor_dielectric_radiation_q_mesh_owner_result_identity"
    ]
    identity.update(
        {
            "branch_generation": "eigenmode-q-380",
            "energy_generation": "eigenmode-q-379",
            "result_generation": "eigenmode-q-378",
            "result_mode_id": "mode-2",
            "result_mode_branch": "stale-crossing-branch",
            "result_frequency_hz": 9.0e9,
            "result_electric_energy_j": 0.1,
            "result_magnetic_energy_j": 0.2,
            "result_stored_energy_j": 2.0,
            "result_q_conductor": 1000.0,
            "result_q_dielectric": 2000.0,
            "result_q_radiation": -5000.0,
            "result_q_total": 10000.0,
            "result_mesh_sha256": "8" * 64,
            "accepted_mode_owner": "cavity/old-mode",
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "eigenmodes_use_current_frequency_branch_energy_q_inverse_sum_mesh_owner_and_result"
    ]


def test_v33_public_tdr_reference_plane_velocity_time_zero_impedance_reflection_causality_energy_mismatch():
    summary = _summary_v33()
    identity = summary["runs"][0][
        "tdr_reference_plane_velocity_time_zero_impedance_arrival_window_causality_energy_owner_result_identity"
    ]
    identity.update(
        {
            "reference_generation": "tdr-reference-380",
            "arrival_generation": "tdr-reference-379",
            "result_generation": "tdr-reference-378",
            "result_reference_plane_m": 0.1,
            "result_propagation_velocity_m_per_s": 3.0e8,
            "result_time_zero_s": -1.0e-9,
            "result_characteristic_impedance_ohm": 75.0,
            "result_reflection_distance_m": 0.5,
            "result_reflection_arrival_s": 1.0e-9,
            "result_time_window_s": [2.0e-9, 1.0e-9],
            "result_time_samples_s": [4.0e-9, 3.0e-9, 2.0e-9],
            "result_reflection_waveform": [0.2, 0.1, 0.0],
            "result_pre_arrival_max_abs": 0.2,
            "result_incident_energy_j": 0.5,
            "result_reflected_energy_j": 0.8,
            "result_accepted_energy_j": -0.3,
            "accepted_waveform_owner": "tdr/old-port",
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "tdr_uses_current_reference_velocity_time_zero_impedance_arrival_causality_energy_owner_and_result"
    ]
