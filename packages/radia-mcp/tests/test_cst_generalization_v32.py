from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v31 import _summary_v31


_PROMOTED_CASE_IDS = (
    "v32_public_dispersive_port_mode_branch_cutoff_normalization_group_delay_mismatch",
    "v32_public_transient_farfield_time_gate_fft_window_phase_center_energy_mismatch",
)


def _summary_v32():
    summary = _summary_v31()
    for index, row in enumerate(summary["runs"]):
        generation = f"dispersive-port-{371 + index}"
        row[
            "dispersive_port_mode_branch_cutoff_normalization_beta_phase_group_delay_mesh_result_identity"
        ] = {
            "port_generation": generation,
            **{
                key: generation
                for key in (
                    "mode_port_generation",
                    "branch_port_generation",
                    "cutoff_port_generation",
                    "normalization_port_generation",
                    "beta_port_generation",
                    "phase_port_generation",
                    "delay_port_generation",
                    "mesh_port_generation",
                    "result_port_generation",
                )
            },
            "mode_id": "port1:TE10",
            "result_mode_id": "port1:TE10",
            "tracked_branch_id": "branch:TE10:forward",
            "result_tracked_branch_id": "branch:TE10:forward",
            "cutoff_frequency_hz": 6.56e9,
            "result_cutoff_frequency_hz": 6.56e9,
            "frequency_hz": [8.0e9, 9.0e9, 10.0e9],
            "result_frequency_hz": [8.0e9, 9.0e9, 10.0e9],
            "modal_normalization": "unit_forward_power",
            "result_modal_normalization": "unit_forward_power",
            "propagation_constant_sign": "positive_forward",
            "result_propagation_constant_sign": "positive_forward",
            "propagation_constant_rad_per_m": [101.0, 132.0, 158.0],
            "result_propagation_constant_rad_per_m": [101.0, 132.0, 158.0],
            "deembedded_phase_rad": [0.4, 0.2, 0.0],
            "result_deembedded_phase_rad": [0.4, 0.2, 0.0],
            "group_delay_s": 3.183098861837907e-11,
            "result_group_delay_s": 3.183098861837907e-11,
            "mesh_sha256": "1" * 64,
            "result_mesh_sha256": "1" * 64,
            "result_owner": "waveguide/case-371/port1-te10",
            "accepted_result_owner": "waveguide/case-371/port1-te10",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        }
        generation = f"transient-farfield-{371 + index}"
        row[
            "transient_farfield_time_gate_fft_window_phase_center_angular_energy_monitor_result_identity"
        ] = {
            "farfield_generation": generation,
            **{
                key: generation
                for key in (
                    "gate_farfield_generation",
                    "fft_farfield_generation",
                    "phase_center_farfield_generation",
                    "angular_farfield_generation",
                    "energy_farfield_generation",
                    "monitor_farfield_generation",
                    "owner_farfield_generation",
                    "result_farfield_generation",
                )
            },
            "time_gate_s": [2.0e-9, 8.0e-9],
            "result_time_gate_s": [2.0e-9, 8.0e-9],
            "fft_window": "hann",
            "result_fft_window": "hann",
            "fft_normalization": "one_sided_energy_preserving",
            "result_fft_normalization": "one_sided_energy_preserving",
            "phase_center_m": [0.0, 0.0, 0.0],
            "result_phase_center_m": [0.0, 0.0, 0.0],
            "theta_deg": [0.0, 45.0, 90.0],
            "result_theta_deg": [0.0, 45.0, 90.0],
            "phi_deg": [0.0, 90.0, 180.0, 270.0],
            "result_phi_deg": [0.0, 90.0, 180.0, 270.0],
            "accepted_energy_j": 1.0,
            "result_accepted_energy_j": 1.0,
            "radiated_energy_j": 0.82,
            "result_radiated_energy_j": 0.82,
            "monitor_owner": "project-371:monitor-transient-ff",
            "result_monitor_owner": "project-371:monitor-transient-ff",
            "result_owner": "farfield/case-371/time-domain",
            "accepted_result_owner": "farfield/case-371/time-domain",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        }
    return summary


def test_v32_public_positive_dispersive_port_and_transient_farfield():
    assert nonlinear_inductance_sweep_gate(_summary_v32())["status"] == "ok"


def test_v32_public_dispersive_port_mode_branch_cutoff_normalization_group_delay_mismatch():
    summary = _summary_v32()
    identity = summary["runs"][0][
        "dispersive_port_mode_branch_cutoff_normalization_beta_phase_group_delay_mesh_result_identity"
    ]
    identity.update(
        {
            "branch_port_generation": "dispersive-port-370",
            "mesh_port_generation": "dispersive-port-369",
            "result_port_generation": "dispersive-port-368",
            "result_mode_id": "port1:TM01",
            "result_tracked_branch_id": "branch:TM01:backward",
            "result_cutoff_frequency_hz": 7.2e9,
            "result_modal_normalization": "unit_voltage",
            "result_propagation_constant_sign": "negative_forward",
            "result_propagation_constant_rad_per_m": [-101.0, -132.0, -158.0],
            "result_deembedded_phase_rad": [0.0, 0.2, 0.4],
            "result_group_delay_s": -3.183098861837907e-11,
            "result_mesh_sha256": "8" * 64,
            "accepted_result_owner": "waveguide/old-port",
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "dispersive_ports_use_current_mode_branch_cutoff_power_normalization_beta_phase_group_delay_mesh_and_result"
    ]


def test_v32_public_transient_farfield_time_gate_fft_window_phase_center_energy_mismatch():
    summary = _summary_v32()
    identity = summary["runs"][0][
        "transient_farfield_time_gate_fft_window_phase_center_angular_energy_monitor_result_identity"
    ]
    identity.update(
        {
            "gate_farfield_generation": "transient-farfield-370",
            "monitor_farfield_generation": "transient-farfield-369",
            "result_farfield_generation": "transient-farfield-368",
            "result_time_gate_s": [0.0, 4.0e-9],
            "result_fft_window": "rectangular",
            "result_fft_normalization": "raw_fft",
            "result_phase_center_m": [0.01, 0.0, 0.0],
            "result_theta_deg": [90.0, 45.0, 0.0],
            "result_phi_deg": [0.0, 180.0],
            "result_accepted_energy_j": 0.5,
            "result_radiated_energy_j": 1.2,
            "result_monitor_owner": "project-old:monitor-old",
            "accepted_result_owner": "farfield/old-result",
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "transient_farfields_use_current_time_gate_fft_phase_center_angles_energy_monitor_owner_and_result"
    ]
