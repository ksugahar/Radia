from __future__ import annotations

from test_femm_generalization_v19 import _gate
from test_femm_generalization_v22 import _identity_v22
from test_force_coenergy_gate import _quadratic_case


def _identity_v23(sample_count):
    identity = _identity_v22(sample_count)
    identity["nonlinear_bh_branch_operating_point_force_mesh_generation_identity"] = {
        "solve_generation": "nonlinear-force-51",
        "branch_solve_generation": "nonlinear-force-51",
        "operating_point_solve_generation": "nonlinear-force-51",
        "permeability_solve_generation": "nonlinear-force-51",
        "force_mesh_solve_generation": "nonlinear-force-51",
        "force_result_solve_generation": "nonlinear-force-51",
        "branch_id": "ascending:23",
        "force_branch_id": "ascending:23",
        "operating_point_current_a": 7.5,
        "force_operating_point_current_a": 7.5,
        "operating_point_flux_density_t": [1.35, 1.42],
        "force_operating_point_flux_density_t": [1.35, 1.42],
        "permeability_state_sha256": "1" * 64,
        "force_permeability_state_sha256": "1" * 64,
        "force_mesh_sha256": "2" * 64,
        "integrated_force_mesh_sha256": "2" * 64,
        "force_n": [18.2, -0.4],
        "reported_force_n": [18.2, -0.4],
    }
    identity["sliding_band_angle_mesh_harmonic_torque_generation_identity"] = {
        "sweep_generation": "sliding-sweep-51",
        "angle_sweep_generation": "sliding-sweep-51",
        "airgap_mesh_sweep_generation": "sliding-sweep-51",
        "phase_current_sweep_generation": "sliding-sweep-51",
        "torque_sample_sweep_generation": "sliding-sweep-51",
        "harmonic_sweep_generation": "sliding-sweep-51",
        "rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "torque_rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "airgap_mesh_sha256": "3" * 64,
        "torque_airgap_mesh_sha256": "3" * 64,
        "phase_current_table_sha256": "4" * 64,
        "torque_phase_current_table_sha256": "4" * 64,
        "torque_samples_nm": [1.0, 1.2, 0.9, 1.1],
        "harmonic_torque_samples_nm": [1.0, 1.2, 0.9, 1.1],
        "harmonic_orders": [0, 1, 2],
        "reported_harmonic_orders": [0, 1, 2],
        "harmonic_amplitudes_nm": [1.05, 0.12, 0.04],
        "reported_harmonic_amplitudes_nm": [1.05, 0.12, 0.04],
        "torque_sample_table_sha256": "5" * 64,
        "harmonic_sample_table_sha256": "5" * 64,
    }
    return identity


def test_v23_public_positive_nonlinear_force_and_sliding_harmonic_identity():
    positions, _, _ = _quadratic_case()
    assert _gate(_identity_v23(len(positions)))["status"] == "ok"


def test_v23_public_nonlinear_bh_branch_operating_point_force_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v23(len(positions))
    identity[
        "nonlinear_bh_branch_operating_point_force_mesh_generation_identity"
    ].update(
        {
            "branch_solve_generation": "nonlinear-force-50",
            "operating_point_solve_generation": "nonlinear-force-49",
            "permeability_solve_generation": "nonlinear-force-48",
            "force_mesh_solve_generation": "nonlinear-force-47",
            "force_result_solve_generation": "nonlinear-force-46",
            "force_branch_id": "descending:23",
            "force_operating_point_current_a": 6.0,
            "force_operating_point_flux_density_t": [1.1, 1.2],
            "force_permeability_state_sha256": "b" * 64,
            "integrated_force_mesh_sha256": "c" * 64,
            "reported_force_n": [15.0, 0.7],
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_force_uses_current_bh_branch_operating_point_mu_and_mesh"
    ]


def test_v23_public_sliding_band_angle_mesh_harmonic_torque_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v23(len(positions))
    identity["sliding_band_angle_mesh_harmonic_torque_generation_identity"].update(
        {
            "angle_sweep_generation": "sliding-sweep-50",
            "airgap_mesh_sweep_generation": "sliding-sweep-49",
            "phase_current_sweep_generation": "sliding-sweep-48",
            "torque_sample_sweep_generation": "sliding-sweep-47",
            "harmonic_sweep_generation": "sliding-sweep-46",
            "torque_rotor_angles_deg": [0.0, 10.0, 5.0, 15.0],
            "torque_airgap_mesh_sha256": "d" * 64,
            "torque_phase_current_table_sha256": "e" * 64,
            "harmonic_torque_samples_nm": [1.0, 0.9, 1.2, 1.1],
            "reported_harmonic_orders": [0, 2, 1],
            "reported_harmonic_amplitudes_nm": [1.05, 0.04, 0.12],
            "harmonic_sample_table_sha256": "f" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "sliding_band_harmonics_use_current_angles_mesh_currents_and_samples"
    ]
