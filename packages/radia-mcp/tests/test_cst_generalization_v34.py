from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v33 import _summary_v33


_PROMOTED_CASE_IDS = (
    "v34_public_sparameter_reference_plane_time_gate_passivity_causality_energy_mismatch",
    "v34_public_eigenmode_degeneracy_subspace_orthogonality_tracking_mesh_convergence_mismatch",
)


def _summary_v34():
    summary = _summary_v33()
    for index, row in enumerate(summary["runs"]):
        generation = f"sparameter-gated-{391 + index}"
        row[
            "sparameter_reference_plane_time_gate_causality_passivity_energy_port_frequency_owner_result_identity"
        ] = {
            "sparameter_generation": generation,
            **{
                key: generation
                for key in (
                    "reference_generation",
                    "gate_generation",
                    "causality_generation",
                    "passivity_generation",
                    "energy_generation",
                    "port_generation",
                    "frequency_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "reference_plane_shift_m": 0.01,
            "result_reference_plane_shift_m": 0.01,
            "time_gate_window_s": [0.0, 4.0e-9],
            "result_time_gate_window_s": [0.0, 4.0e-9],
            "impulse_time_s": [-1.0e-9, 0.0, 1.0e-9, 2.0e-9, 3.0e-9],
            "result_impulse_time_s": [-1.0e-9, 0.0, 1.0e-9, 2.0e-9, 3.0e-9],
            "impulse_response": [0.0, 0.2, 0.1, 0.05, 0.0],
            "result_impulse_response": [0.0, 0.2, 0.1, 0.05, 0.0],
            "pre_zero_max_abs": 0.0,
            "result_pre_zero_max_abs": 0.0,
            "maximum_singular_values": [0.8, 0.9, 0.85],
            "result_maximum_singular_values": [0.8, 0.9, 0.85],
            "incident_energy_j": 1.0,
            "result_incident_energy_j": 1.0,
            "reflected_energy_j": 0.1,
            "result_reflected_energy_j": 0.1,
            "transmitted_energy_j": 0.8,
            "result_transmitted_energy_j": 0.8,
            "absorbed_energy_j": 0.1,
            "result_absorbed_energy_j": 0.1,
            "port_impedance_ohm": [50.0, 50.0],
            "result_port_impedance_ohm": [50.0, 50.0],
            "frequency_grid_hz": [1.0e9, 2.0e9, 3.0e9],
            "result_frequency_grid_hz": [1.0e9, 2.0e9, 3.0e9],
            "sparameter_owner": "network/case-391/gated",
            "accepted_sparameter_owner": "network/case-391/gated",
            "sparameter_sha256": "1" * 64,
            "accepted_sparameter_sha256": "1" * 64,
        }

        generation = f"degenerate-eigenmode-{391 + index}"
        row[
            "eigenmode_degenerate_subspace_principal_angle_mass_orthogonality_phase_tracking_residual_mesh_owner_result_identity"
        ] = {
            "degenerate_mode_generation": generation,
            **{
                key: generation
                for key in (
                    "subspace_generation",
                    "principal_angle_generation",
                    "mass_generation",
                    "phase_generation",
                    "tracking_generation",
                    "residual_generation",
                    "mesh_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "mode_frequencies_hz": [10.0e9, 10.0e9 + 1000.0],
            "result_mode_frequencies_hz": [10.0e9, 10.0e9 + 1000.0],
            "principal_angles_rad": [0.0, 0.001],
            "result_principal_angles_rad": [0.0, 0.001],
            "mass_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "result_mass_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "mass_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "result_mass_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "phase_anchor_complex": [[1.0, 0.0], [1.0, 0.0]],
            "result_phase_anchor_complex": [[1.0, 0.0], [1.0, 0.0]],
            "tracking_subspace_ids": ["subspace-391", "subspace-391"],
            "result_tracking_subspace_ids": ["subspace-391", "subspace-391"],
            "residual_norms": [1.0e-9, 2.0e-9],
            "result_residual_norms": [1.0e-9, 2.0e-9],
            "mesh_cell_counts": [1000, 8000, 64000],
            "result_mesh_cell_counts": [1000, 8000, 64000],
            "mesh_converged_frequency_hz": [9.8e9, 9.98e9, 10.0e9],
            "result_mesh_converged_frequency_hz": [9.8e9, 9.98e9, 10.0e9],
            "eigenmode_mesh_sha256": "2" * 64,
            "result_eigenmode_mesh_sha256": "2" * 64,
            "field_owner": "eigenmode/case-391/subspace",
            "accepted_field_owner": "eigenmode/case-391/subspace",
            "field_sha256": "3" * 64,
            "accepted_field_sha256": "3" * 64,
        }
    return summary


def test_v34_public_positive_gated_sparameters_and_degenerate_subspace():
    assert nonlinear_inductance_sweep_gate(_summary_v34())["status"] == "ok"


def test_v34_public_sparameter_reference_plane_time_gate_passivity_causality_energy_mismatch():
    summary = _summary_v34()
    identity = summary["runs"][0][
        "sparameter_reference_plane_time_gate_causality_passivity_energy_port_frequency_owner_result_identity"
    ]
    identity.update(
        {
            "reference_generation": "sparameter-gated-390",
            "energy_generation": "sparameter-gated-389",
            "result_generation": "sparameter-gated-388",
            "result_reference_plane_shift_m": -0.02,
            "result_time_gate_window_s": [5.0e-9, 1.0e-9],
            "result_impulse_time_s": [2.0e-9, 1.0e-9, -1.0e-9],
            "result_impulse_response": [0.1, 0.2, 0.5],
            "result_pre_zero_max_abs": 0.5,
            "result_maximum_singular_values": [1.2, 1.5],
            "result_incident_energy_j": 0.5,
            "result_reflected_energy_j": 0.8,
            "result_transmitted_energy_j": 0.6,
            "result_absorbed_energy_j": -0.9,
            "result_port_impedance_ohm": [75.0],
            "result_frequency_grid_hz": [3.0e9, 2.0e9, 1.0e9],
            "accepted_sparameter_owner": "network/old",
            "accepted_sparameter_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "sparameters_use_current_reference_gate_causality_passivity_energy_ports_frequency_owner_and_result"
    ]


def test_v34_public_eigenmode_degeneracy_subspace_orthogonality_tracking_mesh_convergence_mismatch():
    summary = _summary_v34()
    identity = summary["runs"][0][
        "eigenmode_degenerate_subspace_principal_angle_mass_orthogonality_phase_tracking_residual_mesh_owner_result_identity"
    ]
    identity.update(
        {
            "subspace_generation": "degenerate-eigenmode-390",
            "tracking_generation": "degenerate-eigenmode-389",
            "result_generation": "degenerate-eigenmode-388",
            "result_mode_frequencies_hz": [9.0e9, 11.0e9],
            "result_principal_angles_rad": [0.5, 1.0],
            "result_mass_gram_real": [[1.0, 0.5], [0.2, 0.1]],
            "result_mass_gram_imag": [[0.0, 1.0], [-1.0, 0.0]],
            "result_phase_anchor_complex": [[-1.0, 1.0], [0.0, 0.0]],
            "result_tracking_subspace_ids": ["mode-a", "mode-b"],
            "result_residual_norms": [1.0, 2.0],
            "result_mesh_cell_counts": [64000, 8000, 1000],
            "result_mesh_converged_frequency_hz": [10.0e9, 9.0e9, 11.0e9],
            "result_eigenmode_mesh_sha256": "b" * 64,
            "accepted_field_owner": "eigenmode/old",
            "accepted_field_sha256": "c" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "degenerate_eigenmodes_use_current_subspace_angles_mass_orthogonality_phase_tracking_residual_mesh_owner_and_result"
    ]


def test_v34_public_self_consistent_nonpassive_sparameters_are_rejected():
    summary = _summary_v34()
    identity = summary["runs"][0][
        "sparameter_reference_plane_time_gate_causality_passivity_energy_port_frequency_owner_result_identity"
    ]
    identity["maximum_singular_values"] = [0.8, 1.01, 0.85]
    identity["result_maximum_singular_values"] = [0.8, 1.01, 0.85]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v34_public_self_consistent_vector_labels_are_not_subspace_tracking():
    summary = _summary_v34()
    identity = summary["runs"][0][
        "eigenmode_degenerate_subspace_principal_angle_mass_orthogonality_phase_tracking_residual_mesh_owner_result_identity"
    ]
    identity["tracking_subspace_ids"] = ["mode-a", "mode-b"]
    identity["result_tracking_subspace_ids"] = ["mode-a", "mode-b"]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"
