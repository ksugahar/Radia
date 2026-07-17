from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v24 import _summary_v24


_PROMOTED_CASE_IDS = (
    "v25_public_adaptive_mesh_pass_sparameter_energy_convergence_grid_generation_mismatch",
    "v25_public_eigenmode_tracking_phase_normalization_port_coupling_mesh_mismatch",
)


def _summary_v25():
    summary = _summary_v24()
    for index, row in enumerate(summary["runs"]):
        generation = f"adaptive-pass-{201 + index}"
        row[
            "adaptive_mesh_pass_sparameter_energy_convergence_grid_generation_identity"
        ] = {
            "adaptive_generation": generation,
            "mesh_pass_adaptive_generation": generation,
            "sparameter_adaptive_generation": generation,
            "energy_adaptive_generation": generation,
            "frequency_grid_adaptive_generation": generation,
            "stopping_rule_adaptive_generation": generation,
            "result_adaptive_generation": generation,
            "mesh_pass_ids": [0, 1, 2],
            "result_mesh_pass_ids": [0, 1, 2],
            "mesh_cell_counts": [10000, 18000, 29000],
            "result_mesh_cell_counts": [10000, 18000, 29000],
            "frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "result_frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "maximum_sparameter_delta": [0.1, 0.03, 0.005],
            "result_maximum_sparameter_delta": [0.1, 0.03, 0.005],
            "stored_energy_closure_residual": [0.05, 0.01, 0.001],
            "result_stored_energy_closure_residual": [0.05, 0.01, 0.001],
            "sparameter_delta_tolerance": 0.01,
            "energy_closure_tolerance": 0.005,
            "converged_pass_id": 2,
            "result_converged_pass_id": 2,
            "adaptive_result_sha256": "1" * 64,
            "reported_adaptive_result_sha256": "1" * 64,
        }
        generation = f"eigenmode-track-{201 + index}"
        row[
            "eigenmode_tracking_phase_normalization_port_coupling_mesh_generation_identity"
        ] = {
            "tracking_generation": generation,
            "modal_subspace_tracking_generation": generation,
            "phase_tracking_generation": generation,
            "normalization_tracking_generation": generation,
            "port_coupling_tracking_generation": generation,
            "mesh_tracking_generation": generation,
            "result_tracking_generation": generation,
            "sweep_parameters": [0.0, 0.5, 1.0],
            "result_sweep_parameters": [0.0, 0.5, 1.0],
            "tracked_mode_ids": ["mode-1", "mode-2"],
            "result_tracked_mode_ids": ["mode-1", "mode-2"],
            "modal_subspace_sha256": ["2" * 64, "3" * 64, "4" * 64],
            "result_modal_subspace_sha256": ["2" * 64, "3" * 64, "4" * 64],
            "phase_anchor_ids": ["probe-ez", "probe-hy"],
            "result_phase_anchor_ids": ["probe-ez", "probe-hy"],
            "normalization": "stored_energy_1j",
            "result_normalization": "stored_energy_1j",
            "port_coupling_magnitudes": [[0.8, 0.1], [0.78, 0.12], [0.75, 0.15]],
            "result_port_coupling_magnitudes": [[0.8, 0.1], [0.78, 0.12], [0.75, 0.15]],
            "mesh_sha256": ["5" * 64, "6" * 64, "7" * 64],
            "result_mesh_sha256": ["5" * 64, "6" * 64, "7" * 64],
            "eigenmode_track_sha256": "8" * 64,
            "reported_eigenmode_track_sha256": "8" * 64,
        }
    return summary


def test_v25_public_positive_adaptive_and_eigenmode_identity() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v25())["status"] == "ok"


def test_v25_public_adaptive_mesh_convergence_generation_mismatch() -> None:
    summary = _summary_v25()
    identity = summary["runs"][0][
        "adaptive_mesh_pass_sparameter_energy_convergence_grid_generation_identity"
    ]
    identity.update(
        {
            "mesh_pass_adaptive_generation": "adaptive-pass-200",
            "sparameter_adaptive_generation": "adaptive-pass-199",
            "energy_adaptive_generation": "adaptive-pass-198",
            "frequency_grid_adaptive_generation": "adaptive-pass-197",
            "stopping_rule_adaptive_generation": "adaptive-pass-196",
            "result_mesh_pass_ids": [0, 2, 1],
            "result_mesh_cell_counts": [10000, 29000, 18000],
            "result_frequency_grid_hz": [1.0e9, 1.6e9, 2.0e9],
            "result_maximum_sparameter_delta": [0.1, 0.03, 0.02],
            "result_stored_energy_closure_residual": [0.05, 0.01, 0.02],
            "result_converged_pass_id": 1,
            "reported_adaptive_result_sha256": "d" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "adaptive_results_use_current_mesh_pass_sparameter_energy_grid_and_stop_rule"
    ]


def test_v25_public_eigenmode_tracking_generation_mismatch() -> None:
    summary = _summary_v25()
    identity = summary["runs"][0][
        "eigenmode_tracking_phase_normalization_port_coupling_mesh_generation_identity"
    ]
    identity.update(
        {
            "modal_subspace_tracking_generation": "eigenmode-track-200",
            "phase_tracking_generation": "eigenmode-track-199",
            "normalization_tracking_generation": "eigenmode-track-198",
            "port_coupling_tracking_generation": "eigenmode-track-197",
            "mesh_tracking_generation": "eigenmode-track-196",
            "result_sweep_parameters": [1.0, 0.5, 0.0],
            "result_tracked_mode_ids": ["mode-2", "mode-1"],
            "result_modal_subspace_sha256": ["4" * 64, "3" * 64, "2" * 64],
            "result_phase_anchor_ids": ["probe-hy", "probe-ez"],
            "result_normalization": "peak_field_1",
            "result_port_coupling_magnitudes": [[0.1, 0.8]],
            "result_mesh_sha256": ["7" * 64, "6" * 64, "5" * 64],
            "reported_eigenmode_track_sha256": "e" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "eigenmodes_use_current_subspace_phase_normalization_ports_and_mesh"
    ]
