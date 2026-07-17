from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v25 import _summary_v25


_PROMOTED_CASE_IDS = (
    "v26_public_port_deembedding_reference_plane_impedance_mode_normalization_smatrix_mismatch",
    "v26_public_farfield_angular_grid_polarization_coordinate_power_normalization_mesh_mismatch",
)


def _summary_v26():
    summary = _summary_v25()
    for index, row in enumerate(summary["runs"]):
        generation = f"port-network-{301 + index}"
        row[
            "port_deembedding_reference_plane_impedance_mode_normalization_smatrix_generation_identity"
        ] = {
            "network_generation": generation,
            "port_mode_network_generation": generation,
            "deembedding_network_generation": generation,
            "reference_impedance_network_generation": generation,
            "normalization_network_generation": generation,
            "frequency_grid_network_generation": generation,
            "result_network_generation": generation,
            "port_mode_ids": ["P1:M1", "P2:M1"],
            "result_port_mode_ids": ["P1:M1", "P2:M1"],
            "reference_plane_offsets_m": [0.001, 0.0015],
            "result_reference_plane_offsets_m": [0.001, 0.0015],
            "reference_impedance_ri_ohm": [[50.0, 0.0], [50.0, 0.0]],
            "result_reference_impedance_ri_ohm": [[50.0, 0.0], [50.0, 0.0]],
            "wave_normalization": "power_wave",
            "result_wave_normalization": "power_wave",
            "frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "result_frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "smatrix_ri": [[[0.1, -0.01], [0.8, -0.1]], [[0.79, -0.12], [0.11, -0.02]]],
            "result_smatrix_ri": [[[0.1, -0.01], [0.8, -0.1]], [[0.79, -0.12], [0.11, -0.02]]],
            "smatrix_sha256": "1" * 64,
            "reported_smatrix_sha256": "1" * 64,
        }
        generation = f"farfield-{301 + index}"
        row[
            "farfield_angular_grid_polarization_coordinate_power_normalization_mesh_generation_identity"
        ] = {
            "farfield_generation": generation,
            "angular_grid_farfield_generation": generation,
            "polarization_farfield_generation": generation,
            "coordinate_farfield_generation": generation,
            "power_farfield_generation": generation,
            "mesh_farfield_generation": generation,
            "result_farfield_generation": generation,
            "theta_deg": [0.0, 45.0, 90.0],
            "result_theta_deg": [0.0, 45.0, 90.0],
            "phi_deg": [0.0, 90.0],
            "result_phi_deg": [0.0, 90.0],
            "polarization_basis": "ludwig3_co_cross",
            "result_polarization_basis": "ludwig3_co_cross",
            "coordinate_frame": "global_xyz_z_up",
            "result_coordinate_frame": "global_xyz_z_up",
            "radiated_power_w": 0.8,
            "result_radiated_power_w": 0.8,
            "field_normalization": "sqrt_radiated_power",
            "result_field_normalization": "sqrt_radiated_power",
            "mesh_sha256": "2" * 64,
            "result_mesh_sha256": "2" * 64,
            "farfield_sha256": "3" * 64,
            "reported_farfield_sha256": "3" * 64,
        }
    return summary


def test_v26_public_positive_port_network_and_farfield_identity() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v26())["status"] == "ok"


def test_v26_public_rejects_port_network_identity_mismatch() -> None:
    summary = _summary_v26()
    identity = summary["runs"][0][
        "port_deembedding_reference_plane_impedance_mode_normalization_smatrix_generation_identity"
    ]
    identity.update(
        {
            "port_mode_network_generation": "port-network-300",
            "result_port_mode_ids": ["P2:M1", "P1:M2"],
            "result_reference_plane_offsets_m": [0.0, 0.002],
            "result_reference_impedance_ri_ohm": [[75.0, 0.0], [50.0, 5.0]],
            "result_wave_normalization": "voltage_wave",
            "result_frequency_grid_hz": [1.0e9, 1.6e9, 2.0e9],
            "reported_smatrix_sha256": "8" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "sparameters_use_current_port_modes_planes_impedances_normalization_grid_and_result"
    ]


def test_v26_public_rejects_farfield_identity_mismatch() -> None:
    summary = _summary_v26()
    identity = summary["runs"][0][
        "farfield_angular_grid_polarization_coordinate_power_normalization_mesh_generation_identity"
    ]
    identity.update(
        {
            "angular_grid_farfield_generation": "farfield-300",
            "result_theta_deg": [90.0, 45.0, 0.0],
            "result_polarization_basis": "spherical_theta_phi",
            "result_coordinate_frame": "local_wcs_y_up",
            "result_radiated_power_w": 1.2,
            "result_field_normalization": "peak_field",
            "result_mesh_sha256": "9" * 64,
            "reported_farfield_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "farfields_use_current_angular_grid_polarization_coordinates_power_mesh_and_result"
    ]
