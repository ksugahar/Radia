from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import nonlinear_inductance_sweep_gate
from test_cst_generalization_v29 import _summary_v29

_PROMOTED_CASE_IDS = (
    "v30_public_waveguide_port_mode_power_normalization_deembed_plane_smatrix_reference_mismatch",
    "v30_public_sar_mass_average_density_voxel_frequency_field_tissue_result_mismatch",
)

def _summary_v30():
    summary = _summary_v29()
    for index, row in enumerate(summary["runs"]):
        generation = f"waveguide-port-{351 + index}"
        row["waveguide_port_mode_power_deembed_impedance_frequency_port_smatrix_result_identity"] = {
            "port_generation": generation, "mode_port_generation": generation, "power_port_generation": generation,
            "deembed_port_generation": generation, "impedance_port_generation": generation,
            "frequency_port_generation": generation, "order_port_generation": generation, "result_port_generation": generation,
            "mode_ids": ["port1:TE10", "port2:TE10"], "result_mode_ids": ["port1:TE10", "port2:TE10"],
            "power_normalization_w": [1.0, 1.0], "result_power_normalization_w": [1.0, 1.0],
            "deembed_plane_m": [0.0, 0.1], "result_deembed_plane_m": [0.0, 0.1],
            "reference_impedance_ohm": [50.0, 50.0], "result_reference_impedance_ohm": [50.0, 50.0],
            "frequency_hz": [8.0e9, 9.0e9, 10.0e9], "result_frequency_hz": [8.0e9, 9.0e9, 10.0e9],
            "port_order": [1, 2], "result_port_order": [1, 2],
            "smatrix_ri": [[[0.1, 0.0], [0.9, 0.0]], [[0.9, 0.0], [0.1, 0.0]]],
            "result_smatrix_ri": [[[0.1, 0.0], [0.9, 0.0]], [[0.9, 0.0], [0.1, 0.0]]],
            "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        }
        generation = f"sar-mass-{351 + index}"
        row["sar_mass_density_voxel_frequency_field_mesh_result_identity"] = {
            "sar_generation": generation, "mass_sar_generation": generation, "density_sar_generation": generation,
            "voxel_sar_generation": generation, "frequency_sar_generation": generation,
            "field_sar_generation": generation, "mesh_sar_generation": generation, "result_sar_generation": generation,
            "averaging_mass_kg": 0.01, "result_averaging_mass_kg": 0.01,
            "tissue_density_kg_m3": 1000.0, "result_tissue_density_kg_m3": 1000.0,
            "voxel_ids": [101, 102, 103], "result_voxel_ids": [101, 102, 103],
            "voxel_mass_kg": [0.003, 0.004, 0.003], "result_voxel_mass_kg": [0.003, 0.004, 0.003],
            "frequency_hz": 2.45e9, "result_frequency_hz": 2.45e9,
            "field_normalization": "accepted_power_1w", "result_field_normalization": "accepted_power_1w",
            "sar_w_kg": 1.2, "result_sar_w_kg": 1.2,
            "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        }
    return summary

def test_v30_public_positive_waveguide_and_sar_identities():
    assert nonlinear_inductance_sweep_gate(_summary_v30())["status"] == "ok"

def test_v30_public_waveguide_port_mode_power_normalization_deembed_plane_smatrix_reference_mismatch():
    summary = _summary_v30(); identity = summary["runs"][0]["waveguide_port_mode_power_deembed_impedance_frequency_port_smatrix_result_identity"]
    identity.update({
        "mode_port_generation": "waveguide-port-350", "frequency_port_generation": "waveguide-port-349",
        "result_mode_ids": ["port1:TM01"], "result_power_normalization_w": [0.5, 2.0],
        "result_deembed_plane_m": [0.02, 0.08], "result_reference_impedance_ohm": [75.0, 50.0],
        "result_frequency_hz": [8.5e9, 9.5e9], "result_port_order": [2, 1],
        "result_smatrix_ri": [[[1.2, 0.0]]], "result_mesh_sha256": "7" * 64,
        "accepted_result_sha256": "8" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary); assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["waveguide_ports_use_current_modes_power_deembed_impedance_frequency_order_smatrix_mesh_and_result"]

def test_v30_public_sar_mass_average_density_voxel_frequency_field_tissue_result_mismatch():
    summary = _summary_v30(); identity = summary["runs"][0]["sar_mass_density_voxel_frequency_field_mesh_result_identity"]
    identity.update({
        "mass_sar_generation": "sar-mass-350", "voxel_sar_generation": "sar-mass-349",
        "result_averaging_mass_kg": 0.001, "result_tissue_density_kg_m3": 800.0,
        "result_voxel_ids": [103, 101], "result_voxel_mass_kg": [0.01, 0.01],
        "result_frequency_hz": 5.8e9, "result_field_normalization": "peak_field_1v_m",
        "result_sar_w_kg": -1.0, "result_mesh_sha256": "9" * 64,
        "accepted_result_sha256": "a" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary); assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["sar_uses_current_average_mass_density_voxels_frequency_field_mesh_and_result"]
