from __future__ import annotations

import cmath
import math

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v35 import _summary_v35


_PROMOTED_CASE_IDS = (
    "v36_public_cavity_eigenmode_frequency_q_energy_normalization_orthogonality_mesh_mismatch",
    "v36_public_active_sparameter_embedded_pattern_scan_impedance_power_closure_mismatch",
)


def _summary_v36():
    summary = _summary_v35()
    for index, row in enumerate(summary["runs"]):
        generation = f"cavity-mode-{412 + index}"
        q_external, q_dielectric, q_conductor = 20000.0, 50000.0, 40000.0
        q_total = 1.0 / (1.0 / q_external + 1.0 / q_dielectric + 1.0 / q_conductor)
        row["cavity_eigenmode_frequency_q_energy_orthogonality_degeneracy_mesh_owner_result_identity"] = {
            "cavity_generation": generation,
            **{key: generation for key in (
                "frequency_generation", "q_generation", "energy_generation",
                "orthogonality_generation", "degeneracy_generation", "mesh_generation",
                "owner_generation", "result_generation",
            )},
            "mode_ids": [1, 2], "result_mode_ids": [1, 2],
            "mode_frequency_hz": [10.0e9, 10.5e9], "result_mode_frequency_hz": [10.0e9, 10.5e9],
            "q_external": q_external, "result_q_external": q_external,
            "q_dielectric": q_dielectric, "result_q_dielectric": q_dielectric,
            "q_conductor": q_conductor, "result_q_conductor": q_conductor,
            "q_total": q_total, "result_q_total": q_total,
            "electric_energy_j": [0.5, 0.5], "result_electric_energy_j": [0.5, 0.5],
            "magnetic_energy_j": [0.5, 0.5], "result_magnetic_energy_j": [0.5, 0.5],
            "normalization": "unit_total_energy_1j", "result_normalization": "unit_total_energy_1j",
            "mode_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "result_mode_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "mode_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "result_mode_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "degeneracy_order": [1, 2], "result_degeneracy_order": [1, 2],
            "mesh_dof": [10000, 20000, 40000], "result_mesh_dof": [10000, 20000, 40000],
            "mesh_frequency_hz": [9.9e9, 9.98e9, 10.0e9],
            "result_mesh_frequency_hz": [9.9e9, 9.98e9, 10.0e9],
            "cavity_mesh_sha256": "1" * 64, "result_cavity_mesh_sha256": "1" * 64,
            "cavity_owner": "cavity/case-412", "accepted_cavity_owner": "cavity/case-412",
            "cavity_result_sha256": "2" * 64, "accepted_cavity_result_sha256": "2" * 64,
        }
        generation = f"active-array-{412 + index}"
        excitation = [1.0 + 0.0j, cmath.exp(1j * math.pi / 2.0)]
        matrix = [[0.1 + 0.0j, 0.2 + 0.0j], [0.2 + 0.0j, 0.1 + 0.0j]]
        active = [sum(matrix[r][c] * excitation[c] for c in range(2)) / excitation[r] for r in range(2)]
        pair = lambda value: [value.real, value.imag]
        row["active_sparameter_embedded_pattern_scan_impedance_power_frequency_mesh_owner_result_identity"] = {
            "array_generation": generation,
            **{key: generation for key in (
                "sparameter_generation", "pattern_generation", "scan_generation",
                "impedance_generation", "power_generation", "frequency_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )},
            "frequency_hz": 10.0e9, "result_frequency_hz": 10.0e9,
            "port_impedance_ohm": [50.0, 50.0], "result_port_impedance_ohm": [50.0, 50.0],
            "scan_phase_rad": [0.0, math.pi / 2.0], "result_scan_phase_rad": [0.0, math.pi / 2.0],
            "excitation_complex": [pair(value) for value in excitation],
            "result_excitation_complex": [pair(value) for value in excitation],
            "s_matrix_complex": [[pair(value) for value in item] for item in matrix],
            "result_s_matrix_complex": [[pair(value) for value in item] for item in matrix],
            "active_s_complex": [pair(value) for value in active],
            "result_active_s_complex": [pair(value) for value in active],
            "embedded_pattern_complex": [[[1.0, 0.0], [0.0, 1.0]], [[0.5, 0.0], [0.0, 0.5]]],
            "result_embedded_pattern_complex": [[[1.0, 0.0], [0.0, 1.0]], [[0.5, 0.0], [0.0, 0.5]]],
            "incident_power_w": 2.0, "result_incident_power_w": 2.0,
            "reflected_power_w": 0.2, "result_reflected_power_w": 0.2,
            "accepted_power_w": 1.8, "result_accepted_power_w": 1.8,
            "radiated_power_w": 1.5, "result_radiated_power_w": 1.5,
            "dissipated_power_w": 0.3, "result_dissipated_power_w": 0.3,
            "array_mesh_sha256": "3" * 64, "result_array_mesh_sha256": "3" * 64,
            "array_owner": "array/case-412", "accepted_array_owner": "array/case-412",
            "array_result_sha256": "4" * 64, "accepted_array_result_sha256": "4" * 64,
        }
    return summary


def test_v36_public_positive_cavity_and_active_array_closure():
    assert nonlinear_inductance_sweep_gate(_summary_v36())["status"] == "ok"


def test_v36_public_cavity_eigenmode_frequency_q_energy_normalization_orthogonality_mesh_mismatch():
    summary = _summary_v36()
    record = summary["runs"][0]["cavity_eigenmode_frequency_q_energy_orthogonality_degeneracy_mesh_owner_result_identity"]
    record.update({
        "frequency_generation": "cavity-mode-411", "q_generation": "cavity-mode-410",
        "result_generation": "cavity-mode-409", "result_mode_frequency_hz": [1.0, 1.0],
        "result_q_total": -1.0, "result_electric_energy_j": [2.0],
        "result_magnetic_energy_j": [-1.0], "result_normalization": "peak_field",
        "result_mode_gram_real": [[1.0, 1.0], [1.0, 1.0]],
        "result_degeneracy_order": [2, 1], "result_mesh_dof": [40000, 10000],
        "result_mesh_frequency_hz": [8.0e9, 7.0e9], "result_cavity_mesh_sha256": "a" * 64,
        "accepted_cavity_owner": "cavity/old", "accepted_cavity_result_sha256": "b" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["cavity_eigenmodes_use_current_frequency_q_energy_orthogonality_degeneracy_mesh_owner_and_result"]


def test_v36_public_active_sparameter_embedded_pattern_scan_impedance_power_closure_mismatch():
    summary = _summary_v36()
    record = summary["runs"][0]["active_sparameter_embedded_pattern_scan_impedance_power_frequency_mesh_owner_result_identity"]
    record.update({
        "sparameter_generation": "active-array-411", "power_generation": "active-array-410",
        "result_generation": "active-array-409", "result_frequency_hz": 9.0e9,
        "result_port_impedance_ohm": [-50.0], "result_scan_phase_rad": [math.pi],
        "result_excitation_complex": [[0.0, 0.0]], "result_s_matrix_complex": [[[2.0, 0.0]]],
        "result_active_s_complex": [[9.0, 9.0]], "result_embedded_pattern_complex": [],
        "result_incident_power_w": -2.0, "result_reflected_power_w": 3.0,
        "result_accepted_power_w": -1.0, "result_radiated_power_w": 2.0,
        "result_dissipated_power_w": -3.0, "result_array_mesh_sha256": "c" * 64,
        "accepted_array_owner": "array/old", "accepted_array_result_sha256": "d" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["active_arrays_use_current_sparameters_patterns_scan_impedance_power_frequency_mesh_owner_and_result"]


def test_v36_public_rejects_self_consistent_wrong_cavity_q_sum():
    summary = _summary_v36()
    for row in summary["runs"]:
        record = row["cavity_eigenmode_frequency_q_energy_orthogonality_degeneracy_mesh_owner_result_identity"]
        record["q_total"] = 100.0
        record["result_q_total"] = 100.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v36_public_rejects_self_consistent_active_array_power_imbalance():
    summary = _summary_v36()
    for row in summary["runs"]:
        record = row["active_sparameter_embedded_pattern_scan_impedance_power_frequency_mesh_owner_result_identity"]
        record["radiated_power_w"] = 1.0
        record["result_radiated_power_w"] = 1.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"
