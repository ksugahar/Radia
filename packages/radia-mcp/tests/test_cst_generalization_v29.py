from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v28 import _summary_v28


_PROMOTED_CASE_IDS = (
    "v29_public_dispersive_vector_fit_pole_residue_passivity_causality_temperature_mismatch",
    "v29_public_array_embedded_element_pattern_feed_phase_active_reflection_scan_angle_mismatch",
)


def _summary_v29():
    summary = _summary_v28()
    for index, row in enumerate(summary["runs"]):
        generation = f"dispersive-fit-{341 + index}"
        row["dispersive_vector_fit_passivity_causality_temperature_generation_identity"] = {
            "fit_generation": generation, "pole_fit_generation": generation,
            "residue_fit_generation": generation, "passivity_fit_generation": generation,
            "causality_fit_generation": generation, "temperature_fit_generation": generation,
            "frequency_fit_generation": generation, "material_fit_generation": generation,
            "result_fit_generation": generation,
            "temperature_c": 25.0, "result_temperature_c": 25.0,
            "frequency_grid_hz": [1.0e9, 2.0e9, 4.0e9, 8.0e9],
            "result_frequency_grid_hz": [1.0e9, 2.0e9, 4.0e9, 8.0e9],
            "poles_rad_s": [[-1.0e9, 1.0e10], [-2.0e9, 2.0e10]],
            "result_poles_rad_s": [[-1.0e9, 1.0e10], [-2.0e9, 2.0e10]],
            "residues": [[1.0e9, 2.0e8], [5.0e8, 1.0e8]],
            "result_residues": [[1.0e9, 2.0e8], [5.0e8, 1.0e8]],
            "passivity_enforced": True, "result_passivity_enforced": True,
            "minimum_dissipation": 0.01, "result_minimum_dissipation": 0.01,
            "causality_residual": 1.0e-8, "result_causality_residual": 1.0e-8,
            "causality_residual_limit": 1.0e-6,
            "material_table_sha256": "1" * 64, "result_material_table_sha256": "1" * 64,
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        }
        generation = f"array-scan-{341 + index}"
        row["array_embedded_pattern_feed_phase_active_reflection_scan_generation_identity"] = {
            "array_generation": generation, "pattern_array_generation": generation,
            "element_array_generation": generation, "phase_array_generation": generation,
            "reflection_array_generation": generation, "scan_array_generation": generation,
            "power_array_generation": generation, "mesh_array_generation": generation,
            "result_array_generation": generation,
            "element_order": [1, 2, 3, 4], "result_element_order": [1, 2, 3, 4],
            "embedded_pattern_sha256": ["3" * 64, "4" * 64, "5" * 64, "6" * 64],
            "result_embedded_pattern_sha256": ["3" * 64, "4" * 64, "5" * 64, "6" * 64],
            "scan_angles_deg": [-30.0, 0.0, 30.0], "result_scan_angles_deg": [-30.0, 0.0, 30.0],
            "feed_phase_deg": [[0.0, -45.0, -90.0, -135.0], [0.0, 0.0, 0.0, 0.0], [0.0, 45.0, 90.0, 135.0]],
            "result_feed_phase_deg": [[0.0, -45.0, -90.0, -135.0], [0.0, 0.0, 0.0, 0.0], [0.0, 45.0, 90.0, 135.0]],
            "active_reflection_magnitude": [0.2, 0.1, 0.2],
            "result_active_reflection_magnitude": [0.2, 0.1, 0.2],
            "accepted_power_fraction": [0.96, 0.99, 0.96],
            "result_accepted_power_fraction": [0.96, 0.99, 0.96],
            "array_mesh_sha256": "7" * 64, "result_array_mesh_sha256": "7" * 64,
            "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
        }
    return summary


def test_v29_public_positive_dispersive_fit_and_array_scan_identities():
    assert nonlinear_inductance_sweep_gate(_summary_v29())["status"] == "ok"


def test_v29_public_rejects_dispersive_fit_identity_mismatch():
    summary = _summary_v29()
    identity = summary["runs"][0]["dispersive_vector_fit_passivity_causality_temperature_generation_identity"]
    identity.update({
        "pole_fit_generation": "dispersive-fit-340", "temperature_fit_generation": "dispersive-fit-339",
        "result_temperature_c": 100.0, "result_frequency_grid_hz": [1.0e9, 3.0e9, 9.0e9],
        "result_poles_rad_s": [[1.0e9, 1.0e10]], "result_residues": [[-1.0e9, 0.0]],
        "result_passivity_enforced": False, "result_minimum_dissipation": -0.1,
        "result_causality_residual": 1.0e-2, "result_material_table_sha256": "e" * 64,
        "accepted_result_sha256": "f" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "dispersive_vector_fit_uses_current_stable_poles_residues_passivity_causality_temperature_and_result"
    ]


def test_v29_public_rejects_array_scan_identity_mismatch():
    summary = _summary_v29()
    identity = summary["runs"][0]["array_embedded_pattern_feed_phase_active_reflection_scan_generation_identity"]
    identity.update({
        "pattern_array_generation": "array-scan-340", "phase_array_generation": "array-scan-339",
        "result_element_order": [4, 3, 2, 1], "result_embedded_pattern_sha256": ["e" * 64],
        "result_scan_angles_deg": [30.0, 0.0, -30.0], "result_feed_phase_deg": [[0.0, 90.0]],
        "result_active_reflection_magnitude": [1.2, 0.1], "result_accepted_power_fraction": [0.2, 1.1],
        "result_array_mesh_sha256": "f" * 64, "accepted_result_sha256": "a" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "array_scan_uses_current_embedded_patterns_element_order_phases_reflection_power_mesh_and_result"
    ]
