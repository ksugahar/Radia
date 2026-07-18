from __future__ import annotations

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v34 import _identity_v34


_PROMOTED_CASE_IDS = (
    "v35_public_inductance_matrix_reciprocity_psd_fluxlinkage_current_energy_mismatch",
    "v35_public_nonlinear_bh_interpolation_differential_permeability_energy_coenergy_branch_mismatch",
)


def _identity_v35():
    identity = _identity_v34()
    generation = "inductance-matrix-235"
    identity[
        "inductance_matrix_reciprocity_psd_fluxlinkage_current_energy_coil_mesh_owner_result_identity"
    ] = {
        "inductance_generation": generation,
        **{key: generation for key in (
            "reciprocity_generation", "psd_generation", "flux_generation",
            "current_generation", "energy_generation", "coil_generation",
            "mesh_generation", "owner_generation", "result_generation")},
        "coil_order": ["phase_a", "phase_b"], "result_coil_order": ["phase_a", "phase_b"],
        "currents_a": [2.0, -1.0], "result_currents_a": [2.0, -1.0],
        "inductance_matrix_h": [[0.008, 0.002], [0.002, 0.006]],
        "result_inductance_matrix_h": [[0.008, 0.002], [0.002, 0.006]],
        "flux_linkages_wb_turn": [0.014, -0.002],
        "result_flux_linkages_wb_turn": [0.014, -0.002],
        "stored_energy_j": 0.015, "result_stored_energy_j": 0.015,
        "minimum_eigenvalue_h": 0.00476393202250021,
        "result_minimum_eigenvalue_h": 0.00476393202250021,
        "coil_mesh_sha256": "1" * 64, "result_coil_mesh_sha256": "1" * 64,
        "inductance_result_owner": "magnetic/coils-235",
        "accepted_inductance_result_owner": "magnetic/coils-235",
        "inductance_result_sha256": "2" * 64,
        "accepted_inductance_result_sha256": "2" * 64,
    }
    generation = "nonlinear-bh-235"
    identity[
        "nonlinear_bh_interpolation_differential_permeability_branch_energy_coenergy_operating_material_solution_identity"
    ] = {
        "bh_generation": generation,
        **{key: generation for key in (
            "interpolation_generation", "differential_generation", "branch_generation",
            "energy_generation", "coenergy_generation", "operating_generation",
            "material_generation", "solution_generation", "result_generation")},
        "b_samples_t": [0.0, 0.5, 1.0, 1.4], "result_b_samples_t": [0.0, 0.5, 1.0, 1.4],
        "h_samples_a_m": [0.0, 100.0, 400.0, 1200.0],
        "result_h_samples_a_m": [0.0, 100.0, 400.0, 1200.0],
        "differential_permeability_h_m": [0.005, 1.0 / 600.0, 0.0005],
        "result_differential_permeability_h_m": [0.005, 1.0 / 600.0, 0.0005],
        "branch": "ascending", "result_branch": "ascending",
        "operating_point_b_t": 1.4, "result_operating_point_b_t": 1.4,
        "operating_point_h_a_m": 1200.0, "result_operating_point_h_a_m": 1200.0,
        "magnetic_energy_density_j_m3": 470.0, "result_magnetic_energy_density_j_m3": 470.0,
        "magnetic_coenergy_density_j_m3": 1210.0,
        "result_magnetic_coenergy_density_j_m3": 1210.0,
        "material_owner": "materials/nonlinear-core-235",
        "accepted_material_owner": "materials/nonlinear-core-235",
        "nonlinear_solution_sha256": "3" * 64,
        "accepted_nonlinear_solution_sha256": "3" * 64,
    }
    return identity


def test_v35_public_positive_inductance_and_nonlinear_bh_closure():
    assert _gate(_identity_v35())["status"] == "ok"


def test_v35_public_inductance_matrix_reciprocity_psd_fluxlinkage_current_energy_mismatch():
    identity = _identity_v35()
    identity["inductance_matrix_reciprocity_psd_fluxlinkage_current_energy_coil_mesh_owner_result_identity"].update({
        "reciprocity_generation": "inductance-matrix-234", "energy_generation": "inductance-matrix-233",
        "result_generation": "inductance-matrix-232", "result_coil_order": ["phase_b", "phase_a"],
        "result_currents_a": [-1.0, 2.0], "result_inductance_matrix_h": [[0.008, 0.004], [-0.003, -0.006]],
        "result_flux_linkages_wb_turn": [0.2, 0.1], "result_stored_energy_j": -0.015,
        "result_minimum_eigenvalue_h": -0.01, "result_coil_mesh_sha256": "9" * 64,
        "accepted_inductance_result_owner": "stale/coils", "accepted_inductance_result_sha256": "a" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["inductance_matrix_closes_reciprocity_psd_flux_current_energy_coils_mesh_owner_and_result"]


def test_v35_public_nonlinear_bh_interpolation_differential_permeability_energy_coenergy_branch_mismatch():
    identity = _identity_v35()
    identity["nonlinear_bh_interpolation_differential_permeability_branch_energy_coenergy_operating_material_solution_identity"].update({
        "interpolation_generation": "nonlinear-bh-234", "energy_generation": "nonlinear-bh-233",
        "result_generation": "nonlinear-bh-232", "result_b_samples_t": [0.0, 0.5, 0.4, 1.4],
        "result_h_samples_a_m": [0.0, 100.0, 80.0, 1200.0],
        "result_differential_permeability_h_m": [0.005, -0.005, 0.001],
        "result_branch": "descending", "result_operating_point_b_t": 1.0,
        "result_operating_point_h_a_m": 400.0, "result_magnetic_energy_density_j_m3": -470.0,
        "result_magnetic_coenergy_density_j_m3": 470.0,
        "accepted_material_owner": "stale/material", "accepted_nonlinear_solution_sha256": "b" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["nonlinear_bh_closes_interpolation_differential_mu_branch_energy_coenergy_operating_material_and_solution"]


def test_v35_rejects_self_consistent_inductance_energy_mismatch():
    identity = _identity_v35()
    row = identity["inductance_matrix_reciprocity_psd_fluxlinkage_current_energy_coil_mesh_owner_result_identity"]
    row["stored_energy_j"] = row["result_stored_energy_j"] = 0.03
    assert _gate(identity)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_nonmonotone_bh_branch():
    identity = _identity_v35()
    row = identity["nonlinear_bh_interpolation_differential_permeability_branch_energy_coenergy_operating_material_solution_identity"]
    row["h_samples_a_m"] = row["result_h_samples_a_m"] = [0.0, 100.0, 80.0, 1200.0]
    assert _gate(identity)["status"] == "needs_attention"
