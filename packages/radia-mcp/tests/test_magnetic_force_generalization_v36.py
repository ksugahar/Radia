from __future__ import annotations

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v35 import _identity_v35


_PROMOTED_CASE_IDS = (
    "v36_public_nonlinear_bh_incremental_permeability_energy_coenergy_differential_inductance_mismatch",
    "v36_public_weighted_stress_force_region_weighting_contour_independence_mesh_convergence_mismatch",
)


def _identity_v36():
    identity = _identity_v35()
    generation = "nonlinear-incremental-236"
    identity[
        "nonlinear_bh_incremental_permeability_energy_coenergy_differential_inductance_current_mesh_owner_solution_identity"
    ] = {
        "nonlinear_generation": generation,
        **{key: generation for key in (
            "branch_generation", "incremental_generation", "energy_generation",
            "coenergy_generation", "inductance_generation", "current_generation",
            "mesh_generation", "owner_generation", "solution_generation", "result_generation")},
        "bh_branch": "ascending", "result_bh_branch": "ascending",
        "current_points_a": [1.0, 2.0, 3.0], "result_current_points_a": [1.0, 2.0, 3.0],
        "flux_linkages_wb_turn": [0.01, 0.018, 0.024],
        "result_flux_linkages_wb_turn": [0.01, 0.018, 0.024],
        "incremental_permeability_h_m": 0.0012, "result_incremental_permeability_h_m": 0.0012,
        "differential_inductance_h": 0.006, "result_differential_inductance_h": 0.006,
        "magnetic_energy_j": 0.03, "result_magnetic_energy_j": 0.03,
        "magnetic_coenergy_j": 0.042, "result_magnetic_coenergy_j": 0.042,
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "field_owner": "nonlinear/core-236", "accepted_field_owner": "nonlinear/core-236",
        "solution_sha256": "2" * 64, "accepted_solution_sha256": "2" * 64,
    }
    generation = "weighted-stress-force-236"
    identity[
        "weighted_stress_force_region_air_contour_mesh_convergence_direction_owner_result_identity"
    ] = {
        "force_generation": generation,
        **{key: generation for key in (
            "weighting_generation", "air_generation", "contour_generation",
            "convergence_generation", "direction_generation", "field_generation",
            "owner_generation", "result_generation")},
        "weighting_region_id": "air:force-mask", "result_weighting_region_id": "air:force-mask",
        "air_enclosure_id": "air:enclosure", "result_air_enclosure_id": "air:enclosure",
        "weighted_force_n": [12.0, -0.2], "result_weighted_force_n": [12.0, -0.2],
        "contour_force_samples_n": [[12.0, -0.2], [12.01, -0.19], [11.99, -0.21]],
        "result_contour_force_samples_n": [[12.0, -0.2], [12.01, -0.19], [11.99, -0.21]],
        "mesh_sizes_m": [0.004, 0.002, 0.001], "result_mesh_sizes_m": [0.004, 0.002, 0.001],
        "mesh_force_sequence_n": [[11.5, -0.3], [11.9, -0.22], [12.0, -0.2]],
        "result_mesh_force_sequence_n": [[11.5, -0.3], [11.9, -0.22], [12.0, -0.2]],
        "force_direction_unit": [1.0, 0.0], "result_force_direction_unit": [1.0, 0.0],
        "field_owner": "force/field-236", "accepted_field_owner": "force/field-236",
        "force_result_sha256": "3" * 64, "accepted_force_result_sha256": "3" * 64,
    }
    return identity


def test_v36_public_positive_nonlinear_incremental_and_weighted_force_closure():
    assert _gate(_identity_v36())["status"] == "ok"


def test_v36_public_nonlinear_bh_incremental_permeability_energy_coenergy_differential_inductance_mismatch():
    identity = _identity_v36()
    identity["nonlinear_bh_incremental_permeability_energy_coenergy_differential_inductance_current_mesh_owner_solution_identity"].update({
        "branch_generation": "nonlinear-incremental-235", "mesh_generation": "nonlinear-incremental-234",
        "result_generation": "nonlinear-incremental-233", "result_bh_branch": "descending",
        "result_current_points_a": [3.0, 2.0, 1.0], "result_flux_linkages_wb_turn": [0.024, 0.018, 0.01],
        "result_incremental_permeability_h_m": -0.0012, "result_differential_inductance_h": -0.006,
        "result_magnetic_energy_j": -0.03, "result_magnetic_coenergy_j": 0.01,
        "result_mesh_sha256": "a" * 64, "accepted_field_owner": "stale/core",
        "accepted_solution_sha256": "b" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["nonlinear_incremental_solution_closes_branch_mu_energy_coenergy_differential_inductance_current_mesh_owner_and_result"]


def test_v36_public_weighted_stress_force_region_weighting_contour_independence_mesh_convergence_mismatch():
    identity = _identity_v36()
    identity["weighted_stress_force_region_air_contour_mesh_convergence_direction_owner_result_identity"].update({
        "weighting_generation": "weighted-stress-force-235", "direction_generation": "weighted-stress-force-234",
        "result_generation": "weighted-stress-force-233", "result_weighting_region_id": "iron:body",
        "result_air_enclosure_id": "air:old", "result_weighted_force_n": [-12.0, 2.0],
        "result_contour_force_samples_n": [[3.0, 4.0], [-5.0, 1.0]],
        "result_mesh_sizes_m": [0.001, 0.004], "result_mesh_force_sequence_n": [[3.0, 0.0], [20.0, 0.0]],
        "result_force_direction_unit": [-1.0, 2.0], "accepted_field_owner": "stale/field",
        "accepted_force_result_sha256": "c" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["weighted_stress_force_closes_region_air_contours_mesh_convergence_direction_owner_and_result"]


def test_v36_rejects_self_consistent_incremental_energy_partition_error():
    identity = _identity_v36()
    row = identity["nonlinear_bh_incremental_permeability_energy_coenergy_differential_inductance_current_mesh_owner_solution_identity"]
    row["magnetic_coenergy_j"] = row["result_magnetic_coenergy_j"] = 0.02
    assert _gate(identity)["status"] == "needs_attention"


def test_v36_rejects_self_consistent_nonconvergent_weighted_force_sequence():
    identity = _identity_v36()
    row = identity["weighted_stress_force_region_air_contour_mesh_convergence_direction_owner_result_identity"]
    row["mesh_force_sequence_n"] = row["result_mesh_force_sequence_n"] = [[11.5, -0.3], [10.0, -0.2], [12.0, -0.2]]
    assert _gate(identity)["status"] == "needs_attention"
