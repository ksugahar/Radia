from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.network_semantic_identity_v48 import FARFIELD, THERMAL, validate_public_v48_identity


PROMOTED_CASE_IDS = {
    "v48_public_farfield_polarization_basis_normalization_radiated_power_angular_grid_owner_mismatch",
    "v48_public_em_thermal_loss_mapping_mesh_interpolation_time_average_frequency_owner_mismatch",
}


def _payload() -> dict[str, object]:
    farfield_generation = "farfield-v48"
    thermal_generation = "em-thermal-v48"
    theta = [0.0, 45.0, 90.0]
    phi = [0.0, 90.0, 180.0]
    samples = [[1.0, 0.2], [0.8, 0.3], [0.4, 0.1]]
    components = ["conductor", "dielectric"]
    losses = [2.5, 0.5]
    return {"runs": [{
        FARFIELD: {
            "generation": farfield_generation,
            **{key: farfield_generation for key in ("polarization_generation", "normalization_generation", "power_generation", "angular_grid_generation", "monitor_generation", "result_generation")},
            "polarization_basis": "spherical_theta_phi", "result_polarization_basis": "spherical_theta_phi",
            "normalization": "accepted_radiated_power", "result_normalization": "accepted_radiated_power",
            "radiated_power_w": 4.0, "result_radiated_power_w": 4.0,
            "theta_deg": theta, "result_theta_deg": theta, "phi_deg": phi, "result_phi_deg": phi,
            "field_samples": samples, "result_field_samples": samples,
            "monitor_owner": "monitor:farfield-v48", "result_monitor_owner": "monitor:farfield-v48",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        THERMAL: {
            "generation": thermal_generation,
            **{key: thermal_generation for key in ("loss_mapping_generation", "source_mesh_generation", "target_mesh_generation", "interpolation_generation", "time_average_generation", "frequency_generation", "task_generation", "result_generation")},
            "loss_component_ids": components, "mapped_loss_component_ids": components,
            "loss_w": losses, "mapped_loss_w": losses,
            "source_mesh_sha256": "2" * 64, "mapped_source_mesh_sha256": "2" * 64,
            "target_mesh_sha256": "3" * 64, "mapped_target_mesh_sha256": "3" * 64,
            "interpolation_method": "conservative_nodal", "mapped_interpolation_method": "conservative_nodal",
            "time_average": "cycle_average", "mapped_time_average": "cycle_average",
            "frequency_hz": 2.45e9, "mapped_frequency_hz": 2.45e9,
            "task_owner": "task:em-thermal-v48", "mapped_task_owner": "task:em-thermal-v48",
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        },
    }]}


def test_v48_positive_network_artifacts_are_accepted() -> None:
    assert all(validate_public_v48_identity(_payload()).values())


def test_v48_farfield_identity_mutation_is_rejected() -> None:
    payload = _payload()
    payload["runs"][0][FARFIELD].update({"result_polarization_basis": "cartesian_xy", "result_normalization": "peak_field", "result_theta_deg": [90.0, 45.0, 0.0], "result_monitor_owner": "monitor:old"})
    assert not all(validate_public_v48_identity(payload).values())


def test_v48_em_thermal_mapping_mutation_is_rejected() -> None:
    payload = _payload()
    payload["runs"][0][THERMAL].update({"mapped_loss_component_ids": ["dielectric", "conductor"], "mapped_interpolation_method": "nearest", "mapped_time_average": "peak", "mapped_frequency_hz": 2.4e9})
    assert not all(validate_public_v48_identity(payload).values())


def test_v48_self_consistent_nonphysical_conventions_are_rejected() -> None:
    payload = deepcopy(_payload())
    farfield = payload["runs"][0][FARFIELD]
    farfield["polarization_basis"] = farfield["result_polarization_basis"] = "cartesian_xy"
    thermal = payload["runs"][0][THERMAL]
    thermal["interpolation_method"] = thermal["mapped_interpolation_method"] = "nearest"
    assert not all(validate_public_v48_identity(payload).values())
