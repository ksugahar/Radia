from copy import deepcopy

from radia_mcp.radia_ngsolve.femm_force_heat_identity_v56 import FORCE, HEAT, validate_public_identity


CASE_IDS = {
    "v56_public_magnetostatic_energy_coenergy_force_displacement_derivative_owner_mismatch",
    "v56_public_heatflow_joulesource_temperature_flux_balance_region_owner_mismatch",
}


def _identity() -> dict[str, object]:
    generation = "femm-public-v56-test"
    generations = lambda fields: {field: generation for field in fields}
    samples = [{"displacement_m": -1.0e-4, "coenergy_j": 0.9998}, {"displacement_m": 1.0e-4, "coenergy_j": 1.0002}]
    temperatures = {"minimum_k": 293.15, "maximum_k": 331.5, "mean_k": 307.2}
    regions = {"region:coil": 32.0, "region:core": 18.0}
    return {
        FORCE: {"generation": generation, **generations(("energy_generation", "coenergy_generation", "displacement_generation", "derivative_generation", "force_generation", "owner_generation", "result_generation")), "magnetic_energy_j": 0.98, "result_magnetic_energy_j": 0.98, "coenergy_samples": samples, "result_coenergy_samples": samples, "displacement_step_m": 1.0e-4, "result_displacement_step_m": 1.0e-4, "coenergy_derivative_n": 2.0, "result_coenergy_derivative_n": 2.0, "force_n": 2.0, "result_force_n": 2.0, "solution_owner": "solution:force-v56", "result_solution_owner": "solution:force-v56", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64},
        HEAT: {"generation": generation, **generations(("source_generation", "temperature_generation", "flux_generation", "balance_generation", "region_generation", "owner_generation", "result_generation")), "joule_source_w": 50.0, "result_joule_source_w": 50.0, "temperature_field_k": temperatures, "result_temperature_field_k": temperatures, "outward_boundary_flux_w": 50.0, "result_outward_boundary_flux_w": 50.0, "region_joule_balance_w": regions, "result_region_joule_balance_w": regions, "balance_residual_w": 0.0, "result_balance_residual_w": 0.0, "solution_owner": "solution:heat-v56", "result_solution_owner": "solution:heat-v56", "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64},
    }


def test_v56_positive_identity_is_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v56_frozen_result_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[FORCE]["result_force_n"] = -5.0
    identity[HEAT]["result_outward_boundary_flux_w"] = -10.0
    assert not all(validate_public_identity(identity).values())


def test_v56_self_consistent_physics_contradictions_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[FORCE]["force_n"] = identity[FORCE]["result_force_n"] = 3.0
    identity[HEAT]["outward_boundary_flux_w"] = identity[HEAT]["result_outward_boundary_flux_w"] = 40.0
    identity[HEAT]["balance_residual_w"] = identity[HEAT]["result_balance_residual_w"] = 10.0
    assert not all(validate_public_identity(identity).values())


def test_v56_malformed_samples_reject_without_raising() -> None:
    identity = deepcopy(_identity())
    identity[FORCE]["coenergy_samples"] = [[0.0]]
    assert not all(validate_public_identity(identity).values())
