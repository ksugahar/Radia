from copy import deepcopy
import math

from radia_mcp.radia_ngsolve.wave_energy_identity_v52 import ENERGY_BALANCE, EIGENMODE_Q, validate_public_v52_identity


PROMOTED_CASE_IDS = {
    "v52_public_energybalance_incident_reflected_transmitted_absorbed_dissipated_owner_mismatch",
    "v52_public_eigenmode_qfactor_storedenergy_boundaryloss_normalization_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _payload() -> dict[str, object]:
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    dissipated = {"ohmic": 40.0, "dielectric": 20.0}
    q_factor = 2.0 * math.pi * 1.0e9 * 0.01 / 3.0
    run = {
        ENERGY_BALANCE: {
            **_generations("energy-balance-v52", ("incident_generation", "scattered_generation", "absorbed_generation", "dissipated_generation", "owner_generation", "result_generation")),
            "incident_power_w": 100.0, "result_incident_power_w": 100.0,
            "reflected_power_w": 10.0, "result_reflected_power_w": 10.0,
            "transmitted_power_w": 30.0, "result_transmitted_power_w": 30.0,
            "absorbed_power_w": 60.0, "result_absorbed_power_w": 60.0,
            "dissipated_power_w": dissipated, "result_dissipated_power_w": dissipated,
            "run_owner": "run:energy-balance-v52", "result_run_owner": "run:energy-balance-v52", **result,
        },
        EIGENMODE_Q: {
            **_generations("eigenmode-q-v52", ("frequency_generation", "energy_generation", "loss_generation", "normalization_generation", "owner_generation", "result_generation")),
            "frequency_hz": 1.0e9, "result_frequency_hz": 1.0e9,
            "stored_energy_j": 0.01, "result_stored_energy_j": 0.01,
            "boundary_loss_w": 2.0, "result_boundary_loss_w": 2.0,
            "volume_loss_w": 1.0, "result_volume_loss_w": 1.0,
            "q_factor": q_factor, "result_q_factor": q_factor,
            "normalization": "physical_stored_energy", "result_normalization": "physical_stored_energy",
            "mode_owner": "mode:eigenmode-v52", "result_mode_owner": "mode:eigenmode-v52", **result,
        },
    }
    return {"runs": [deepcopy(run), deepcopy(run)]}


def test_v52_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_v52_identity(_payload()).values())


def test_v52_frozen_counterfactuals_are_rejected() -> None:
    payload = _payload()
    payload["runs"][0][ENERGY_BALANCE]["result_absorbed_power_w"] = 80.0
    payload["runs"][0][EIGENMODE_Q]["result_q_factor"] = 1.0
    assert not all(validate_public_v52_identity(payload).values())


def test_v52_self_consistent_wrong_physics_are_rejected() -> None:
    payload = _payload()
    for run in payload["runs"]:
        run[ENERGY_BALANCE]["absorbed_power_w"] = run[ENERGY_BALANCE]["result_absorbed_power_w"] = 80.0
        run[EIGENMODE_Q]["q_factor"] = run[EIGENMODE_Q]["result_q_factor"] = 1.0
    assert not all(validate_public_v52_identity(payload).values())
