from copy import deepcopy

from radia_mcp.radia_ngsolve.multiphysics_identity_v56 import INDUCTION, MODAL, validate_public_v56_identity


CASE_IDS = {
    "v56_public_acoustic_modalparticipation_frequencyresponse_normalization_energy_owner_mismatch",
    "v56_public_inductionheating_coilinput_jouleheat_thermalrise_time_owner_mismatch",
}


def _records() -> dict[str, object]:
    modal_generation = "modal-public-v56"
    induction_generation = "induction-public-v56"
    modes = [{"mode": "mode:1", "frequency_hz": 500.0, "participation": 0.8, "normalized_energy": 0.8}, {"mode": "mode:2", "frequency_hz": 1000.0, "participation": 0.2, "normalized_energy": 0.2}]
    response = [{"frequency_hz": 500.0, "pressure_pa": 2.0}, {"frequency_hz": 1000.0, "pressure_pa": 0.5}]
    times = [0.0, 0.5, 1.0]
    temperatures = [293.15, 303.15, 313.15]
    return {
        MODAL: {
            "generation": modal_generation, **{field: modal_generation for field in ("participation_generation", "response_generation", "normalization_generation", "energy_generation", "owner_generation", "result_generation")},
            "modal_rows": modes, "result_modal_rows": modes,
            "frequency_response": response, "result_frequency_response": response,
            "normalization": "unit_total_modal_energy", "result_normalization": "unit_total_modal_energy",
            "total_normalized_energy": 1.0, "result_total_normalized_energy": 1.0,
            "solution_owner": "solution:modal-v56", "result_solution_owner": "solution:modal-v56",
            "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        INDUCTION: {
            "generation": induction_generation, **{field: induction_generation for field in ("input_generation", "joule_generation", "thermal_generation", "time_generation", "owner_generation", "result_generation")},
            "coil_input_energy_j": 100.0, "result_coil_input_energy_j": 100.0,
            "joule_heat_energy_j": 80.0, "result_joule_heat_energy_j": 80.0,
            "stored_thermal_energy_j": 15.0, "result_stored_thermal_energy_j": 15.0,
            "boundary_heat_loss_j": 5.0, "result_boundary_heat_loss_j": 5.0,
            "time_s": times, "result_time_s": times,
            "average_temperature_k": temperatures, "result_average_temperature_k": temperatures,
            "solution_owner": "solution:induction-v56", "result_solution_owner": "solution:induction-v56",
            "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }


def test_v56_public_positive_replay_is_accepted() -> None:
    assert validate_public_v56_identity(_records())["status"] == "ok"


def test_v56_frozen_public_mutations_are_rejected() -> None:
    value = deepcopy(_records())
    value[MODAL].update({"result_modal_rows": [], "result_normalization": "peak_pressure", "result_solution_owner": "solution:stale"})
    value[INDUCTION].update({"result_coil_input_energy_j": 50.0, "result_time_s": [0.0, 2.0], "result_solution_owner": "solution:stale"})
    assert validate_public_v56_identity(value)["status"] == "needs_attention"


def test_v56_self_consistent_energy_and_history_contradictions_are_rejected() -> None:
    value = deepcopy(_records())
    value[MODAL]["total_normalized_energy"] = value[MODAL]["result_total_normalized_energy"] = 2.0
    value[INDUCTION]["average_temperature_k"] = value[INDUCTION]["result_average_temperature_k"] = [313.15, 303.15, 293.15]
    assert validate_public_v56_identity(value)["status"] == "needs_attention"


def test_v56_malformed_values_reject_without_raising() -> None:
    value = deepcopy(_records())
    value[MODAL]["modal_rows"] = [{"mode": ["mode:1"]}]
    value[INDUCTION]["time_s"] = [[0.0]]
    assert validate_public_v56_identity(value)["status"] == "needs_attention"
