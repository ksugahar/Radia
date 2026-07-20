from copy import deepcopy

from radia_mcp.radia_ngsolve.coupled_periodic_identity_v53 import ELECTROTHERMAL, FLOQUET, validate_public_v53_identity


CASE_IDS = {
    "v53_public_electrothermal_contact_resistance_power_heat_reciprocity_time_owner_mismatch",
    "v53_public_floquet_periodic_phase_wavevector_boundary_pair_orientation_owner_mismatch",
}


def _records() -> dict[str, object]:
    generation = "comsol-public-v53"
    generations = lambda names: {name: generation for name in names}
    pair = {"source": "left", "destination": "right"}
    return {
        ELECTROTHERMAL: {
            "generation": generation, **generations(("contact_generation", "electric_generation", "thermal_generation", "time_generation", "owner_generation", "result_generation")),
            "contact_resistance_ohm": 0.02, "result_contact_resistance_ohm": 0.02,
            "contact_current_a": 3.0, "result_contact_current_a": 3.0,
            "electric_power_w": 0.18, "result_electric_power_w": 0.18,
            "deposited_heat_w": 0.18, "result_deposited_heat_w": 0.18,
            "time_s": 0.5, "result_time_s": 0.5,
            "solution_owner": "solution:electrothermal-v53", "result_solution_owner": "solution:electrothermal-v53",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        FLOQUET: {
            "generation": generation, **generations(("phase_generation", "wavevector_generation", "pair_generation", "orientation_generation", "owner_generation", "result_generation")),
            "phase_rad": 0.6, "result_phase_rad": 0.6,
            "wave_vector_per_m": [20.0, 0.0, 0.0], "result_wave_vector_per_m": [20.0, 0.0, 0.0],
            "translation_m": [0.03, 0.0, 0.0], "result_translation_m": [0.03, 0.0, 0.0],
            "boundary_pair": pair, "result_boundary_pair": pair,
            "orientation": "source_to_destination", "result_orientation": "source_to_destination",
            "field_owner": "field:floquet-v53", "result_field_owner": "field:floquet-v53",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def test_v53_public_positive_replay_is_accepted() -> None:
    assert validate_public_v53_identity(_records())["status"] == "ok"


def test_v53_frozen_public_mutations_are_rejected() -> None:
    value = deepcopy(_records())
    value[ELECTROTHERMAL].update({"result_electric_power_w": 0.36, "result_deposited_heat_w": 0.12, "result_solution_owner": "solution:stale"})
    value[FLOQUET].update({"result_phase_rad": -0.6, "result_boundary_pair": {"source": "right", "destination": "left"}, "result_field_owner": "field:stale"})
    assert validate_public_v53_identity(value)["status"] == "needs_attention"


def test_v53_self_consistent_energy_and_phase_errors_are_rejected() -> None:
    value = deepcopy(_records())
    value[ELECTROTHERMAL]["electric_power_w"] = value[ELECTROTHERMAL]["result_electric_power_w"] = 0.12
    value[ELECTROTHERMAL]["deposited_heat_w"] = value[ELECTROTHERMAL]["result_deposited_heat_w"] = 0.12
    value[FLOQUET]["phase_rad"] = value[FLOQUET]["result_phase_rad"] = -0.6
    assert validate_public_v53_identity(value)["status"] == "needs_attention"
