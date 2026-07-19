from copy import deepcopy

from radia_mcp.radia_ngsolve.frequency_contact_identity_v50 import validate_public_v50_identity


PROMOTED_CASE_IDS = {
    "v50_public_frequency_sweep_complex_branch_phase_unit_dataset_interpolation_owner_mismatch",
    "v50_public_contact_pair_augmented_lagrange_penalty_gap_pressure_frame_owner_mismatch",
}


def _records() -> dict[str, object]:
    frequency_generation = "frequency-sweep-v50-1001"
    contact_generation = "contact-v50-1001"
    frequencies = [100.0, 1000.0, 10000.0]
    gaps = [0.0, 1e-6, 2e-6]
    pressure = [2e6, 1e6, 0.0]
    return {
        "frequency_sweep_complex_branch_phase_unit_dataset_interpolation_owner_identity": {
            "generation": frequency_generation,
            **{
                name: frequency_generation
                for name in (
                    "frequency_generation",
                    "branch_generation",
                    "phase_generation",
                    "unit_generation",
                    "dataset_generation",
                    "interpolation_generation",
                    "solution_generation",
                    "result_generation",
                )
            },
            "frequency_hz": frequencies,
            "result_frequency_hz": frequencies,
            "complex_branch": "positive_frequency",
            "result_complex_branch": "positive_frequency",
            "phase_convention": "exp(+jomega_t)",
            "result_phase_convention": "exp(+jomega_t)",
            "field_units": {"electric_field": "V/m", "magnetic_field": "A/m"},
            "result_field_units": {"electric_field": "V/m", "magnetic_field": "A/m"},
            "dataset_tag": "dataset:v50-frequency",
            "result_dataset_tag": "dataset:v50-frequency",
            "dataset_interpolation": "linear_complex",
            "result_dataset_interpolation": "linear_complex",
            "solution_owner": "solution:frequency-v50-1001",
            "result_solution_owner": "solution:frequency-v50-1001",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        "contact_pair_augmented_lagrange_penalty_gap_pressure_frame_owner_identity": {
            "generation": contact_generation,
            **{
                name: contact_generation
                for name in (
                    "pair_generation",
                    "method_generation",
                    "penalty_generation",
                    "gap_generation",
                    "pressure_generation",
                    "frame_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "contact_pair_id": "pair:source-destination-v50",
            "result_contact_pair_id": "pair:source-destination-v50",
            "contact_method": "augmented_lagrange",
            "result_contact_method": "augmented_lagrange",
            "penalty_factor": 1e9,
            "result_penalty_factor": 1e9,
            "gap_m": gaps,
            "result_gap_m": gaps,
            "pressure_pa": pressure,
            "result_pressure_pa": pressure,
            "coordinate_frame": "spatial",
            "result_coordinate_frame": "spatial",
            "contact_owner": "contact:v50-1001",
            "result_contact_owner": "contact:v50-1001",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
    }


def test_v50_public_positive_replay_is_accepted() -> None:
    result = validate_public_v50_identity(_records())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v50_public_mixed_frequency_sweep_is_rejected() -> None:
    records = deepcopy(_records())
    row = records["frequency_sweep_complex_branch_phase_unit_dataset_interpolation_owner_identity"]
    row.update(
        {
            "result_frequency_hz": [10000.0, 1000.0, 100.0],
            "result_complex_branch": "negative_frequency",
            "result_phase_convention": "exp(-jomega_t)",
            "result_field_units": {"electric_field": "mV/m", "magnetic_field": "A/m"},
            "result_dataset_tag": "dataset:old",
            "result_dataset_interpolation": "nearest_magnitude",
            "result_solution_owner": "solution:old",
        }
    )
    assert validate_public_v50_identity(records)["status"] == "needs_attention"


def test_v50_public_mixed_contact_state_is_rejected() -> None:
    records = deepcopy(_records())
    row = records["contact_pair_augmented_lagrange_penalty_gap_pressure_frame_owner_identity"]
    row.update(
        {
            "result_contact_pair_id": "pair:reversed",
            "result_contact_method": "penalty",
            "result_penalty_factor": 1e6,
            "result_gap_m": [2e-6, 1e-6, 0.0],
            "result_pressure_pa": [0.0, 1e6, 2e6],
            "result_coordinate_frame": "material",
            "result_contact_owner": "contact:old",
        }
    )
    assert validate_public_v50_identity(records)["status"] == "needs_attention"


def test_v50_public_invalid_canonical_sequences_are_rejected() -> None:
    records = deepcopy(_records())
    records["frequency_sweep_complex_branch_phase_unit_dataset_interpolation_owner_identity"]["frequency_hz"] = [1000.0, 100.0]
    records["contact_pair_augmented_lagrange_penalty_gap_pressure_frame_owner_identity"]["pressure_pa"] = [0.0, 1e6, 2e6]
    assert validate_public_v50_identity(records)["status"] == "needs_attention"
