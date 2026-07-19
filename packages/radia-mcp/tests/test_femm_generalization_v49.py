from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v49 import (
    HARMONIC,
    NONLINEAR,
    validate_public_identity,
)


PROMOTED_CASE_IDS = {
    "v49_public_nonlinear_bh_branch_incremental_permeability_temperature_lamination_owner_mismatch",
    "v49_public_harmonic_geometry_depth_frequency_phase_circuit_loss_owner_mismatch",
}


def _identity() -> dict[str, object]:
    nonlinear_generation = "nonlinear-material-v49-901"
    harmonic_generation = "harmonic-model-v49-901"
    bh_rows = [[0.0, 0.0], [0.7, 120.0], [1.35, 620.0], [1.62, 2400.0]]
    incremental = [[0.0023, 0.0001], [0.0001, 0.0020]]
    lamination = {"fill_factor": 0.95, "direction": "in-plane", "sheet_thickness_m": 0.00035}
    circuit = {"name": "coil-a", "current_a": [3.0, -1.0], "turns": 120}
    losses = {"joule_w": 4.25, "core_w": 1.75}
    return {
        NONLINEAR: {
            "generation": nonlinear_generation,
            "branch_generation": nonlinear_generation,
            "incremental_generation": nonlinear_generation,
            "temperature_generation": nonlinear_generation,
            "lamination_generation": nonlinear_generation,
            "result_generation": nonlinear_generation,
            "bh_branch": "ascending:first-quadrant",
            "result_bh_branch": "ascending:first-quadrant",
            "bh_rows_t_a_per_m": bh_rows,
            "result_bh_rows_t_a_per_m": bh_rows,
            "incremental_permeability_h_per_m": incremental,
            "result_incremental_permeability_h_per_m": incremental,
            "temperature_c": 80.0,
            "result_temperature_c": 80.0,
            "lamination": lamination,
            "result_lamination": lamination,
            "material_owner": "material:core-v49-901",
            "result_material_owner": "material:core-v49-901",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        HARMONIC: {
            "generation": harmonic_generation,
            "geometry_generation": harmonic_generation,
            "frequency_generation": harmonic_generation,
            "circuit_generation": harmonic_generation,
            "loss_generation": harmonic_generation,
            "result_generation": harmonic_generation,
            "geometry_type": "planar",
            "result_geometry_type": "planar",
            "depth_m": 0.04,
            "result_depth_m": 0.04,
            "frequency_hz": 400.0,
            "result_frequency_hz": 400.0,
            "phase_convention": "exp(+jwt)",
            "result_phase_convention": "exp(+jwt)",
            "circuit": circuit,
            "result_circuit": circuit,
            "losses_w": losses,
            "result_losses_w": losses,
            "result_owner": "harmonic-result:v49-901",
            "accepted_result_owner": "harmonic-result:v49-901",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
    }


def test_v49_positive_nonlinear_and_harmonic_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v49_nonlinear_branch_temperature_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[NONLINEAR]["result_bh_branch"] = "descending:minor-loop"
    identity[NONLINEAR]["result_temperature_c"] = 20.0
    identity[NONLINEAR]["result_material_owner"] = "material:old"
    assert validate_public_identity(identity)["v49_nonlinear_bh_incremental_temperature_lamination_owner"] is False


def test_v49_harmonic_geometry_phase_loss_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[HARMONIC]["result_geometry_type"] = "axisymmetric"
    identity[HARMONIC]["result_phase_convention"] = "exp(-jwt)"
    identity[HARMONIC]["result_losses_w"] = {"joule_w": 1.0, "core_w": 8.0}
    identity[HARMONIC]["accepted_result_owner"] = "harmonic-result:old"
    assert validate_public_identity(identity)["v49_harmonic_geometry_depth_frequency_phase_circuit_loss_owner"] is False
