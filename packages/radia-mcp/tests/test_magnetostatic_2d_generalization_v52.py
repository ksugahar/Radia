from copy import deepcopy

from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v52 import (
    FROZEN_INDUCTANCE,
    MAGNETIC_PRESSURE,
    validate_public_identity,
)


PROMOTED_CASE_IDS = {
    "v52_public_magneticpressure_fieldjump_boundarynormal_traction_owner_mismatch",
    "v52_public_frozen_permeability_incremental_inductance_biaspoint_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity() -> dict[str, object]:
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    field_jump = {"normal_h_a_per_m": 1200.0, "tangential_b_t": 0.42}
    normal = [0.6, 0.8, 0.0]
    traction = [3200.0, 4200.0, 0.0]
    bias = {"current_a": 8.0, "solution_sha256": "b" * 64}
    perturbation = {"delta_current_a": 0.02, "frequency_hz": 400.0}
    return {
        MAGNETIC_PRESSURE: {
            **_generations("magnetic-pressure-v52", ("field_generation", "normal_generation", "traction_generation", "owner_generation", "result_generation")),
            "field_jump": field_jump,
            "result_field_jump": field_jump,
            "boundary_normal": normal,
            "result_boundary_normal": normal,
            "traction_n_per_m2": traction,
            "result_traction_n_per_m2": traction,
            "integration_measure": "surface_area_m2",
            "result_integration_measure": "surface_area_m2",
            "field_owner": "field:magnetic-pressure-v52",
            "result_field_owner": "field:magnetic-pressure-v52",
            **result,
        },
        FROZEN_INDUCTANCE: {
            **_generations("frozen-inductance-v52", ("permeability_generation", "bias_generation", "perturbation_generation", "inductance_generation", "owner_generation", "result_generation")),
            "permeability_mode": "frozen_at_bias",
            "result_permeability_mode": "frozen_at_bias",
            "bias_point": bias,
            "result_bias_point": bias,
            "perturbation": perturbation,
            "result_perturbation": perturbation,
            "incremental_inductance_h": 0.014,
            "result_incremental_inductance_h": 0.014,
            "solution_owner": "solution:frozen-inductance-v52",
            "result_solution_owner": "solution:frozen-inductance-v52",
            **result,
        },
    }


def test_v52_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v52_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[MAGNETIC_PRESSURE]["result_boundary_normal"] = [-0.6, -0.8, 0.0]
    identity[FROZEN_INDUCTANCE]["result_bias_point"] = {"current_a": 0.0, "solution_sha256": "a" * 64}
    assert not all(validate_public_identity(identity).values())


def test_v52_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[MAGNETIC_PRESSURE]["boundary_normal"] = identity[MAGNETIC_PRESSURE]["result_boundary_normal"] = [2.0, 0.0, 0.0]
    identity[FROZEN_INDUCTANCE]["perturbation"] = identity[FROZEN_INDUCTANCE]["result_perturbation"] = {"delta_current_a": 2.0, "frequency_hz": 400.0}
    assert not all(validate_public_identity(identity).values())
