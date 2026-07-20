from copy import deepcopy

from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v51 import (
    INCREMENTAL,
    WEIGHTED_FORCE,
    validate_public_identity,
)


PROMOTED_CASE_IDS = {
    "v51_public_incremental_frozen_permeability_bias_harmonic_tangent_branch_owner_mismatch",
    "v51_public_weighted_stress_tensor_mask_air_elements_axisym_factor_force_frame_owner_mismatch",
}


def _identity() -> dict[str, object]:
    generation = "femm-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    tangent = {"radial": 410.0, "tangential": 370.0}
    force = [12.5, -0.4]
    return {
        INCREMENTAL: {
            "generation": generation, "bias_generation": generation, "harmonic_generation": generation,
            "tangent_generation": generation, "branch_generation": generation, "owner_generation": generation,
            "result_generation": generation, "analysis_mode": "incremental_permeability",
            "result_analysis_mode": "incremental_permeability", "frozen_bias_solution_sha256": "1" * 64,
            "result_frozen_bias_solution_sha256": "1" * 64, "harmonic_frequency_hz": 1000.0,
            "result_harmonic_frequency_hz": 1000.0, "tangent_permeability_relative": tangent,
            "result_tangent_permeability_relative": tangent, "branch_state": "ascending_major_loop",
            "result_branch_state": "ascending_major_loop", "operating_point_owner": "operating-point:v51",
            "result_operating_point_owner": "operating-point:v51", **result,
        },
        WEIGHTED_FORCE: {
            "generation": generation, "mask_generation": generation, "air_generation": generation,
            "axisym_generation": generation, "force_generation": generation, "frame_generation": generation,
            "owner_generation": generation, "result_generation": generation, "weighted_stress_mask_sha256": "2" * 64,
            "result_weighted_stress_mask_sha256": "2" * 64, "air_element_ids": [101, 102, 103],
            "result_air_element_ids": [101, 102, 103], "axisymmetric_radius_m": 0.025,
            "result_axisymmetric_radius_m": 0.025, "axisymmetric_factor_m": 2.0 * 3.141592653589793 * 0.025,
            "result_axisymmetric_factor_m": 2.0 * 3.141592653589793 * 0.025, "force_n": force,
            "result_force_n": force, "force_frame": "global_rz", "result_force_frame": "global_rz",
            "force_owner": "force:weighted-v51", "result_force_owner": "force:weighted-v51", **result,
        },
    }


def test_v51_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v51_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[INCREMENTAL].update({"result_analysis_mode": "frozen_permeability", "result_branch_state": "recoil_branch"})
    identity[WEIGHTED_FORCE].update({"result_air_element_ids": [101, 105], "result_force_frame": "local_xy"})
    assert not all(validate_public_identity(identity).values())


def test_v51_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[INCREMENTAL]["analysis_mode"] = identity[INCREMENTAL]["result_analysis_mode"] = "frozen_permeability"
    identity[WEIGHTED_FORCE]["axisymmetric_factor_m"] = identity[WEIGHTED_FORCE]["result_axisymmetric_factor_m"] = 1.0
    assert not all(validate_public_identity(identity).values())
