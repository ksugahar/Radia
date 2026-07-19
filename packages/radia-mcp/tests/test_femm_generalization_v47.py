from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_artifact_lineage_v47 import BH, FORCE, validate_public_identity


PROMOTED_CASE_IDS = {
    "v47_public_force_method_body_owner_sign_displacement_pair_causal_mismatch",
    "v47_public_nonlinear_bh_operating_point_row_hysteresis_branch_mapping_mismatch",
}


def _identity() -> dict[str, object]:
    force_generation = "force-v47"
    bh_generation = "bh-v47"
    keys = ["current=1", "current=2", "current=3"]
    branches = ["ascending", "ascending", "descending"]
    return {
        FORCE: {
            "generation": force_generation,
            "force_method_generation": force_generation,
            "body_owner_generation": force_generation,
            "displacement_pair_generation": force_generation,
            "result_generation": force_generation,
            "force_method": "weighted_stress_tensor",
            "result_force_method": "weighted_stress_tensor",
            "body_owner": "group:moving",
            "result_body_owner": "group:moving",
            "force_sign_convention": "positive_displacement_direction",
            "result_force_sign_convention": "positive_displacement_direction",
            "displacement_pair_m": [0.0, 0.001],
            "result_displacement_pair_m": [0.0, 0.001],
            "coenergy_pair_j": [1.0, 1.01],
            "result_coenergy_pair_j": [1.0, 1.01],
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        BH: {
            "generation": bh_generation,
            "operating_point_generation": bh_generation,
            "branch_generation": bh_generation,
            "history_generation": bh_generation,
            "result_generation": bh_generation,
            "operating_point_row_keys": keys,
            "result_operating_point_row_keys": keys,
            "hysteresis_branches": branches,
            "result_hysteresis_branches": branches,
            "excitation_history_sha256": "2" * 64,
            "result_excitation_history_sha256": "2" * 64,
            "material_owner": "material:steel",
            "result_material_owner": "material:steel",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def test_v47_positive_replays_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v47_force_owner_sign_pair_mutation_is_rejected() -> None:
    identity = _identity()
    identity[FORCE]["result_force_method"] = "contour_stress"
    identity[FORCE]["result_body_owner"] = "group:fixed"
    identity[FORCE]["result_force_sign_convention"] = "negative_displacement_direction"
    identity[FORCE]["result_displacement_pair_m"] = [0.001, 0.0]
    assert not all(validate_public_identity(identity).values())


def test_v47_bh_row_branch_history_mutation_is_rejected() -> None:
    identity = _identity()
    identity[BH]["result_operating_point_row_keys"] = ["current=3", "current=1", "current=2"]
    identity[BH]["result_hysteresis_branches"] = ["descending", "ascending", "ascending"]
    identity[BH]["result_excitation_history_sha256"] = "a" * 64
    assert not all(validate_public_identity(identity).values())
