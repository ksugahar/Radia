from copy import deepcopy

from radia_mcp.radia_ngsolve.modal_continuation_identity_v51 import validate_public_v51_identity


PROMOTED_CASE_IDS = {
    "v51_public_eigenmode_frequency_normalization_phase_degenerate_subspace_mesh_owner_mismatch",
    "v51_public_continuation_branch_predictor_corrector_loadpath_turningpoint_solution_owner_mismatch",
}


def _records() -> dict[str, object]:
    eigen = "eigenmode-v51-1001"
    continuation = "continuation-v51-1001"
    frequencies = [1240.0, 1240.0]
    phase = {"dof": 17, "component": "real", "sign": "positive"}
    states = ["predictor:42", "corrector:42"]
    load_path = [0.0, 0.35, 0.7, 0.93, 0.88]
    return {
        "eigenmode_frequency_normalization_phase_subspace_mesh_owner_identity": {
            "generation": eigen,
            **{name: eigen for name in ("frequency_generation", "normalization_generation", "phase_generation", "subspace_generation", "mesh_generation", "owner_generation", "result_generation")},
            "frequency_hz": frequencies,
            "result_frequency_hz": frequencies,
            "normalization": "unit_generalized_mass",
            "result_normalization": "unit_generalized_mass",
            "phase_anchor": phase,
            "result_phase_anchor": phase,
            "degenerate_subspace_basis_sha256": "1" * 64,
            "result_degenerate_subspace_basis_sha256": "1" * 64,
            "mesh_revision": "mesh:v51-r3",
            "result_mesh_revision": "mesh:v51-r3",
            "mode_owner": "mode-set:v51-1001",
            "result_mode_owner": "mode-set:v51-1001",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
        "continuation_branch_predictor_corrector_loadpath_turningpoint_owner_identity": {
            "generation": continuation,
            **{name: continuation for name in ("branch_generation", "state_generation", "loadpath_generation", "turningpoint_generation", "owner_generation", "result_generation")},
            "branch_id": "branch:v51-primary",
            "result_branch_id": "branch:v51-primary",
            "predictor_corrector_states": states,
            "result_predictor_corrector_states": states,
            "load_path": load_path,
            "result_load_path": load_path,
            "turning_point_index": 3,
            "result_turning_point_index": 3,
            "solution_owner": "solution:continuation-v51-1001",
            "result_solution_owner": "solution:continuation-v51-1001",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def test_v51_public_positive_replay_is_accepted() -> None:
    result = validate_public_v51_identity(_records())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v51_public_mixed_eigenmode_identity_is_rejected() -> None:
    records = deepcopy(_records())
    row = records["eigenmode_frequency_normalization_phase_subspace_mesh_owner_identity"]
    row.update({"result_frequency_hz": [1235.0, 1245.0], "result_normalization": "unit_peak", "result_phase_anchor": {"dof": 21, "component": "imag", "sign": "negative"}, "result_degenerate_subspace_basis_sha256": "a" * 64, "result_mesh_revision": "mesh:v51-r2", "result_mode_owner": "mode-set:stale"})
    assert validate_public_v51_identity(records)["status"] == "needs_attention"


def test_v51_public_mixed_continuation_identity_is_rejected() -> None:
    records = deepcopy(_records())
    row = records["continuation_branch_predictor_corrector_loadpath_turningpoint_owner_identity"]
    row.update({"result_branch_id": "branch:v51-secondary", "result_predictor_corrector_states": ["predictor:41", "corrector:41"], "result_load_path": [0.0, 0.35, 0.7, 0.88], "result_turning_point_index": 2, "result_solution_owner": "solution:continuation-stale"})
    assert validate_public_v51_identity(records)["status"] == "needs_attention"


def test_v51_public_invalid_canonical_records_are_rejected() -> None:
    records = deepcopy(_records())
    records["eigenmode_frequency_normalization_phase_subspace_mesh_owner_identity"]["phase_anchor"] = {"dof": -1, "component": "real", "sign": "positive"}
    records["continuation_branch_predictor_corrector_loadpath_turningpoint_owner_identity"]["turning_point_index"] = 2
    assert validate_public_v51_identity(records)["status"] == "needs_attention"
