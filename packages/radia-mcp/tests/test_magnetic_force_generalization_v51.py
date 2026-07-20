from copy import deepcopy

from radia_mcp.radia_ngsolve.potential_bem_identity_v51 import BEM, POTENTIAL, validate_public_identity


PROMOTED_CASE_IDS = {
    "v51_public_scalar_vector_potential_gauge_domain_interface_trace_solution_owner_mismatch",
    "v51_public_bem_matrix_reciprocity_symmetry_panel_orientation_cache_revision_owner_mismatch",
}


def _identity() -> dict[str, object]:
    generation = "magnetic-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    scalar = ["air", "iron"]
    vector = ["coil", "magnet"]
    traces = [{"scalar_domain": "air", "vector_domain": "magnet", "trace": "tangential_continuity"}]
    return {
        POTENTIAL: {
            "generation": generation, "gauge_generation": generation, "domain_generation": generation,
            "interface_generation": generation, "trace_generation": generation, "owner_generation": generation,
            "result_generation": generation, "gauge": "coulomb", "result_gauge": "coulomb",
            "scalar_potential_domains": scalar, "result_scalar_potential_domains": scalar,
            "vector_potential_domains": vector, "result_vector_potential_domains": vector,
            "interface_traces": traces, "result_interface_traces": traces, "interface_trace_sha256": "1" * 64,
            "result_interface_trace_sha256": "1" * 64, "solution_owner": "solution:coupled-v51",
            "result_solution_owner": "solution:coupled-v51", **result,
        },
        BEM: {
            "generation": generation, "reciprocity_generation": generation, "symmetry_generation": generation,
            "orientation_generation": generation, "cache_generation": generation, "owner_generation": generation,
            "result_generation": generation, "matrix_shape": [128, 128], "result_matrix_shape": [128, 128],
            "reciprocity_relative_error": 2.0e-13, "result_reciprocity_relative_error": 2.0e-13,
            "symmetry_class": "symmetric", "result_symmetry_class": "symmetric", "panel_orientation": "outward",
            "result_panel_orientation": "outward", "panel_orientation_sha256": "2" * 64,
            "result_panel_orientation_sha256": "2" * 64, "cache_revision": "cache:bem-v51",
            "result_cache_revision": "cache:bem-v51", "matrix_owner": "matrix:bem-v51",
            "result_matrix_owner": "matrix:bem-v51", **result,
        },
    }


def test_v51_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v51_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[POTENTIAL].update({"result_gauge": "tree_cotree", "result_solution_owner": "solution:stale"})
    identity[BEM].update({"result_panel_orientation": "inward", "result_matrix_owner": "matrix:stale"})
    assert not all(validate_public_identity(identity).values())


def test_v51_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[POTENTIAL]["gauge"] = identity[POTENTIAL]["result_gauge"] = "tree_cotree"
    identity[BEM]["panel_orientation"] = identity[BEM]["result_panel_orientation"] = "inward"
    assert not all(validate_public_identity(identity).values())
