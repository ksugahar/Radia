from copy import deepcopy

from radia_mcp.radia_ngsolve.energy_derivative_identity_v53 import MAGLEV, QUADRATURE, validate_public_identity


PROMOTED_CASE_IDS = {
    "v53_public_bem_singular_quadrature_self_near_far_panel_owner_mismatch",
    "v53_public_maglev_equilibrium_force_stiffness_derivative_displacement_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity():
    interactions = [{"source_panel": 11, "target_panel": 11, "separation_over_size": 0.0, "classification": "self", "rule": "duffy_singular"}, {"source_panel": 11, "target_panel": 12, "separation_over_size": 0.4, "classification": "near", "rule": "adaptive_near"}, {"source_panel": 11, "target_panel": 91, "separation_over_size": 4.5, "classification": "far", "rule": "gauss_far"}]
    quadrature = {**_generations("quad-v53", ("classification_generation", "quadrature_generation", "panel_generation", "owner_generation", "result_generation")), "panel_interactions": interactions, "result_panel_interactions": interactions, "panel_owner": "panel-set:v53", "result_panel_owner": "panel-set:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    displacement = [-0.001, 0.0, 0.001]; force = [1.2, 0.0, -1.2]
    maglev = {**_generations("maglev-v53", ("equilibrium_generation", "force_generation", "stiffness_generation", "displacement_generation", "owner_generation", "result_generation")), "displacement_path_m": displacement, "result_displacement_path_m": displacement, "force_path_n": force, "result_force_path_n": force, "equilibrium_index": 1, "result_equilibrium_index": 1, "equilibrium_force_n": 0.0, "result_equilibrium_force_n": 0.0, "stiffness_n_per_m": -1200.0, "result_stiffness_n_per_m": -1200.0, "body_owner": "body:v53", "result_body_owner": "body:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {QUADRATURE: quadrature, MAGLEV: maglev}


def test_v53_positive_public_artifacts_are_accepted():
    assert all(validate_public_identity(_identity()).values())


def test_v53_frozen_counterfactuals_are_rejected():
    identity = deepcopy(_identity())
    identity[QUADRATURE]["result_panel_owner"] = "panel-set:stale"
    identity[MAGLEV]["result_stiffness_n_per_m"] = 1200.0
    assert not all(validate_public_identity(identity).values())


def test_v53_self_consistent_wrong_physics_is_rejected():
    identity = deepcopy(_identity())
    identity[QUADRATURE]["panel_interactions"][0]["classification"] = identity[QUADRATURE]["result_panel_interactions"][0]["classification"] = "far"
    identity[MAGLEV]["stiffness_n_per_m"] = identity[MAGLEV]["result_stiffness_n_per_m"] = 1200.0
    assert not all(validate_public_identity(identity).values())
