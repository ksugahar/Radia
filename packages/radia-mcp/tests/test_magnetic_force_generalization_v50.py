from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.bem_motion_identity_v50 import BEM, MOTION, validate_public_identity


PROMOTED_CASE_IDS = {
    "v50_public_bem_singular_quadrature_self_panel_nearfield_regularization_owner_mismatch",
    "v50_public_motion_emf_velocity_frame_conductor_path_integration_result_owner_mismatch",
}


def _identity() -> dict[str, object]:
    bem_generation = "bem-quadrature-v50-901"
    motion_generation = "motion-emf-v50-901"
    nearfield = {"distance_ratio": 0.15, "regularization": "adaptive-subdivision", "max_depth": 6}
    velocity = [12.0, 0.0, 0.0]
    path = [[0.0, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.2, 0.0]]
    return {
        BEM: {
            "generation": bem_generation, "quadrature_generation": bem_generation, "self_panel_generation": bem_generation,
            "nearfield_generation": bem_generation, "mesh_generation": bem_generation, "result_generation": bem_generation,
            "singular_quadrature": "duffy-order-8", "result_singular_quadrature": "duffy-order-8",
            "self_panel_treatment": "analytic-solid-angle", "result_self_panel_treatment": "analytic-solid-angle",
            "nearfield_regularization": nearfield, "result_nearfield_regularization": nearfield,
            "mesh_owner": "mesh:bem-v50-901", "result_mesh_owner": "mesh:bem-v50-901",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        MOTION: {
            "generation": motion_generation, "velocity_generation": motion_generation, "frame_generation": motion_generation,
            "path_generation": motion_generation, "direction_generation": motion_generation, "result_generation": motion_generation,
            "velocity_m_s": velocity, "result_velocity_m_s": velocity,
            "velocity_frame": "frame:global", "result_velocity_frame": "frame:global",
            "conductor_path_m": path, "result_conductor_path_m": path,
            "integration_direction": "path-forward", "result_integration_direction": "path-forward",
            "emf_result_owner": "emf:conductor-v50-901", "result_emf_owner": "emf:conductor-v50-901",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def test_v50_positive_bem_and_motional_emf_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v50_bem_quadrature_self_panel_nearfield_and_owner_drift_is_rejected() -> None:
    identity = deepcopy(_identity())
    identity[BEM]["result_singular_quadrature"] = "gauss-order-2"
    identity[BEM]["result_self_panel_treatment"] = "centroid-sample"
    identity[BEM]["result_nearfield_regularization"] = {"distance_ratio": 1.0, "regularization": "none", "max_depth": 0}
    identity[BEM]["result_mesh_owner"] = "mesh:foreign"
    assert validate_public_identity(identity)["magnetic_force_v50_bem_singular_self_nearfield_mesh_owner"] is False


def test_v50_motion_emf_velocity_frame_path_direction_and_owner_drift_is_rejected() -> None:
    identity = deepcopy(_identity())
    identity[MOTION]["result_velocity_m_s"] = [-12.0, 0.0, 0.0]
    identity[MOTION]["result_velocity_frame"] = "frame:body"
    identity[MOTION]["result_conductor_path_m"] = list(reversed(identity[MOTION]["conductor_path_m"]))
    identity[MOTION]["result_integration_direction"] = "path-reverse"
    identity[MOTION]["result_emf_owner"] = "emf:foreign"
    assert validate_public_identity(identity)["magnetic_force_v50_motion_emf_velocity_frame_path_direction_owner"] is False
