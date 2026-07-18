from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v31_profile import _summary_v31


_PROMOTED_CASE_IDS = (
    "v32_public_maglev_equilibrium_force_displacement_stiffness_derivative_sign_mesh_mismatch",
    "v32_public_bem_surface_charge_net_zero_gauge_reference_energy_reciprocity_mismatch",
)


def _summary_v32():
    summary = _summary_v31()
    identity = summary["artifact_identity"]
    generation = "maglev-equilibrium-361"
    identity[
        "maglev_equilibrium_force_displacement_derivative_stiffness_gravity_mesh_result_identity"
    ] = {
        "equilibrium_generation": generation,
        **{
            key: generation
            for key in (
                "force_equilibrium_generation",
                "frame_equilibrium_generation",
                "displacement_equilibrium_generation",
                "derivative_equilibrium_generation",
                "stiffness_equilibrium_generation",
                "gravity_equilibrium_generation",
                "mesh_equilibrium_generation",
                "result_equilibrium_generation",
            )
        },
        "force_sign_convention": "positive_up",
        "result_force_sign_convention": "positive_up",
        "displacement_frame": "global_z_up",
        "result_displacement_frame": "global_z_up",
        "displacement_samples_m": [-1.0e-4, 0.0, 1.0e-4],
        "result_displacement_samples_m": [-1.0e-4, 0.0, 1.0e-4],
        "magnetic_force_samples_n": [10.001, 9.81, 9.619],
        "result_magnetic_force_samples_n": [10.001, 9.81, 9.619],
        "derivative_stencil": "symmetric_central_difference",
        "result_derivative_stencil": "symmetric_central_difference",
        "force_derivative_n_per_m": -1910.0,
        "result_force_derivative_n_per_m": -1910.0,
        "vertical_stiffness_n_per_m": 1910.0,
        "result_vertical_stiffness_n_per_m": 1910.0,
        "gravity_force_n": -9.81,
        "result_gravity_force_n": -9.81,
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_owner": "maglev/case-361/equilibrium-z",
        "accepted_result_owner": "maglev/case-361/equilibrium-z",
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "bem-charge-361"
    identity[
        "bem_surface_charge_gauge_normal_energy_reciprocity_geometry_owner_result_identity"
    ] = {
        "bem_generation": generation,
        **{
            key: generation
            for key in (
                "charge_bem_generation",
                "gauge_bem_generation",
                "normal_bem_generation",
                "energy_bem_generation",
                "reciprocity_bem_generation",
                "geometry_bem_generation",
                "owner_bem_generation",
                "result_bem_generation",
            )
        },
        "net_surface_charge": 0.0,
        "result_net_surface_charge": 0.0,
        "charge_balance_tolerance": 1.0e-12,
        "gauge_reference": "mean_zero_scalar_potential",
        "result_gauge_reference": "mean_zero_scalar_potential",
        "source_normal": [0.0, 0.0, 1.0],
        "result_source_normal": [0.0, 0.0, 1.0],
        "target_normal": [0.0, 0.0, -1.0],
        "result_target_normal": [0.0, 0.0, -1.0],
        "field_energy_j": 0.25,
        "result_field_energy_j": 0.25,
        "reciprocity_residual": 1.0e-12,
        "result_reciprocity_residual": 1.0e-12,
        "reciprocity_tolerance": 1.0e-9,
        "geometry_sha256": "3" * 64,
        "result_geometry_sha256": "3" * 64,
        "result_owner": "bem/case-361/scalar-potential",
        "accepted_result_owner": "bem/case-361/scalar-potential",
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v32_public_positive_maglev_equilibrium_and_bem_surface_charge():
    assert magnetic_force_method_profile_gate(_summary_v32())["status"] == "ok"


def test_v32_public_maglev_equilibrium_force_displacement_stiffness_derivative_sign_mesh_mismatch():
    summary = _summary_v32()
    record = summary["artifact_identity"][
        "maglev_equilibrium_force_displacement_derivative_stiffness_gravity_mesh_result_identity"
    ]
    record.update(
        {
            "force_equilibrium_generation": "maglev-equilibrium-360",
            "mesh_equilibrium_generation": "maglev-equilibrium-359",
            "result_equilibrium_generation": "maglev-equilibrium-358",
            "result_force_sign_convention": "positive_down",
            "result_displacement_frame": "local_y_down",
            "result_displacement_samples_m": [0.0, 1.0e-4, 2.0e-4],
            "result_magnetic_force_samples_n": [-9.81, -9.7, -9.5],
            "result_derivative_stencil": "forward_difference",
            "result_force_derivative_n_per_m": 1100.0,
            "result_vertical_stiffness_n_per_m": -1100.0,
            "result_gravity_force_n": 9.81,
            "result_mesh_sha256": "7" * 64,
            "accepted_result_owner": "maglev/old-case",
            "accepted_result_sha256": "8" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "maglev_equilibrium_uses_upward_force_global_displacement_central_stiffness_gravity_mesh_and_result"
    ]


def test_v32_public_bem_surface_charge_net_zero_gauge_reference_energy_reciprocity_mismatch():
    summary = _summary_v32()
    record = summary["artifact_identity"][
        "bem_surface_charge_gauge_normal_energy_reciprocity_geometry_owner_result_identity"
    ]
    record.update(
        {
            "charge_bem_generation": "bem-charge-360",
            "geometry_bem_generation": "bem-charge-359",
            "result_bem_generation": "bem-charge-358",
            "result_net_surface_charge": 0.1,
            "result_gauge_reference": "pin_first_node",
            "result_source_normal": [0.0, 0.0, -1.0],
            "result_target_normal": [0.0, 0.0, -1.0],
            "result_field_energy_j": -0.25,
            "result_reciprocity_residual": 0.2,
            "result_geometry_sha256": "9" * 64,
            "accepted_result_owner": "bem/old-case",
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_surface_charge_uses_neutral_charge_mean_zero_gauge_opposed_normals_energy_reciprocity_geometry_and_result"
    ]
