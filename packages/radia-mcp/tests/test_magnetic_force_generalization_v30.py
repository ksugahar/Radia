from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.force_coenergy_gate import force_coenergy_displacement_gate
from test_force_coenergy_gate import _quadratic_case
from test_magnetic_force_generalization_v29 import _force_summary_v29


_PROMOTED_CASE_IDS = (
    "v30_public_nonlinear_coenergy_force_current_perturbation_remesh_central_difference_mismatch",
    "v30_public_axisymmetric_force_r_weight_jacobian_coordinate_stress_contour_mismatch",
)


def _identity_v30():
    identity = _force_summary_v29()
    generation = "nonlinear-coenergy-171"
    dx = 1.0e-4
    coenergy = [2.005, 2.0, 1.995]
    force = -(coenergy[2] - coenergy[0]) / (2.0 * dx)
    identity[
        "nonlinear_coenergy_force_current_perturbation_remesh_central_difference_frame_result_identity"
    ] = {
        "force_generation": generation,
        "current_force_generation": generation,
        "coenergy_force_generation": generation,
        "remesh_force_generation": generation,
        "difference_force_generation": generation,
        "frame_force_generation": generation,
        "result_force_generation": generation,
        "nonlinear_material": True,
        "result_nonlinear_material": True,
        "current_constraint": "fixed_current",
        "result_current_constraint": "fixed_current",
        "nominal_current_a": 10.0,
        "branch_currents_a": [10.0, 10.0, 10.0],
        "result_branch_currents_a": [10.0, 10.0, 10.0],
        "displacements_m": [-dx, 0.0, dx],
        "result_displacements_m": [-dx, 0.0, dx],
        "coenergy_j": coenergy,
        "result_coenergy_j": coenergy,
        "difference_rule": "negative_central_difference",
        "result_difference_rule": "negative_central_difference",
        "branch_mesh_generations": ["mesh-minus-171", "mesh-center-171", "mesh-plus-171"],
        "result_branch_mesh_generations": ["mesh-minus-171", "mesh-center-171", "mesh-plus-171"],
        "displacement_frame": "global_x",
        "result_displacement_frame": "global_x",
        "force_n": force,
        "result_force_n": force,
        "branch_result_sha256": "1" * 64,
        "accepted_branch_result_sha256": "1" * 64,
    }
    generation = "axisym-stress-171"
    radius, jacobian, stress = 0.025, 0.001, 1200.0
    radial_weight = 2.0 * math.pi * radius
    force = stress * radial_weight * jacobian
    identity[
        "axisymmetric_force_radial_weight_jacobian_coordinate_stress_contour_material_mesh_result_identity"
    ] = {
        "axisymmetric_generation": generation,
        "radial_weight_generation": generation,
        "jacobian_generation": generation,
        "coordinate_generation": generation,
        "stress_contour_generation": generation,
        "material_side_generation": generation,
        "mesh_generation": generation,
        "result_generation": generation,
        "coordinate_convention": "r_z_axisymmetric",
        "result_coordinate_convention": "r_z_axisymmetric",
        "radius_m": radius,
        "result_radius_m": radius,
        "radial_weight": radial_weight,
        "result_radial_weight": radial_weight,
        "line_jacobian_m": jacobian,
        "result_line_jacobian_m": jacobian,
        "stress_normal_pa": stress,
        "result_stress_normal_pa": stress,
        "stress_contour_closed": True,
        "result_stress_contour_closed": True,
        "stress_contour_orientation": "counterclockwise",
        "result_stress_contour_orientation": "counterclockwise",
        "stress_contour_material_side": "air",
        "result_stress_contour_material_side": "air",
        "force_n": force,
        "result_force_n": force,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "force_result_sha256": "3" * 64,
        "accepted_force_result_sha256": "3" * 64,
    }
    return identity


def _gate(identity):
    positions, coenergy, forces = _quadratic_case()
    return force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )


def test_v30_public_positive_nonlinear_coenergy_and_axisymmetric_stress():
    assert _gate(_identity_v30())["status"] == "ok"


def test_v30_public_nonlinear_coenergy_force_current_perturbation_remesh_central_difference_mismatch():
    identity = _identity_v30()
    identity[
        "nonlinear_coenergy_force_current_perturbation_remesh_central_difference_frame_result_identity"
    ].update({
        "current_force_generation": "nonlinear-coenergy-170",
        "remesh_force_generation": "nonlinear-coenergy-169",
        "result_force_generation": "nonlinear-coenergy-168",
        "result_current_constraint": "fixed_flux",
        "result_branch_currents_a": [9.0, 10.0, 11.0],
        "result_displacements_m": [-2.0e-4, 0.0, 1.0e-4],
        "result_coenergy_j": [2.005, 2.0, 2.02],
        "result_difference_rule": "forward_difference",
        "result_branch_mesh_generations": ["mesh-old", "mesh-center-171", "mesh-plus-171"],
        "result_displacement_frame": "local_r",
        "result_force_n": -150.0,
        "accepted_branch_result_sha256": "8" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_coenergy_force_uses_fixed_current_symmetric_displacement_remesh_frame_and_result"
    ]


def test_v30_public_axisymmetric_force_r_weight_jacobian_coordinate_stress_contour_mismatch():
    identity = _identity_v30()
    identity[
        "axisymmetric_force_radial_weight_jacobian_coordinate_stress_contour_material_mesh_result_identity"
    ].update({
        "radial_weight_generation": "axisym-stress-170",
        "stress_contour_generation": "axisym-stress-169",
        "result_generation": "axisym-stress-168",
        "result_coordinate_convention": "x_y_planar",
        "result_radius_m": 0.05,
        "result_radial_weight": 1.0,
        "result_line_jacobian_m": 0.002,
        "result_stress_normal_pa": -1200.0,
        "result_stress_contour_closed": False,
        "result_stress_contour_orientation": "clockwise",
        "result_stress_contour_material_side": "iron",
        "result_force_n": 12.0,
        "result_mesh_sha256": "9" * 64,
        "accepted_force_result_sha256": "a" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_stress_force_uses_two_pi_r_jacobian_air_contour_mesh_and_result"
    ]
