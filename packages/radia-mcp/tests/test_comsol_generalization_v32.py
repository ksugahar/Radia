from __future__ import annotations

from test_comsol_generalization_v31 import (
    _summary,
    _with_v31_transform_and_force_identity,
    gate,
)


_PROMOTED_CASE_IDS = (
    "v32_public_nonlinear_segregated_iteration_relaxation_residual_jacobian_continuation_mismatch",
    "v32_public_degenerate_eigenmode_subspace_phase_normalization_participation_mass_mismatch",
)


def _with_v32_nonlinear_and_eigenmode_identity(summary: dict) -> dict:
    summary = _with_v31_transform_and_force_identity(summary)
    generation = "nonlinear-segregated-closure-191"
    summary[
        "nonlinear_segregated_group_relaxation_residual_jacobian_continuation_mesh_result_generation_identity"
    ] = {
        "nonlinear_generation": generation,
        "segregated_group_generation": generation,
        "relaxation_generation": generation,
        "residual_generation": generation,
        "jacobian_generation": generation,
        "continuation_generation": generation,
        "mesh_generation": generation,
        "result_generation": generation,
        "segregated_group_order": ["magnetic", "thermal", "structural"],
        "result_segregated_group_order": ["magnetic", "thermal", "structural"],
        "relaxation_schedule": [0.5, 0.7, 1.0],
        "result_relaxation_schedule": [0.5, 0.7, 1.0],
        "residual_norms": [1.0e-2, 2.0e-5, 4.0e-9],
        "accepted_residual_norms": [1.0e-2, 2.0e-5, 4.0e-9],
        "residual_relative_tolerance": 1.0e-8,
        "continuation_parameter": 1.0,
        "accepted_continuation_parameter": 1.0,
        "continuation_unit": "1",
        "accepted_continuation_unit": "1",
        "jacobian_sha256": "1" * 64,
        "accepted_jacobian_sha256": "1" * 64,
        "nonlinear_mesh_sha256": "2" * 64,
        "accepted_nonlinear_mesh_sha256": "2" * 64,
        "nonlinear_solution_sha256": "3" * 64,
        "accepted_nonlinear_solution_sha256": "3" * 64,
    }
    generation = "degenerate-eigenmode-closure-191"
    summary[
        "degenerate_eigenmode_subspace_phase_normalization_participation_mass_mesh_owner_result_generation_identity"
    ] = {
        "mode_generation": generation,
        "subspace_generation": generation,
        "phase_generation": generation,
        "normalization_generation": generation,
        "participation_generation": generation,
        "mass_generation": generation,
        "mesh_generation": generation,
        "result_generation": generation,
        "eigenvalues": [100.0, 100.0],
        "result_eigenvalues": [100.0, 100.0],
        "degenerate_subspace_sha256": "4" * 64,
        "result_degenerate_subspace_sha256": "4" * 64,
        "phase_anchor_dofs": [12, 37],
        "result_phase_anchor_dofs": [12, 37],
        "normalization": "mass_orthonormal",
        "result_normalization": "mass_orthonormal",
        "participation_factors": [0.8, 0.6],
        "result_participation_factors": [0.8, 0.6],
        "effective_masses_kg": [0.64, 0.36],
        "result_effective_masses_kg": [0.64, 0.36],
        "mode_owner": "component1/solid/eig1",
        "result_mode_owner": "component1/solid/eig1",
        "eigenmode_mesh_sha256": "5" * 64,
        "result_eigenmode_mesh_sha256": "5" * 64,
        "eigenmode_result_sha256": "6" * 64,
        "accepted_eigenmode_result_sha256": "6" * 64,
    }
    return summary


def test_v32_public_positive_nonlinear_and_eigenmode_contracts() -> None:
    result = gate(_with_v32_nonlinear_and_eigenmode_identity(_summary()))
    assert result["status"] == "ok"


def test_v32_public_nonlinear_segregated_iteration_relaxation_residual_jacobian_continuation_mismatch() -> None:
    summary = _with_v32_nonlinear_and_eigenmode_identity(_summary())
    summary[
        "nonlinear_segregated_group_relaxation_residual_jacobian_continuation_mesh_result_generation_identity"
    ].update(
        {
            "segregated_group_generation": "nonlinear-segregated-190",
            "jacobian_generation": "nonlinear-segregated-189",
            "result_generation": "nonlinear-segregated-188",
            "result_segregated_group_order": ["thermal", "magnetic", "structural"],
            "result_relaxation_schedule": [1.0, 1.0, 1.0],
            "accepted_residual_norms": [1.0e-2, 4.0e-4, 9.0e-5],
            "accepted_continuation_parameter": 0.8,
            "accepted_continuation_unit": "percent",
            "accepted_jacobian_sha256": "c" * 64,
            "accepted_nonlinear_mesh_sha256": "d" * 64,
            "accepted_nonlinear_solution_sha256": "e" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_segregated_solutions_use_current_groups_relaxation_residual_jacobian_continuation_mesh_and_result"
    ]


def test_v32_public_degenerate_eigenmode_subspace_phase_normalization_participation_mass_mismatch() -> None:
    summary = _with_v32_nonlinear_and_eigenmode_identity(_summary())
    summary[
        "degenerate_eigenmode_subspace_phase_normalization_participation_mass_mesh_owner_result_generation_identity"
    ].update(
        {
            "subspace_generation": "degenerate-eigenmode-190",
            "normalization_generation": "degenerate-eigenmode-189",
            "result_generation": "degenerate-eigenmode-188",
            "result_eigenvalues": [100.0, 101.0],
            "result_degenerate_subspace_sha256": "f" * 64,
            "result_phase_anchor_dofs": [37, 12],
            "result_normalization": "max_component",
            "result_participation_factors": [0.6, -0.8],
            "result_effective_masses_kg": [0.36, 0.80],
            "result_mode_owner": "component2/acpr/eig1",
            "result_eigenmode_mesh_sha256": "0" * 64,
            "accepted_eigenmode_result_sha256": "1" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "degenerate_eigenmodes_use_current_subspace_phase_normalization_participation_mass_mesh_owner_and_result"
    ]
