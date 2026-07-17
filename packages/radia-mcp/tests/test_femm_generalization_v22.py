from __future__ import annotations

from test_femm_generalization_v21 import _identity_v21
from test_femm_generalization_v19 import _gate
from test_force_coenergy_gate import _quadratic_case


def _identity_v22(sample_count):
    identity = _identity_v21(sample_count)
    identity["axisymmetric_force_energy_measure_depth_coordinate_generation_identity"] = {
        "solve_generation": "axisym-solve-41",
        "measure_solve_generation": "axisym-solve-41",
        "depth_solve_generation": "axisym-solve-41",
        "coordinate_solve_generation": "axisym-solve-41",
        "result_solve_generation": "axisym-solve-41",
        "problem_type": "axisymmetric",
        "result_problem_type": "axisymmetric",
        "measure_convention": "2*pi*r",
        "result_measure_convention": "2*pi*r",
        "planar_depth_m": 1.0,
        "result_planar_depth_m": 1.0,
        "coordinate_convention": "r_z",
        "result_coordinate_convention": "r_z",
        "force_energy_values": [12.5, 0.031],
        "reported_force_energy_values": [12.5, 0.031],
        "normalization_table_sha256": "1" * 64,
        "result_normalization_table_sha256": "1" * 64,
    }
    identity["nonlinear_incremental_mu_force_branch_perturbation_generation_identity"] = {
        "operating_point_generation": "nonlinear-op-41",
        "branch_operating_point_generation": "nonlinear-op-41",
        "differential_mu_operating_point_generation": "nonlinear-op-41",
        "perturbation_operating_point_generation": "nonlinear-op-41",
        "force_operating_point_generation": "nonlinear-op-41",
        "branch_id": "up-sweep:17",
        "force_branch_id": "up-sweep:17",
        "perturbation_current_a": 0.01,
        "force_perturbation_current_a": 0.01,
        "differential_mu_sha256": "2" * 64,
        "force_differential_mu_sha256": "2" * 64,
        "incremental_force_n": [0.12, -0.03],
        "reported_incremental_force_n": [0.12, -0.03],
        "incremental_state_sha256": "3" * 64,
        "force_incremental_state_sha256": "3" * 64,
    }
    return identity


def test_v22_public_positive_axisymmetric_and_incremental_force_identity():
    positions, _, _ = _quadratic_case()
    assert _gate(_identity_v22(len(positions)))["status"] == "ok"


def test_v22_public_axisymmetric_force_energy_two_pi_r_depth_normalization_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v22(len(positions))
    identity[
        "axisymmetric_force_energy_measure_depth_coordinate_generation_identity"
    ].update(
        {
            "measure_solve_generation": "axisym-solve-40",
            "depth_solve_generation": "axisym-solve-39",
            "coordinate_solve_generation": "axisym-solve-38",
            "result_problem_type": "planar",
            "result_measure_convention": "planar_depth",
            "result_planar_depth_m": 0.05,
            "result_coordinate_convention": "x_y",
            "reported_force_energy_values": [0.625, 0.00155],
            "result_normalization_table_sha256": "9" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_force_energy_uses_current_measure_depth_and_coordinates"
    ]


def test_v22_public_nonlinear_incremental_permeability_force_perturbation_branch_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v22(len(positions))
    identity[
        "nonlinear_incremental_mu_force_branch_perturbation_generation_identity"
    ].update(
        {
            "branch_operating_point_generation": "nonlinear-op-40",
            "differential_mu_operating_point_generation": "nonlinear-op-39",
            "perturbation_operating_point_generation": "nonlinear-op-38",
            "force_branch_id": "down-sweep:17",
            "force_perturbation_current_a": 0.1,
            "force_differential_mu_sha256": "a" * 64,
            "reported_incremental_force_n": [-0.08, 0.04],
            "force_incremental_state_sha256": "b" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_incremental_force_uses_current_branch_mu_and_perturbation"
    ]
