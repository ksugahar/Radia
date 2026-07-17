from __future__ import annotations

import math

from test_femm_generalization_v19 import _gate
from test_femm_generalization_v26 import _identity_v26
from test_force_coenergy_gate import _quadratic_case


def _identity_v27(sample_count):
    identity = _identity_v26(sample_count)
    identity["nonlinear_bh_minor_loop_branch_interpolation_state_coenergy_force_generation_identity"] = {
        "nonlinear_generation": "bh-minor-loop-141",
        "branch_nonlinear_generation": "bh-minor-loop-141",
        "interpolation_nonlinear_generation": "bh-minor-loop-141",
        "state_nonlinear_generation": "bh-minor-loop-141",
        "coenergy_nonlinear_generation": "bh-minor-loop-141",
        "mesh_nonlinear_generation": "bh-minor-loop-141",
        "force_nonlinear_generation": "bh-minor-loop-141",
        "result_nonlinear_generation": "bh-minor-loop-141",
        "bh_branch": "ascending-minor-loop",
        "result_bh_branch": "ascending-minor-loop",
        "interpolation_rule": "monotone-cubic-h",
        "result_interpolation_rule": "monotone-cubic-h",
        "state_point_am": [120.0, 800.0],
        "result_state_point_am": [120.0, 800.0],
        "magnetic_coenergy_j": 0.37,
        "result_magnetic_coenergy_j": 0.37,
        "force_n": [18.0, -0.4],
        "result_force_n": [18.0, -0.4],
        "bh_table_sha256": "1" * 64,
        "result_bh_table_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "solution_sha256": "3" * 64,
        "accepted_solution_sha256": "3" * 64,
    }
    mu0 = 4.0e-7 * math.pi
    sigma = 5.8e7
    frequency = 2000.0
    skin_depth = math.sqrt(2.0 / (2.0 * math.pi * frequency * mu0 * sigma))
    identity["harmonic_eddy_phasor_conductivity_skin_depth_frequency_loss_mesh_generation_identity"] = {
        "eddy_generation": "harmonic-eddy-141",
        "phasor_eddy_generation": "harmonic-eddy-141",
        "conductivity_eddy_generation": "harmonic-eddy-141",
        "skin_eddy_generation": "harmonic-eddy-141",
        "frequency_eddy_generation": "harmonic-eddy-141",
        "loss_eddy_generation": "harmonic-eddy-141",
        "mesh_eddy_generation": "harmonic-eddy-141",
        "result_eddy_generation": "harmonic-eddy-141",
        "phasor_convention": "exp(+jwt)",
        "result_phasor_convention": "exp(+jwt)",
        "conductivity_s_m": sigma,
        "result_conductivity_s_m": sigma,
        "relative_permeability": 1.0,
        "result_relative_permeability": 1.0,
        "frequency_hz": frequency,
        "result_frequency_hz": frequency,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "minimum_elements_per_skin_depth": 4.0,
        "result_minimum_elements_per_skin_depth": 4.0,
        "joule_loss_w": 42.0,
        "result_joule_loss_w": 42.0,
        "mesh_sha256": "4" * 64,
        "result_mesh_sha256": "4" * 64,
        "loss_result_sha256": "5" * 64,
        "accepted_loss_result_sha256": "5" * 64,
    }
    return identity


def test_v27_public_positive_nonlinear_and_harmonic_eddy_identity():
    positions, _, _ = _quadratic_case()
    assert _gate(_identity_v27(len(positions)))["status"] == "ok"


def test_v27_public_nonlinear_bh_minor_loop_branch_interpolation_state_coenergy_force_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v27(len(positions))
    identity[
        "nonlinear_bh_minor_loop_branch_interpolation_state_coenergy_force_generation_identity"
    ].update({
        "branch_nonlinear_generation": "bh-minor-loop-140",
        "state_nonlinear_generation": "bh-minor-loop-139",
        "result_bh_branch": "descending-major-loop",
        "result_interpolation_rule": "linear-b",
        "result_state_point_am": [-120.0, 700.0],
        "result_magnetic_coenergy_j": 0.28,
        "result_force_n": [12.0, 1.0],
        "result_bh_table_sha256": "a" * 64,
        "result_mesh_sha256": "b" * 64,
        "accepted_solution_sha256": "c" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_force_uses_current_minor_loop_branch_interpolation_state_coenergy_and_mesh"
    ]


def test_v27_public_harmonic_eddy_phasor_convention_conductivity_skin_depth_loss_mesh_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v27(len(positions))
    identity[
        "harmonic_eddy_phasor_conductivity_skin_depth_frequency_loss_mesh_generation_identity"
    ].update({
        "phasor_eddy_generation": "harmonic-eddy-140",
        "conductivity_eddy_generation": "harmonic-eddy-139",
        "mesh_eddy_generation": "harmonic-eddy-138",
        "result_phasor_convention": "exp(-jwt)",
        "result_conductivity_s_m": 5.8e4,
        "result_frequency_hz": 50.0,
        "result_skin_depth_m": 0.02,
        "result_minimum_elements_per_skin_depth": 0.5,
        "result_joule_loss_w": 4.2,
        "result_mesh_sha256": "d" * 64,
        "accepted_loss_result_sha256": "e" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "harmonic_eddy_loss_uses_current_phasor_conductivity_skin_depth_frequency_and_mesh"
    ]
