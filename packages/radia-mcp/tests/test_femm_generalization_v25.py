from __future__ import annotations

from test_femm_generalization_v19 import _gate
from test_femm_generalization_v24 import _identity_v24
from test_force_coenergy_gate import _quadratic_case


def _identity_v25(sample_count):
    identity = _identity_v24(sample_count)
    identity[
        "nonlinear_bh_incremental_energy_coenergy_force_branch_mesh_generation_identity"
    ] = {
        "nonlinear_generation": "nonlinear-111",
        "bh_curve_nonlinear_generation": "nonlinear-111",
        "material_nonlinear_generation": "nonlinear-111",
        "branch_nonlinear_generation": "nonlinear-111",
        "incremental_state_nonlinear_generation": "nonlinear-111",
        "mesh_nonlinear_generation": "nonlinear-111",
        "energy_nonlinear_generation": "nonlinear-111",
        "coenergy_nonlinear_generation": "nonlinear-111",
        "force_nonlinear_generation": "nonlinear-111",
        "result_nonlinear_generation": "nonlinear-111",
        "nonlinear_material_ids": ["core"],
        "result_nonlinear_material_ids": ["core"],
        "bh_curve_sha256": "1" * 64,
        "result_bh_curve_sha256": "1" * 64,
        "material_map_sha256": "2" * 64,
        "result_material_map_sha256": "2" * 64,
        "branch_id": "ascending:step-12",
        "result_branch_id": "ascending:step-12",
        "load_current_a": 8.0,
        "result_load_current_a": 8.0,
        "incremental_state_sha256": "3" * 64,
        "result_incremental_state_sha256": "3" * 64,
        "mesh_sha256": "4" * 64,
        "result_mesh_sha256": "4" * 64,
        "magnetic_energy_j": 0.42,
        "result_magnetic_energy_j": 0.42,
        "magnetic_coenergy_j": 0.58,
        "result_magnetic_coenergy_j": 0.58,
        "incremental_force_n": [12.2, -0.1],
        "result_incremental_force_n": [12.2, -0.1],
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    identity[
        "open_boundary_domain_decay_multipole_moment_material_generation_identity"
    ] = {
        "boundary_generation": "open-boundary-111",
        "domain_boundary_generation": "open-boundary-111",
        "mesh_boundary_generation": "open-boundary-111",
        "material_boundary_generation": "open-boundary-111",
        "multipole_boundary_generation": "open-boundary-111",
        "decay_boundary_generation": "open-boundary-111",
        "result_boundary_generation": "open-boundary-111",
        "boundary_type": "asymptotic_multipole",
        "result_boundary_type": "asymptotic_multipole",
        "source_radius_m": 0.05,
        "result_source_radius_m": 0.05,
        "outer_radius_m": 0.4,
        "result_outer_radius_m": 0.4,
        "multipole_order": 3,
        "result_multipole_order": 3,
        "material_map_sha256": "6" * 64,
        "result_material_map_sha256": "6" * 64,
        "mesh_sha256": "7" * 64,
        "result_mesh_sha256": "7" * 64,
        "multipole_moment_sha256": "8" * 64,
        "result_multipole_moment_sha256": "8" * 64,
        "decay_sample_radii_m": [0.2, 0.3, 0.4],
        "result_decay_sample_radii_m": [0.2, 0.3, 0.4],
        "decay_flux_density_t": [1.0e-3, 3.0e-4, 1.2e-4],
        "result_decay_flux_density_t": [1.0e-3, 3.0e-4, 1.2e-4],
        "result_sha256": "9" * 64,
        "accepted_result_sha256": "9" * 64,
    }
    return identity


def test_v25_public_positive_nonlinear_and_open_boundary_identity():
    positions, _, _ = _quadratic_case()
    assert _gate(_identity_v25(len(positions)))["status"] == "ok"


def test_v25_public_nonlinear_bh_incremental_energy_coenergy_force_branch_mesh_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v25(len(positions))
    identity[
        "nonlinear_bh_incremental_energy_coenergy_force_branch_mesh_generation_identity"
    ].update(
        {
            "bh_curve_nonlinear_generation": "nonlinear-110",
            "material_nonlinear_generation": "nonlinear-109",
            "branch_nonlinear_generation": "nonlinear-108",
            "incremental_state_nonlinear_generation": "nonlinear-107",
            "mesh_nonlinear_generation": "nonlinear-106",
            "energy_nonlinear_generation": "nonlinear-105",
            "coenergy_nonlinear_generation": "nonlinear-104",
            "force_nonlinear_generation": "nonlinear-103",
            "result_nonlinear_material_ids": ["core-old"],
            "result_bh_curve_sha256": "d" * 64,
            "result_material_map_sha256": "e" * 64,
            "result_branch_id": "descending:step-5",
            "result_load_current_a": 6.0,
            "result_incremental_state_sha256": "f" * 64,
            "result_mesh_sha256": "0" * 64,
            "result_magnetic_energy_j": 0.31,
            "result_magnetic_coenergy_j": 0.45,
            "result_incremental_force_n": [9.8, 0.2],
            "accepted_result_sha256": "1" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_bh_incremental_force_uses_current_branch_material_mesh_energy_and_coenergy"
    ]


def test_v25_public_open_boundary_domain_decay_multipole_moment_material_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v25(len(positions))
    identity[
        "open_boundary_domain_decay_multipole_moment_material_generation_identity"
    ].update(
        {
            "domain_boundary_generation": "open-boundary-110",
            "mesh_boundary_generation": "open-boundary-109",
            "material_boundary_generation": "open-boundary-108",
            "multipole_boundary_generation": "open-boundary-107",
            "decay_boundary_generation": "open-boundary-106",
            "result_boundary_type": "dirichlet_zero",
            "result_source_radius_m": 0.08,
            "result_outer_radius_m": 0.15,
            "result_multipole_order": 1,
            "result_material_map_sha256": "2" * 64,
            "result_mesh_sha256": "3" * 64,
            "result_multipole_moment_sha256": "4" * 64,
            "result_decay_sample_radii_m": [0.1, 0.12, 0.15],
            "result_decay_flux_density_t": [1.0e-3, 1.2e-3, 1.1e-3],
            "accepted_result_sha256": "5" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "open_boundary_uses_current_domain_decay_multipole_material_and_mesh"
    ]
