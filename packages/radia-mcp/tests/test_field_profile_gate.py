from __future__ import annotations

import json

from radia_mcp.radia_ngsolve.field_profile_gate import (
    dual_formulation_symmetric_field_profile_gate,
    nonlinear_magnetic_field_energy_parity_gate,
    nonlinear_magnetic_refinement_energy_gate,
    nonlinear_magnetic_spatial_evidence_gate,
)
from radia_mcp.radia_ngsolve.server import (
    nonlinear_magnetic_field_energy_parity_gate as mcp_nonlinear_parity_gate,
    nonlinear_magnetic_refinement_energy_gate as mcp_nonlinear_refinement_gate,
    nonlinear_magnetic_spatial_evidence_gate as mcp_nonlinear_gate,
)


def _summary():
    return {
        "component_id": "B_axis",
        "field_unit": "T",
        "axis_unit": "m",
        "sample_count": 201,
        "axis_min": -1.0,
        "axis_max": 1.0,
        "center_axis": 0.0,
        "center_value_a": -5.61756e-9,
        "center_value_b": -5.61627e-9,
        "profile_relative_l2_difference": 9.42e-4,
        "center_relative_difference": 2.31e-4,
        "symmetry_relative_a": 2.51e-4,
        "symmetry_relative_b": 2.63e-4,
    }


def test_dual_formulation_profile_gate_accepts_agreement_and_symmetry():
    result = dual_formulation_symmetric_field_profile_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_dual_formulation_profile_gate_rejects_single_point_only_or_asymmetry():
    bad = _summary()
    bad["sample_count"] = 1
    bad["profile_relative_l2_difference"] = 0.2
    bad["symmetry_relative_b"] = 0.3
    result = dual_formulation_symmetric_field_profile_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "sample_count_sufficient",
        "profile_formulations_agree",
        "formulation_b_is_symmetric",
    }


def test_dual_formulation_profile_gate_rejects_trivial_zero_field():
    bad = _summary()
    bad["center_value_a"] = 0.0
    bad["center_value_b"] = 0.0
    assert dual_formulation_symmetric_field_profile_gate(bad)["checks"]["center_field_nonzero"] is False


def _nonlinear_volume_summary():
    identity = {
        "observable_id": "iron_volume_B",
        "field_unit": "T",
        "coordinate_system": "right-handed Cartesian",
    }
    return {
        **identity,
        "nonlinear": True,
        "solver_converged": True,
        "linear_reference_only": False,
        "response_order": 1,
        "material_update_order": 0,
        "spatial_observable": "volume_integral",
        "integration_order": 8,
        "sample_count": 0,
        "volume_m3": 1.0e-6,
        "average_field_T": [0.0, 0.0, -1.50],
        "rms_magnitude_T": 1.62,
        "reference": {
            **identity,
            "average_field_T": [0.0, 0.0, -1.48],
            "rms_magnitude_T": 1.60,
        },
    }


def test_nonlinear_magnetic_spatial_gate_accepts_matched_volume_evidence():
    result = nonlinear_magnetic_spatial_evidence_gate(_nonlinear_volume_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert result["metrics"]["average_vector_relative_difference"] < 0.02


def test_nonlinear_magnetic_spatial_gate_rejects_point_only_and_order_mismatch():
    bad = _nonlinear_volume_summary()
    bad.update(
        response_order=2,
        material_update_order=0,
        spatial_observable="point",
        integration_order=0,
        sample_count=1,
    )
    result = nonlinear_magnetic_spatial_evidence_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "response_material_orders_compatible",
        "spatial_observable_supported",
        "spatial_sampling_sufficient",
    }


def test_nonlinear_magnetic_spatial_gate_rejects_unconverged_or_linear_only():
    bad = _nonlinear_volume_summary()
    bad["solver_converged"] = False
    bad["linear_reference_only"] = True
    result = nonlinear_magnetic_spatial_evidence_gate(bad)
    assert set(result["issues"]) >= {
        "solver_converged",
        "not_linear_only_evidence",
    }


def test_nonlinear_magnetic_spatial_gate_rejects_rms_disagreement_and_wraps_mcp():
    bad = _nonlinear_volume_summary()
    bad["rms_magnitude_T"] = 2.1
    result = nonlinear_magnetic_spatial_evidence_gate(bad)
    assert result["checks"]["rms_magnitude_matches_reference"] is False
    wrapped = json.loads(mcp_nonlinear_gate(json.dumps(bad)))
    assert wrapped["status"] == "needs_attention"


def _nonlinear_refinement_summary():
    def level(mesh_size, average_z, rms, energy):
        identity = {
            "material_domain": "air|iron",
            "coordinate_system": "right-handed Cartesian",
            "unit_system": "SI",
            "nonlinear_state_id": "bh-state-17",
            "mesh_topology_geometry_sha256": "1" * 64,
            "bh_table_sha256": "2" * 64,
            "material_state_identity_sha256": "3" * 64,
            "solution_sha256": "4" * 64,
            "refinement_parent_identity_sha256": "5" * 64,
        }
        return {
            "mesh_size_m": mesh_size,
            "solver_converged": True,
            "response_order": 2,
            "material_update_order": 1,
            "physical_relative_permeability_bounds": [1.0, 2100.0],
            "volume_m3": 1.0e-6,
            "average_field_T": [0.0, 0.0, average_z],
            "rms_magnitude_T": rms,
            "magnetic_energy_J": energy,
            "magnetic_coenergy_J": 1.1 * energy,
            "h_dot_b_integral_J": 2.1 * energy,
            "legendre_residual_relative": 0.0,
            "field_identity": {**identity, "unit": "T"},
            "energy_identity": {**identity, "unit": "J"},
        }

    return {
        "levels": [
            level(8.0e-3, -1.42, 1.66, 1.20e-3),
            level(4.0e-3, -1.48, 1.62, 1.24e-3),
            level(2.0e-3, -1.50, 1.61, 1.25e-3),
        ]
    }


def test_nonlinear_refinement_energy_gate_accepts_contracting_matched_ladder():
    result = nonlinear_magnetic_refinement_energy_gate(_nonlinear_refinement_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert result["pair_relative_changes"][1]["combined"] < result["pair_relative_changes"][0]["combined"]


def test_nonlinear_refinement_energy_gate_rejects_nonpositive_material_state():
    bad = _nonlinear_refinement_summary()
    bad["levels"][1]["physical_relative_permeability_bounds"] = [0.0, 2100.0]
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_levels_valid"] is False
    assert result["level_checks"][1]["physical_permeability_positive"] is False


def test_nonlinear_refinement_energy_gate_rejects_noncontracting_ladder():
    bad = _nonlinear_refinement_summary()
    bad["levels"][2]["average_field_T"] = [0.0, 0.0, -1.25]
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["checks"]["field_rms_energy_changes_contract"] is False


def test_nonlinear_refinement_energy_gate_rejects_field_energy_identity_mismatch_and_wraps_mcp():
    bad = _nonlinear_refinement_summary()
    bad["levels"][2]["energy_identity"]["nonlinear_state_id"] = "stale-state"
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["checks"]["all_levels_valid"] is False
    assert result["level_checks"][2]["field_energy_identity_matches"] is False
    wrapped = json.loads(mcp_nonlinear_refinement_gate(json.dumps(bad)))
    assert wrapped["status"] == "needs_attention"


def test_nonlinear_refinement_energy_gate_rejects_stale_solution_digest():
    bad = _nonlinear_refinement_summary()
    bad["levels"][1]["energy_identity"]["solution_sha256"] = "5" * 64
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["level_checks"][1]["field_energy_identity_matches"] is False


def test_nonlinear_refinement_energy_gate_rejects_stale_parent_lineage():
    bad = _nonlinear_refinement_summary()
    for identity_name in ("field_identity", "energy_identity"):
        bad["levels"][1][identity_name][
            "refinement_parent_identity_sha256"
        ] = "6" * 64
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["refinement_parent_identity_consistent"] is False
    assert result["level_checks"][1]["field_energy_identity_matches"] is True


def test_nonlinear_refinement_energy_gate_rejects_legendre_residual():
    bad = _nonlinear_refinement_summary()
    bad["levels"][2]["legendre_residual_relative"] = 2.0e-4
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["checks"]["all_levels_valid"] is False
    assert result["level_checks"][2]["legendre_energy_identity_satisfied"] is False


def _nonlinear_parity_summary():
    identity = {
        "geometry_identity_sha256": "1" * 64,
        "material_identity_sha256": "2" * 64,
        "excitation_identity_sha256": "3" * 64,
        "bh_table_sha256": "4" * 64,
        "constitutive_interpolation": "monotone_pchip",
        "constitutive_extrapolation": "vacuum_slope",
        "magnetic_anisotropy": "isotropic",
        "region_labels": ["iron"],
        "coordinate_system": "right-handed Cartesian",
        "unit_system": "SI",
        "analysis_kind": "magnetostatic",
        "case_index": 1,
        "time_semantics": "static",
    }
    return {
        "candidate": {
            "identity": dict(identity),
            "average_field_T": [0.0, 0.0, 1.02],
            "rms_magnitude_T": 1.12,
            "magnetic_energy_J": 2.04,
            "magnetic_coenergy_J": 2.97,
        },
        "reference": {
            "identity": dict(identity),
            "average_field_T": [0.0, 0.0, 1.0],
            "rms_magnitude_T": 1.1,
            "magnetic_energy_J": 2.0,
            "magnetic_coenergy_J": 3.0,
        },
    }


def test_nonlinear_field_energy_parity_gate_accepts_identity_matched_observables():
    summary = _nonlinear_parity_summary()
    result = nonlinear_magnetic_field_energy_parity_gate(summary)
    assert result["status"] == "ok"
    assert result["comparison_performed"] is True
    wrapped = json.loads(mcp_nonlinear_parity_gate(json.dumps(summary)))
    assert wrapped["status"] == "ok"


def test_nonlinear_field_energy_parity_gate_rejects_material_law_before_numbers():
    bad = _nonlinear_parity_summary()
    bad["reference"]["identity"]["constitutive_interpolation"] = "piecewise_linear"
    result = nonlinear_magnetic_field_energy_parity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["physical_identity_matches"] is False
    assert result["identity_checks"]["constitutive_interpolation"] is False
    assert result["comparison_performed"] is False
    assert all(value is None for value in result["relative_differences"].values())
