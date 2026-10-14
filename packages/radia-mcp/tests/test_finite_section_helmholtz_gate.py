import copy
import asyncio
import json

import pytest

from radia_mcp.radia_ngsolve.finite_section_helmholtz_gate import (
    divergence_free_source_assembly_gate,
    finite_section_helmholtz_axis_bz,
    finite_section_helmholtz_mixed_omega_gate,
)
from radia_mcp.radia_ngsolve.server import (
    divergence_free_source_assembly_gate as mcp_source_gate,
    finite_section_helmholtz_mixed_omega_gate as mcp_gate,
    mcp,
)


def test_accepts_oriented_divergence_free_surface_source_assembly():
    result = divergence_free_source_assembly_gate({
        "schema": "radia.divergence-free-source-assembly.v1",
        "scope": "affine_p1_constant_mu",
        "source_identity": "divergence_free",
        "boundary_orientation": "opposite_vertex",
        "manufactured_load_relative_error": 4e-15,
        "solved_field_relative_difference": 4.7e-4,
        "allowed_solved_field_relative_difference": 1e-3,
        "volume_rhs_seconds": 7.9,
        "surface_rhs_seconds": 0.55,
    })
    assert result["status"] == "ok"
    assert result["speedup"] > 10
    assert json.loads(mcp_source_gate(json.dumps({
        "schema": "radia.divergence-free-source-assembly.v1",
        "scope": "affine_p1_constant_mu",
        "source_identity": "divergence_free",
        "boundary_orientation": "opposite_vertex",
        "manufactured_load_relative_error": 4e-15,
        "solved_field_relative_difference": 4.7e-4,
        "allowed_solved_field_relative_difference": 1e-3,
        "volume_rhs_seconds": 7.9,
        "surface_rhs_seconds": 0.55,
    })))['status'] == 'ok'


def test_rejects_surface_source_without_orientation_or_solution_parity():
    result = divergence_free_source_assembly_gate({
        "schema": "radia.divergence-free-source-assembly.v1",
        "scope": "affine_p1_constant_mu",
        "source_identity": "divergence_free",
        "boundary_orientation": "boundary_label",
        "manufactured_load_relative_error": 1e-3,
        "solved_field_relative_difference": 0.4,
        "allowed_solved_field_relative_difference": 1e-3,
        "volume_rhs_seconds": 8.0,
        "surface_rhs_seconds": 1.0,
    })
    assert result["status"] == "needs_attention"
    assert "boundary_orientation_is_element_derived" in result["issues"]
    assert "solved_field_matches" in result["issues"]


def test_accepts_hierarchical_p2_surface_source_evidence():
    result = divergence_free_source_assembly_gate({
        "schema": "radia.divergence-free-source-assembly.v1",
        "scope": "affine_p2_constant_mu",
        "response_order": 2,
        "source_identity": "divergence_free",
        "boundary_orientation": "opposite_vertex",
        "manufactured_load_relative_error": 8e-15,
        "solved_field_relative_difference": 8e-4,
        "allowed_solved_field_relative_difference": 1e-3,
        "volume_rhs_seconds": 390.0,
        "surface_rhs_seconds": 6.0,
    })
    assert result["status"] == "ok"
    assert result["response_order"] == 2


def test_rejects_response_order_scope_mismatch():
    summary = {
        "schema": "radia.divergence-free-source-assembly.v1",
        "scope": "affine_p1_constant_mu",
        "response_order": 2,
        "source_identity": "divergence_free",
        "boundary_orientation": "opposite_vertex",
        "manufactured_load_relative_error": 1e-15,
        "solved_field_relative_difference": 1e-4,
        "allowed_solved_field_relative_difference": 1e-3,
        "volume_rhs_seconds": 10.0,
        "surface_rhs_seconds": 1.0,
    }
    result = divergence_free_source_assembly_gate(summary)
    assert "scope_matches_response_order" in result["issues"]


def good_summary():
    exact = [
        0.00792907237221318,
        0.007902349246165159,
        0.007818835918418848,
        0.007670109587620525,
    ]
    return {
        "schema": "radia.finite-section-helmholtz-mixed-omega.v1",
        "units": {"length": "m", "field": "T", "current_density": "A/m^2"},
        "geometry": {
            "radial_bounds_m": [10.0, 13.0],
            "positive_axial_bounds_m": [3.0, 6.0],
            "mirror_pair": True,
            "current_density_A_per_m2": 10000.0,
        },
        "source": {
            "representation": "coil_builder_finite_cross_section",
            "closed_paths": True,
        },
        "open_boundary": {"method": "kelvin_transform", "radius_m": 20.0},
        "axis_z_m": [0.0, 1.0, 2.0, 3.0],
        "reference_bz_T": exact,
        "coil_builder_field_T": [
            [-2.42e-14, 0.0, 0.007929072379766044],
            [-2.40e-14, 0.0, 0.00790234925275819],
            [-2.25e-14, 0.0, 0.007818835921249668],
            [-1.96e-14, 0.0, 0.007670109583852553],
        ],
        "mixed_total_reduced": {
            "formulation": "mixed_total_reduced_omega",
            "field_T": [
                [2.5060801242412015e-06, -3.971699417892497e-06, 0.007931620193590618],
                [2.6535937266985635e-06, -4.469148949450773e-06, 0.007904718445317082],
                [2.5097592752474656e-06, -4.73770926732656e-06, 0.007820880092642974],
                [2.245642341972054e-06, -4.7301580172209686e-06, 0.007671789916063336],
            ],
            "refinement_levels": [
                {"mesh_maxh_m": 3.0, "order": 2, "cell_count": 8011, "ndof": 14226, "maximum_relative_error": 0.00614167},
                {"mesh_maxh_m": 2.0, "order": 3, "cell_count": 23274, "ndof": 121637, "maximum_relative_error": 0.000635068},
                {"mesh_maxh_m": 2.0, "order": 4, "cell_count": 19429, "ndof": 226944, "maximum_relative_error": 0.0003213265383183086},
            ],
        },
        "gate_tolerances": {
            "maximum_reference_relative_error": 1.0e-11,
            "maximum_coil_builder_relative_error": 1.0e-8,
            "maximum_mixed_relative_error": 4.0e-4,
            "maximum_coil_transverse_ratio": 1.0e-10,
            "maximum_mixed_transverse_ratio": 1.0e-3,
        },
    }


def test_axis_integral_matches_high_order_reference():
    value = finite_section_helmholtz_axis_bz(
        0.0,
        radial_bounds_m=[10.0, 13.0],
        positive_axial_bounds_m=[3.0, 6.0],
        current_density_A_per_m2=10000.0,
    )
    assert abs(value - 0.00792907237221318) / value < 1.0e-12


def test_accepts_coil_builder_and_converged_mixed_omega_result():
    result = finite_section_helmholtz_mixed_omega_gate(good_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_mixed_relative_error"] < 4.0e-4
    assert json.loads(mcp_gate(json.dumps(good_summary())))["status"] == "ok"


def test_core_mcp_catalog_and_runner_expose_gate():
    catalog = mcp._tool_manager._tools["radia_ngsolve_validation_catalog"].fn(
        query="finite_section_helmholtz"
    )
    assert [item["name"] for item in catalog["operations"]] == [
        "finite_section_helmholtz_mixed_omega_gate"
    ]
    result = asyncio.run(
        mcp._tool_manager._tools["radia_ngsolve_validation_run"].fn(
            name="finite_section_helmholtz_mixed_omega_gate",
            arguments={"summary_json": json.dumps(good_summary())},
        )
    )
    assert json.loads(result)["status"] == "ok"


def test_rejects_source_identity_and_false_refinement():
    bad = copy.deepcopy(good_summary())
    bad["source"]["representation"] = "filament_only"
    bad["mixed_total_reduced"]["refinement_levels"][-1]["maximum_relative_error"] = 0.01
    result = finite_section_helmholtz_mixed_omega_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_is_coil_builder_finite_section"] is False
    assert result["checks"]["refinement_error_decreases"] is False
    assert result["checks"]["final_refinement_matches_field"] is False


def test_rejects_axis_symmetry_leakage():
    bad = copy.deepcopy(good_summary())
    bad["mixed_total_reduced"]["field_T"][0][0] = 1.0e-3
    result = finite_section_helmholtz_mixed_omega_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["mixed_solution_preserves_axis_symmetry"] is False


def test_rejects_fractional_quadrature_order():
    with pytest.raises(ValueError, match="quadrature_order must be an integer"):
        finite_section_helmholtz_axis_bz(
            0.0,
            radial_bounds_m=[10.0, 13.0],
            positive_axial_bounds_m=[3.0, 6.0],
            current_density_A_per_m2=10000.0,
            quadrature_order=48.5,
        )


def test_rejects_fractional_refinement_metadata():
    bad = copy.deepcopy(good_summary())
    bad["mixed_total_reduced"]["refinement_levels"][0]["order"] = 2.5
    with pytest.raises(ValueError, match=r"refinement_levels\[0\]\.order must be an integer"):
        finite_section_helmholtz_mixed_omega_gate(bad)
