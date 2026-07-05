import pytest

from radia_mcp.radia_ngsolve.axisymmetric_3d_validation import (
    axisymmetric_to_3d_force_gate,
    axisymmetric_to_3d_validation_plan,
)
from radia_mcp.radia_ngsolve.knowledge.force_validation import get_force_validation_documentation


def test_axisymmetric_to_3d_force_gate_accepts_full_revolution_force_vector():
    gate = axisymmetric_to_3d_force_gate(
        -12.0,
        [1.0e-5, -2.0e-5, -12.01],
        case_id="axi_full_3d_fixture",
        axial_axis="z",
        axial_rtol=0.002,
        transverse_rtol=1.0e-4,
    )

    assert gate["status"] == "ok"
    assert gate["axisymmetric_reference"]["quantity_basis"] == "full_3d_revolution_from_2pi_r_weight"
    assert gate["three_d_result"]["result_basis"] == "full_revolution"
    assert gate["checks"]["transverse_components_cancel"] is True
    assert gate["errors"]["axial_rel_error"] < 0.002


def test_axisymmetric_to_3d_force_gate_scales_symmetry_sector_without_overclaiming_transverse_cancellation():
    gate = axisymmetric_to_3d_force_gate(
        40.0,
        [0.25, 0.5, 10.0],
        case_id="axi_sector_fixture",
        result_basis="symmetry_sector",
        sector_angle_deg=90.0,
        axial_rtol=1.0e-12,
    )

    assert gate["status"] == "ok"
    assert gate["three_d_result"]["scale_to_full_revolution"] == pytest.approx(4.0)
    assert gate["three_d_result"]["axial_component_N"] == pytest.approx(40.0)
    assert gate["checks"]["transverse_components_cancel"] == "not_checked_for_symmetry_sector"


def test_axisymmetric_to_3d_force_gate_rejects_wrong_transverse_force_for_full_revolution():
    gate = axisymmetric_to_3d_force_gate(
        10.0,
        [1.0, 0.0, 10.0],
        case_id="bad_transverse_fixture",
        transverse_rtol=0.01,
    )

    assert gate["status"] == "needs_attention"
    assert gate["checks"]["axial_component_matches_axisymmetric_reference"] is True
    assert gate["checks"]["transverse_components_cancel"] is False


def test_axisymmetric_to_3d_validation_plan_names_validation_test_artifact_contract():
    plan = axisymmetric_to_3d_validation_plan("coaxial_loop_force")

    assert plan["policy"] == "axisymmetric_to_3d_validation_plan"
    assert plan["reference_quantity"]["recommended_extractor"].endswith("eggshell_force_axi")
    assert plan["three_d_quantity"]["recommended_extractor"].endswith("eggshell_force")
    assert any("validation_test/force_validation" in item for item in plan["required_artifacts"])
    assert "transverse force components cancel" in " ".join(plan["required_checks"])


def test_force_validation_method_map_exposes_axisymmetric_to_3d_gate():
    doc = get_force_validation_documentation("method_map")

    assert "Axisymmetric reference -> 3D force-vector validation" in doc
    assert "axisymmetric_to_3d_force_gate" in doc
    assert "full-revolution 3D axial force" in doc
    assert "validation_axisymmetric_to_3d_vol_force.py" in doc
