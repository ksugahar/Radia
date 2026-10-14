from __future__ import annotations

import copy
import json

import pytest

from radia_mcp.motor.pm_armature_reaction_gate import (
    pm_armature_reaction_hdiv_hex_gate,
)
from radia_mcp.motor.server import motor_pm_armature_reaction_hdiv_hex_gate


SPACE = {
    "family": "HDiv",
    "order": 1,
    "cell_family": "HEX",
    "project_lane": "BDM1",
    "strict_name": "tensor_product_hdiv_order1",
    "simplex_analogue": "BDM1",
}


def _fixture() -> dict[str, object]:
    identity = "a" * 64
    levels = []
    for count, incremental, correlation, absolute, knee in (
        (24, 0.115, 0.9993, 0.180, False),
        (192, 0.045, 0.9998, 0.110, True),
        (648, 0.025, 0.9999, 0.072, True),
    ):
        levels.append({
            "hex_element_count": count,
            "discrete_space": copy.deepcopy(SPACE),
            "delta_B_normalized_rms_difference": incremental,
            "delta_B_waveform_correlation": correlation,
            "zero_current_B_normalized_rms_difference": absolute,
            "loaded_B_normalized_rms_difference": absolute * 0.95,
            "zero_current_knee_classification_match": knee,
            "loaded_knee_classification_match": True,
            "coil": {
                "representation": "finite_section_gauss_filaments",
                "quadrature": [8, 8],
                "filament_count": 64,
            },
        })
    return {
        "schema": "radia.hdiv-hex-pm-armature-reaction-evidence.v1",
        "lane": "hdiv_mmm",
        "physics_identity_sha256": identity,
        "levels": levels,
        "reference": {"executed": True, "physics_identity_sha256": identity},
        "research_lab_retirement_ready": False,
        "product_or_market_retirement_ready": False,
    }


def test_accepts_incremental_response_but_keeps_absolute_demag_pending() -> None:
    result = pm_armature_reaction_hdiv_hex_gate(_fixture())
    assert result["status"] == "validated_partial"
    assert result["armature_reaction_increment_validated"] is True
    assert result["absolute_self_demagnetizing_field_validated"] is False
    assert result["research_lab_retirement_ready"] is False


def test_accepts_absolute_scope_only_when_its_independent_tolerance_is_met() -> None:
    artifact = _fixture()
    for level, absolute in zip(artifact["levels"], (0.08, 0.04, 0.02)):
        level["zero_current_B_normalized_rms_difference"] = absolute
        level["loaded_B_normalized_rms_difference"] = absolute * 0.95
    result = pm_armature_reaction_hdiv_hex_gate(artifact)
    assert result["status"] == "validated"
    assert result["absolute_self_demagnetizing_field_validated"] is True
    assert result["irreversible_demagnetization_state_validated"] is False


def test_rejects_non_hex_or_filament_reconstruction_drift() -> None:
    artifact = _fixture()
    artifact["levels"][1]["discrete_space"]["cell_family"] = "TET"
    artifact["levels"][2]["coil"]["filament_count"] = 1
    result = pm_armature_reaction_hdiv_hex_gate(artifact)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_levels_use_bdm1_project_lane_on_hex"] is False
    assert result["checks"][
        "finite_section_coil_is_reconstructed_by_tensor_quadrature"
    ] is False


def test_mcp_tool_preserves_partial_status_and_rejects_overclaim() -> None:
    payload = _fixture()
    response = json.loads(motor_pm_armature_reaction_hdiv_hex_gate(json.dumps(payload)))
    assert response["status"] == "validated_partial"

    payload["research_lab_retirement_ready"] = True
    response = json.loads(motor_pm_armature_reaction_hdiv_hex_gate(json.dumps(payload)))
    assert response["status"] == "needs_attention"
    assert "scope_does_not_overclaim_retirement" in response["issues"]


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("delta_B_normalized_rms_difference", -0.01, "nonnegative"),
        ("delta_B_waveform_correlation", 1.01, r"\[-1, 1\]"),
    ],
)
def test_rejects_nonphysical_error_metrics(field, value, message) -> None:
    artifact = _fixture()
    artifact["levels"][-1][field] = value

    with pytest.raises(ValueError, match=message):
        pm_armature_reaction_hdiv_hex_gate(artifact)
