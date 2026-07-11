import json
import math

from radia_mcp.radia_ngsolve.field_profile_gate import symmetric_axial_field_profile_gate
from radia_mcp.radia_ngsolve.server import symmetric_axial_field_profile_gate as mcp_profile_gate


POSITIONS = [0.04472135955 * index for index in range(-10, 11)]
AXIAL_FIELD = [
    1.188564073157e-6, 1.209449943741e-6, 1.225699725776e-6,
    1.237706903688e-6, 1.246035332249e-6, 1.251364220809e-6,
    1.254423025513e-6, 1.255923332067e-6, 1.256494439400e-6,
    1.256628188252e-6, 1.256637069998e-6, 1.256628188252e-6,
    1.256494439400e-6, 1.255923332067e-6, 1.254423025513e-6,
    1.251364220809e-6, 1.246035332249e-6, 1.237706903688e-6,
    1.225699725776e-6, 1.209449943741e-6, 1.188564073157e-6,
]


def test_symmetric_axial_profile_accepts_full_profile_and_analytic_center():
    result = symmetric_axial_field_profile_gate(
        POSITIONS,
        AXIAL_FIELD,
        expected_center_field=4.0e-7 * math.pi,
        transverse_field_1=[0.0] * 21,
        transverse_field_2=[0.0] * 21,
    )
    assert result["status"] == "ok"
    assert result["metrics"]["center_relative_error"] < 1.0e-8
    assert result["metrics"]["field_symmetry_relative"] == 0.0


def test_symmetric_axial_profile_rejects_center_only_match_with_bad_profile():
    bad = AXIAL_FIELD.copy()
    bad[3] *= 0.8
    result = symmetric_axial_field_profile_gate(
        POSITIONS,
        bad,
        expected_center_field=AXIAL_FIELD[10],
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["axial_field_is_symmetric"] is False
    assert result["checks"]["field_increases_toward_center"] is False


def test_symmetric_axial_profile_mcp_tool_rejects_transverse_component():
    transverse = [0.0] * 21
    transverse[10] = 1.0e-7
    result = json.loads(mcp_profile_gate(
        POSITIONS,
        AXIAL_FIELD,
        AXIAL_FIELD[10],
        transverse_field_1=transverse,
    ))
    assert result["status"] == "needs_attention"
    assert result["checks"]["transverse_field_is_negligible"] is False
