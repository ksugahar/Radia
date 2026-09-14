"""Regression coverage for public paths identified by the September review."""

import pytest
from radia_mcp.presentation.cross_lint import presentation_lint_bedrock
from radia_mcp.radia_ngsolve import slot_gates


def test_presentation_bedrock_uses_shared_hedge_scanner():
    result = presentation_lint_bedrock("改善すると考えられる。効果があると思われる。")
    hedges = [issue for issue in result["issues"] if issue["rule"] == "hedging_expressions"]
    assert len(hedges) == 1
    assert hedges[0]["count"] == 2


@pytest.mark.parametrize("coordinates", ["1, 2, 3", "1;2;3", "1 2 3"])
@pytest.mark.parametrize("gate_name,packet,keyword,sequence", [
    ("cst_touchstone_solver_ready_manifest_gate", [{}], "expected_port_face_centers_xyz_m", True),
    ("jmag_force_table_metadata_gate", {"columns": ["Fx"]}, "expected_target_region_centroid_xyz_m", False),
    ("pm_demag_margin_screening_package_gate", [{}], "expected_field_probe_point_xyz_m", False),
])
def test_public_gate_string_coordinates_match_numeric_coordinates(
    coordinates, gate_name, packet, keyword, sequence,
):
    gate = getattr(slot_gates, gate_name)
    numeric = [[1, 2, 3]] if sequence else [1, 2, 3]
    encoded = [coordinates] if sequence else coordinates
    # Incomplete packets compare normalization, not numerical acceptance.
    expected = gate(packet, **{keyword: numeric})
    assert gate(packet, **{keyword: encoded}) == expected


@pytest.mark.parametrize("error", [SystemExit(3), KeyboardInterrupt(), RuntimeError("broken")])
def test_background_failure_reaches_terminal_state_and_can_restart(error):
    from radia_mcp.common.async_runner import AsyncRunner

    def fail():
        raise error

    runner = AsyncRunner()
    assert runner.start(fail)
    result = runner.wait(timeout=5)
    assert result["status"] == "failed"
    assert result["end_time"] is not None
    assert not runner.is_running
    assert runner.start(lambda: 42)
    assert runner.wait(timeout=5)["result"] == 42
