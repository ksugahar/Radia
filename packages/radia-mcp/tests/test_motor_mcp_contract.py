from radia_mcp.common.mcp_contract import SCHEMA
from radia_mcp.motor.server import mcp


def test_motor_server_has_mathworks_style_runtime_contract():
    low_level = mcp._mcp_server
    assert low_level.version == "1.4.19"
    assert low_level.instructions

    tools = mcp._tool_manager._tools
    assert tools
    for name, tool in tools.items():
        assert tool.title, name
        assert tool.annotations is not None, name
        assert tool.annotations.readOnlyHint is not None, name
        assert tool.annotations.destructiveHint is not None, name
        assert tool.annotations.idempotentHint is not None, name
        assert tool.annotations.openWorldHint is not None, name
        assert tool.meta["caeai.contract"] == SCHEMA


def test_circuit_age_plan_is_read_only_and_closed_world():
    tool = mcp._tool_manager._tools["motor_circuit_age_application_plan"]
    assert tool.annotations.readOnlyHint is True
    assert tool.annotations.destructiveHint is False
    assert tool.annotations.idempotentHint is True
    assert tool.annotations.openWorldHint is False

    analysis = mcp._tool_manager._tools["motor_circuit_field_analysis"]
    assert analysis.annotations.readOnlyHint is True
    assert analysis.annotations.destructiveHint is False
    assert analysis.annotations.idempotentHint is True
    assert analysis.annotations.openWorldHint is False

    vol_analysis = mcp._tool_manager._tools["motor_vol2d_circuit_analysis"]
    assert vol_analysis.annotations.readOnlyHint is True
    assert vol_analysis.annotations.destructiveHint is False
    assert vol_analysis.annotations.idempotentHint is True
    assert vol_analysis.annotations.openWorldHint is False

    dynamic = mcp._tool_manager._tools["motor_vol2d_dynamic_analysis"]
    assert dynamic.annotations.readOnlyHint is True
    assert dynamic.annotations.destructiveHint is False
    assert dynamic.annotations.idempotentHint is True
    assert dynamic.annotations.openWorldHint is False
