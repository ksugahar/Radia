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
