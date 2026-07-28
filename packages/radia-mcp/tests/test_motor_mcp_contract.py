from radia_mcp.common.mcp_contract import SCHEMA
import pytest

from radia_mcp.motor.server import _decode_owned_worker_json, mcp


def test_owned_worker_decoder_tolerates_native_diagnostics_before_json():
    payload = _decode_owned_worker_json(
        b"native diagnostic line\nsecond line\n{\"status\":\"solved\",\"value\":3}\n"
    )
    assert payload == {"status": "solved", "value": 3}


@pytest.mark.parametrize("payload", [b"", b"native diagnostic only\n"])
def test_owned_worker_decoder_rejects_missing_json(payload):
    with pytest.raises(ValueError, match="JSON"):
        _decode_owned_worker_json(payload)


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

    force = mcp._tool_manager._tools["motor_vol2d_force_analysis"]
    assert force.annotations.readOnlyHint is True
    assert force.annotations.destructiveHint is False
    assert force.annotations.idempotentHint is True
    assert force.annotations.openWorldHint is False
