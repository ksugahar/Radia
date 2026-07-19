import asyncio
import json

from radia_mcp.matlab import (
    matlab_extension_contract,
    matlab_extension_path,
    matlab_official_server_config,
)


def test_generic_extension():
    contract = matlab_extension_contract()

    assert contract["ok"]
    assert contract["tool_count"] == 43
    assert contract["matlab_function_count"] == 86


def test_official_server_composition():
    config = matlab_official_server_config(
        "new_nodesktop",
        include_generic_extension=True,
    )

    assert len(config["extension_files"]) == 1
    assert config["matlab_setup_code"].startswith("addpath(")


def test_server_tools():
    from radia_mcp.matlab.server import mcp, matlab_extension_contract as tool

    tool_names = {item.name for item in asyncio.run(mcp.list_tools())}

    assert "matlab_extension_contract" in tool_names
    assert json.loads(tool())["tool_count"] == 43


def test_extension_schema_matches_matlab_function_contract():
    extension_path = matlab_extension_path()
    payload = json.loads(extension_path.read_text(encoding="utf-8"))
    tools = payload["tools"]
    signatures = payload["signatures"]
    tool_names = [tool["name"] for tool in tools]
    matlab_root = extension_path.parent.parent / "matlab" / "+radia_mcp_matlab"
    matlab_functions = {path.stem for path in matlab_root.glob("*.m")}

    assert len(tool_names) == len(set(tool_names))
    assert set(tool_names) == set(signatures)

    for tool in tools:
        name = tool["name"]
        schema = tool["inputSchema"]
        properties = schema["properties"]
        required = schema["required"]
        argument_order = signatures[name]["input"]["order"]
        function_name = signatures[name]["function"].rsplit(".", maxsplit=1)[-1]

        assert argument_order == required
        assert set(argument_order) == set(properties)
        assert function_name in matlab_functions

        seed_schema = properties["seed"]
        assert seed_schema["type"] == "integer"
        assert seed_schema["minimum"] == 0
