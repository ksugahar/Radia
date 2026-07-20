"""Regression lock for the public radia-mcp server naming contract."""

from pathlib import Path

import tomllib


def test_all_mcp_server_scripts_are_cataloged_with_canonical_names():
    from radia_mcp.meta import catalog

    package_root = Path(__file__).resolve().parents[1]
    project = tomllib.loads((package_root / "pyproject.toml").read_text(encoding="utf-8"))
    scripts = {name for name in project["project"]["scripts"] if name.startswith("mcp-server-")}
    catalog_commands = {entry["command"] for entry in catalog.list_all()}

    assert scripts == catalog_commands
    for entry in catalog.list_all():
        server_id = entry["command"].removeprefix("mcp-server-")
        assert entry["server_id"] == server_id
        assert entry["product_name"] == f"radia-mcp.{server_id}"
        assert entry["python_name"] == entry["subpackage"]
        assert catalog.get(entry["command"])["name"] == entry["name"]
        assert catalog.get(server_id)["name"] == entry["name"]
        assert catalog.get(server_id.replace("-", "_"))["name"] == entry["name"]


def test_naming_conventions_tool_exposes_contract():
    from radia_mcp.meta.server import radia_mcp_naming_conventions

    result = radia_mcp_naming_conventions()
    assert result["product"] == "radia-mcp.<server-id>"
    assert result["command"] == "mcp-server-<server-id>"
    assert result["python"] == "radia_mcp.<python_module>"
    assert result["servers"]
