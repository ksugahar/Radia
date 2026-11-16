"""Regression lock for the public radia-mcp server naming contract."""

import re
import tomllib
from pathlib import Path


def test_single_command_serves_every_hyphenated_catalog_key():
    from radia_mcp.meta import catalog

    package_root = Path(__file__).resolve().parents[1]
    project = tomllib.loads((package_root / "pyproject.toml").read_text(encoding="utf-8"))
    servers = [name for name in project["project"]["scripts"] if "mcp-server" in name]
    assert servers == [], servers
    assert project["project"]["scripts"]["radia-mcp"] == "radia_mcp.maintenance:main"
    for entry in catalog.list_all():
        assert re.fullmatch(r"[a-z0-9]+(-[a-z0-9]+)*", entry["name"]), entry["name"]
        assert entry["python_name"] == entry["subpackage"]
        if entry.get("kind") == "domain":
            assert "command" not in entry, entry["name"]
        else:
            assert entry["command"] == f"radia-mcp serve {entry['name']}"
            assert catalog.get(entry["name"])["name"] == entry["name"]


def test_retired_names_do_not_resolve():
    from radia_mcp.meta import catalog

    for old in ("mcp-server-radia-meta", "radia-" + "analysis", "radia-ltspice", "panel-review",
                "mcp-ltspice", "radia-" + "matlab", "radia_mcp_meta",
                "magnetic_materials", "gmsh_mcp"):
        assert catalog.get(old) is None, old


def test_naming_conventions_tool_exposes_contract():
    from radia_mcp.meta.server import radia_mcp_naming_conventions

    result = radia_mcp_naming_conventions()
    assert result["command"] == "radia-mcp serve <catalog-key>"
    assert result["python"] == "radia_mcp.<python_module>"
    assert result["servers"]
    assert "cubit_show" in result["compatibility"]
    assert "No aliases" in result["compatibility"]
