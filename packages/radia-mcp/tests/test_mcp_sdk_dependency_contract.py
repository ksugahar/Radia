"""Keep package metadata and minimal CI on the supported FastMCP SDK line."""

import asyncio
import importlib.metadata
import json
import os
from pathlib import Path
import re

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    from importlib import import_module

    tomllib = import_module("tomli")


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SDK_REQUIREMENT = "mcp>=1.20.0,<2"


def test_mcp_sdk_dependency_declares_supported_floor():
    metadata = tomllib.loads(
        (PACKAGE_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )
    project = metadata["project"]
    assert SDK_REQUIREMENT in project["dependencies"]
    assert metadata["build-system"]["requires"] == ["setuptools>=77.0"]
    assert project["license"] == "BSD-3-Clause"
    assert project["license-files"] == ["LICENSE"]
    from packaging.requirements import Requirement
    names = {Requirement(item).name for item in project['dependencies']}
    assert not {'cae-mcp-core', 'cubit-mesh-export', 'radia'}.intersection(names)


def test_public_description_and_optional_cubit_ownership():
    from packaging.requirements import Requirement

    project = tomllib.loads((PACKAGE_ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    assert len(project["description"]) <= 512
    assert "distributed separately in cubit-mesh-export" in project["description"]
    assert "mcp-server-cubit" not in project["scripts"]
    requirement, = map(Requirement, project["optional-dependencies"]["cubit"])
    assert requirement.name == "cubit-mesh-export"
    assert requirement.marker.evaluate({"python_version": "3.12", "sys_platform": "win32"})
    assert not requirement.marker.evaluate({"python_version": "3.12", "sys_platform": "linux"})
    readme = (PACKAGE_ROOT / "README.md").read_text(encoding="utf-8")
    assert "core (Cubit" not in readme
    assert "| **Cubit** | `mcp-server-cubit`" not in readme
    discovery_surfaces = [
        PACKAGE_ROOT / "README.md",
        PACKAGE_ROOT / "REGISTRY_SUBMISSION.md",
        PACKAGE_ROOT / "src" / "radia_mcp" / "common" / "status.py",
        *(PACKAGE_ROOT / "src" / "radia_mcp" / "meta").glob("*.py"),
    ]
    for path in discovery_surfaces:
        if path.exists():
            text = path.read_text(encoding="utf-8")
            assert re.search(r"\b\d+[- ]servers?\b", text, re.IGNORECASE) is None, path
    from radia_mcp.meta.catalog import CATALOG
    from radia_mcp.meta.server import radia_mcp_overview

    overview = radia_mcp_overview()
    assert overview["n_servers"] == len(CATALOG)
    assert len(overview["servers"]) == len(CATALOG)
    registry_path = PACKAGE_ROOT / "REGISTRY_SUBMISSION.md"
    if registry_path.exists():
        registry = registry_path.read_text(encoding="utf-8")
        assert "Cubit MCP is\ndistributed separately" in registry
        assert "Entry points (3 MCP servers shipped in one wheel)" not in registry
        assert "mcp-server-cubit         # Coreform Cubit" not in registry
    assert not (PACKAGE_ROOT / "docs" / "SHARED_LIB_DESIGN.md").exists()


def test_retired_cae_mcp_core_cannot_return_as_shared_policy() -> None:
    readme = (PACKAGE_ROOT / "README.md").read_text(encoding="utf-8")
    assert "`cae-mcp-core` foundation is retired" in readme
    assert "not a dependency or\nshared runtime" in readme



def test_supported_sdk_registers_metadata_lists_schema_and_calls_tool():
    from mcp.server.fastmcp import FastMCP
    from radia_mcp.common.server_hardening import ANN_READONLY

    expected = os.environ.get("RADIA_MCP_EXPECTED_SDK")
    if expected:
        assert importlib.metadata.version("mcp") == expected
    server = FastMCP("sdk-contract")

    def increment(value: int, step: int = 1) -> dict:
        """Increment a number without external state."""
        return {"value": value + step}

    # Same metadata arguments used when Radia refreshes registered tools.
    server.add_tool(increment, annotations=ANN_READONLY, icons=None,
                    meta={"radia": {"contract": "sdk-floor"}})

    async def exercise():
        tools = await server.list_tools()
        tool = next(item for item in tools if item.name == "increment")
        assert tool.meta == {"radia": {"contract": "sdk-floor"}}
        assert tool.annotations.readOnlyHint is True
        assert tool.inputSchema["required"] == ["value"]
        assert tool.inputSchema["properties"]["step"]["default"] == 1
        result = await server.call_tool("increment", {"value": "4"})
        # FastMCP may return content alone or (content, structured output).
        content = result[0] if isinstance(result, tuple) else result
        assert json.loads(content[0].text) == {"value": 5}

    asyncio.run(exercise())
