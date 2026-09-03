"""Contracts for the thin-server, coarse-workflow MCP surface."""

from __future__ import annotations

import asyncio
import ast
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from mcp.server.fastmcp import FastMCP

from radia_mcp.common.lazy_call import lazy_callable
from radia_mcp.common.status import register_status_tool
from radia_mcp.common.tool_group import CoarseToolRegistry, selected_tool_profile
from radia_mcp.meta.catalog import CATALOG


def _tool_functions(mcp: FastMCP) -> dict[str, object]:
    return {
        name: tool.fn
        for name, tool in mcp._tool_manager._tools.items()
    }


def test_common_package_does_not_eagerly_import_optional_subsystems():
    package_src = Path(__file__).parents[1] / "src"
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_src)
    code = (
        "import json,sys; import radia_mcp.common; "
        "print(json.dumps(sorted(n for n in sys.modules "
        "if n.startswith('radia_mcp.common.'))))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
    )
    assert "chroma_retriever" not in result.stdout
    assert "async_runner" not in result.stdout


def test_small_server_status_import_does_not_preload_chroma():
    package_src = Path(__file__).parents[1] / "src"
    env = dict(os.environ)
    env["PYTHONPATH"] = str(package_src)
    code = (
        "import json,sys; import radia_mcp.accelerator.server; "
        "print(json.dumps(sorted(n for n in sys.modules "
        "if n.startswith('radia_mcp.common.'))))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
    )
    loaded = set(json.loads(result.stdout.splitlines()[-1]))
    assert "radia_mcp.common.status" in loaded
    assert "radia_mcp.common.chroma_retriever" not in loaded
    assert "radia_mcp.common.async_runner" not in loaded


def test_selected_tool_profile_defaults_to_core(monkeypatch):
    monkeypatch.delenv("RADIA_MCP_TOOL_PROFILE", raising=False)
    assert selected_tool_profile([]) == "core"
    assert selected_tool_profile(["--tool-profile", "full"]) == "full"
    assert selected_tool_profile(["--tool-profile=full"]) == "full"


def test_selected_tool_profile_rejects_unknown_value(monkeypatch):
    monkeypatch.setenv("RADIA_MCP_TOOL_PROFILE", "wide")
    with pytest.raises(ValueError, match="invalid MCP tool profile"):
        selected_tool_profile([])


def test_core_profile_exposes_catalog_and_runner_only():
    mcp = FastMCP("test")
    registry = CoarseToolRegistry(
        mcp,
        namespace="demo",
        category="validation",
        profile="core",
        min_group_size=1,
    )

    @registry.tool()
    def demo_energy_gate(value: float, limit: float = 1.0) -> dict:
        """Check an energy limit."""
        return {"passed": value <= limit}

    registry.install()
    functions = _tool_functions(mcp)
    assert set(functions) == {
        "demo_validation_catalog",
        "demo_validation_run",
    }
    catalog_tool = mcp._tool_manager._tools["demo_validation_catalog"]
    runner_tool = mcp._tool_manager._tools["demo_validation_run"]
    assert catalog_tool.annotations.readOnlyHint is True
    assert runner_tool.annotations.readOnlyHint is False
    assert runner_tool.annotations.destructiveHint is False
    assert runner_tool.annotations.idempotentHint is False
    catalog = functions["demo_validation_catalog"](query="energy")
    assert catalog["matched"] == 1
    assert catalog["operations"][0]["name"] == "demo_energy_gate"
    result = asyncio.run(
        functions["demo_validation_run"](
            name="demo_energy_gate",
            arguments={"value": 2.0},
        )
    )
    assert result == {"passed": False}


def test_full_profile_keeps_individual_tools_for_compatibility():
    mcp = FastMCP("test")
    registry = CoarseToolRegistry(
        mcp, namespace="demo", profile="full", min_group_size=1
    )

    @registry.tool()
    def demo_identity_gate(value: int) -> int:
        return value

    registry.install()
    assert set(_tool_functions(mcp)) == {
        "demo_identity_gate",
        "demo_validation_catalog",
        "demo_validation_run",
    }


def test_lazy_callable_resolves_only_on_call_and_tracks_reloaded_attribute(
    monkeypatch,
):
    calls: list[str] = []

    class Target:
        operation = staticmethod(lambda value: value + 1)

    target = Target()

    def fake_import(name: str):
        calls.append(name)
        return target

    monkeypatch.setattr("radia_mcp.common.lazy_call.import_module", fake_import)
    proxy = lazy_callable(".checks", "operation", "radia_mcp.demo")
    assert calls == []
    assert proxy(2) == 3
    target.operation = lambda value: value + 10
    assert proxy(2) == 12
    assert calls == ["radia_mcp.demo.checks", "radia_mcp.demo.checks"]


def test_status_reports_grouped_surface():
    mcp = FastMCP("test")
    registry = CoarseToolRegistry(
        mcp, namespace="demo", profile="core", min_group_size=1
    )

    @registry.tool()
    def demo_gate() -> bool:
        return True

    registry.install()
    register_status_tool(
        mcp,
        server_name="mcp-server-demo",
        description="Demo server",
        subpackage="radia_mcp.demo",
    )
    status = _tool_functions(mcp)["demo_status"]()
    assert status["tool_profile"] == "core"
    assert status["n_grouped_operations"] == 1
    assert status["tool_groups"][0]["namespace"] == "demo"
    assert status["tool_groups"][0]["individual_tools_exposed"] == 0


def test_small_family_remains_direct_when_pair_would_be_larger():
    mcp = FastMCP("test")
    registry = CoarseToolRegistry(mcp, namespace="demo", profile="core")

    @registry.tool()
    def demo_single_gate() -> bool:
        return True

    registry.install()
    assert set(_tool_functions(mcp)) == {"demo_single_gate"}
    assert mcp._radia_tool_groups[0]["mode"] == "direct-small-group"
    assert mcp._radia_tool_groups[0]["grouped_operations"] == 0


def test_primary_catalog_tools_remain_directly_registered():
    package_root = Path(__file__).parents[1] / "src" / "radia_mcp"
    primary_tools = {
        name
        for server in CATALOG.values()
        for name in server.get("primary_tools", ())
    }
    grouped: set[str] = set()
    for server_path in package_root.glob("*/server.py"):
        tree = ast.parse(server_path.read_text(encoding="utf-8"))
        for node in tree.body:
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for decorator in node.decorator_list:
                target = decorator.func if isinstance(decorator, ast.Call) else decorator
                if (
                    isinstance(target, ast.Attribute)
                    and isinstance(target.value, ast.Name)
                    and target.value.id == "_validation"
                    and target.attr == "tool"
                ):
                    grouped.add(node.name)
    assert not grouped.intersection(primary_tools)


def test_each_grouped_server_installs_its_registry_once():
    package_root = Path(__file__).parents[1] / "src" / "radia_mcp"
    for server_path in package_root.glob("*/server.py"):
        source = server_path.read_text(encoding="utf-8")
        grouped_count = source.count("@_validation.tool(")
        if grouped_count:
            assert source.count("_validation.install()") == 1, server_path
            assert "CoarseToolRegistry" in source, server_path
