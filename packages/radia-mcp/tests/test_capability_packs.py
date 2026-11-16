"""Protect discovery and protocol behavior while reducing MCP processes."""

import asyncio
import importlib
import os
import sys
from pathlib import Path

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.types import Implementation

from radia_mcp.capability_packs import PACKS, modules_for
from radia_mcp.capability_packs.server import CapabilityServer


@pytest.mark.parametrize("pack", PACKS)
def test_pack_preserves_domain_schemas_prompts_resources(pack):
    async def check():
        app = await CapabilityServer(pack).initialize()
        expected = {}
        prompts = {}
        resources = {}
        for module in modules_for(pack):
            source = importlib.import_module(f"radia_mcp.{module}.server").mcp
            expected.update({t.name: t.model_dump() for t in await source.list_tools()})
            prompts.update({p.name: p.model_dump() for p in await source.list_prompts()})
            resources.update({str(r.uri): r.model_dump() for r in await source.list_resources()})
        actual = {t.name: t.model_dump() for t in await app.list_tools()}
        actual.pop("capability_pack_status")
        assert actual == expected
        status = app.status()
        assert status["runtime_contract"]["complete"]
        assert status["runtime_contract"]["n_tools"] == len(expected) + 1
        assert status["runtime_provenance"]["module_sha256_at_registration"]
        assert not status["runtime_provenance"]["source_changed_since_registration"]
        assert {p.name: p.model_dump() for p in await app.list_prompts()} == prompts
        assert {str(r.uri): r.model_dump() for r in await app.list_resources()} == resources
    asyncio.run(check())


def test_unknown_profile_fails_before_import():
    with pytest.raises(ValueError, match="Unknown profile"):
        CapabilityServer("radia-mcp-optimize", "typo")


def test_duplicate_tool_is_rejected(monkeypatch):
    monkeypatch.setitem(PACKS, "duplicate", {
        "description": "test", "profiles": {"test": ("bayesian_opt", "bayesian_opt")},
    })
    with pytest.raises(ValueError, match="Duplicate capability"):
        asyncio.run(CapabilityServer("duplicate", "test").initialize())


@pytest.mark.parametrize("pack,profile", [("radia-mcp-analysis", "field")])
def test_real_stdio_discovery_call_and_errors(pack, profile, tmp_path):
    async def check():
        env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1] / "src"))
        params = StdioServerParameters(command=sys.executable, args=[
            "-X", "utf8", "-m", f"{PACKS[pack]['subpackage']}.server",
            "--profile", profile,
        ], env=env)
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write, client_info=Implementation(
                    name="pack-contract", version="1.0")) as client:
                await client.initialize()
                listing = await client.list_tools()
                assert listing.tools
                status = await client.call_tool("capability_pack_status", {})
                assert not status.isError
                assert status.structuredContent["modules"] == list(modules_for(pack, profile))
                again = await client.call_tool("capability_pack_status", {})
                connection = again.structuredContent["client_connection"]
                assert connection["connection_count"] == 1
                assert connection["latest"]["name"] == "pack-contract"
                assert connection["latest"]["version"] == "1.0"
                bad = await client.call_tool("no_such_tool", {})
                assert bad.isError
                resources = await client.list_resources()
                assert resources.resources
                content = await client.read_resource(resources.resources[0].uri)
                assert content.contents
                prompts = await client.list_prompts()
                assert prompts.prompts
                prompt = prompts.prompts[0]
                result = await client.get_prompt(prompt.name, {
                    arg.name: "test" for arg in (prompt.arguments or []) if arg.required
                })
                assert result.messages
    asyncio.run(check())


def test_pack_tool_list_follows_a_domain_registry_change():
    async def check():
        app = await CapabilityServer("radia-mcp-optimize", "search").initialize()
        module = modules_for("radia-mcp-optimize", "search")[0]
        source = importlib.import_module(f"radia_mcp.{module}.server").mcp

        def pack_probe_tool() -> str:
            return "probe"

        source.add_tool(pack_probe_tool)
        try:
            assert "pack_probe_tool" in {t.name for t in await app.list_tools()}
        finally:
            source.remove_tool("pack_probe_tool")
        assert "pack_probe_tool" not in {t.name for t in await app.list_tools()}
        with pytest.raises(ValueError, match="Unknown prompt"):
            await app.get_prompt("no_such_prompt", {})
    asyncio.run(check())
