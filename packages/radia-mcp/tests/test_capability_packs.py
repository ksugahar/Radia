"""Protect discovery and protocol behavior while reducing MCP processes."""

import asyncio
import importlib
import os
import sys
from pathlib import Path

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

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
        assert {p.name: p.model_dump() for p in await app.list_prompts()} == prompts
        assert {str(r.uri): r.model_dump() for r in await app.list_resources()} == resources
    asyncio.run(check())


def test_unknown_profile_fails_before_import():
    with pytest.raises(ValueError, match="Unknown profile"):
        CapabilityServer("radia-design", "typo")


@pytest.mark.parametrize("language", ["ja", "en"])
def test_writing_route_separates_scoring_objectives(language):
    from radia_mcp.paper_writing.review_route import paper_writing_review_route
    paper = paper_writing_review_route("paper", language)
    grant = paper_writing_review_route("grant_proposal", language)
    assert set(paper["tools"]).isdisjoint(grant["tools"])
    assert paper["aggregate_score"] is grant["aggregate_score"] is None
    assert not paper["venue_requirements_verified"]
    if language == "en":
        assert "No validated English grant score" in grant["language_scoring"]


def test_duplicate_tool_is_rejected(monkeypatch):
    monkeypatch.setitem(PACKS, "duplicate", {
        "description": "test", "profiles": {"test": ("bayesian_opt", "bayesian_opt")},
    })
    with pytest.raises(ValueError, match="Duplicate capability"):
        asyncio.run(CapabilityServer("duplicate", "test").initialize())


def test_grant_and_poster_are_available_in_paper_writing():
    from radia_mcp.paper_writing.server import mcp
    from radia_mcp.grant_writing.server import mcp as grant
    from radia_mcp.poster.server import mcp as poster

    async def check():
        tools = {t.name: t for t in await mcp.list_tools()}
        for source, prefix in ((grant, "grant_writing_"), (poster, "poster_")):
            for tool in await source.list_tools():
                if tool.name.startswith(prefix) and not tool.name.endswith(("_status", "_reload_code")):
                    assert tools[tool.name].inputSchema == tool.inputSchema
                    assert tools[tool.name].outputSchema == tool.outputSchema
                    if prefix == "grant_writing_":
                        assert tools[tool.name].annotations.readOnlyHint
        assert "never average scores" in mcp.instructions
    asyncio.run(check())


@pytest.mark.parametrize("pack,profile", [("document-ops", "convert"), ("radia-analysis", "field")])
def test_real_stdio_discovery_call_and_errors(pack, profile, tmp_path):
    async def check():
        env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1] / "src"))
        params = StdioServerParameters(command=sys.executable, args=[
            "-X", "utf8", "-m", f"radia_mcp.{pack.replace('-', '_')}.server",
            "--profile", profile,
        ], env=env)
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as client:
                await client.initialize()
                listing = await client.list_tools()
                assert listing.tools
                status = await client.call_tool("capability_pack_status", {})
                assert not status.isError
                assert status.structuredContent["modules"] == list(modules_for(pack, profile))
                bad = await client.call_tool("no_such_tool", {})
                assert bad.isError
                if pack == "document-ops":
                    invalid = await client.call_tool("md2html_convert", {})
                    assert invalid.isError
                    if importlib.util.find_spec("markdown") is None:
                        return  # Conversion is optional; protocol/error checks still ran.
                    source = tmp_path / "input.md"
                    target = tmp_path / "output.html"
                    source.write_text("# Capability test\n", encoding="utf-8")
                    converted = await client.call_tool("md2html_convert", {
                        "md_file": str(source), "output_file": str(target),
                    })
                    assert not converted.isError
                    assert "Capability test" in target.read_text(encoding="utf-8-sig")
                else:
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
