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


# Literal imports keep each domain server statically resolvable for the
# impact selection; the coverage test fails when a pack adds a domain.
_SERVER_IMPORTS = {
    "accelerator": lambda: importlib.import_module("radia_mcp.accelerator.server"),
    "acoustic_fembem": lambda: importlib.import_module("radia_mcp.acoustic_fembem.server"),
    "bayesian_opt": lambda: importlib.import_module("radia_mcp.bayesian_opt.server"),
    "bem": lambda: importlib.import_module("radia_mcp.bem.server"),
    "blender": lambda: importlib.import_module("radia_mcp.blender.server"),
    "data_assimilation": lambda: importlib.import_module("radia_mcp.data_assimilation.server"),
    "differential_forms": lambda: importlib.import_module("radia_mcp.differential_forms.server"),
    "electromagnet": lambda: importlib.import_module("radia_mcp.electromagnet.server"),
    "evolutionary": lambda: importlib.import_module("radia_mcp.evolutionary.server"),
    "fem": lambda: importlib.import_module("radia_mcp.fem.server"),
    "force": lambda: importlib.import_module("radia_mcp.force.server"),
    "fusion_reactor": lambda: importlib.import_module("radia_mcp.fusion_reactor.server"),
    "gnn": lambda: importlib.import_module("radia_mcp.gnn.server"),
    "ih": lambda: importlib.import_module("radia_mcp.ih.server"),
    "litz_transmission": lambda: importlib.import_module("radia_mcp.litz_transmission.server"),
    "maglev": lambda: importlib.import_module("radia_mcp.maglev.server"),
    "magnetic_materials": lambda: importlib.import_module("radia_mcp.magnetic_materials.server"),
    "matrix_solvers": lambda: importlib.import_module("radia_mcp.matrix_solvers.server"),
    "metamaterial": lambda: importlib.import_module("radia_mcp.metamaterial.server"),
    "mor": lambda: importlib.import_module("radia_mcp.mor.server"),
    "motor": lambda: importlib.import_module("radia_mcp.motor.server"),
    "optimization": lambda: importlib.import_module("radia_mcp.optimization.server"),
    "pcb": lambda: importlib.import_module("radia_mcp.pcb.server"),
    "peec": lambda: importlib.import_module("radia_mcp.peec.server"),
    "pinn": lambda: importlib.import_module("radia_mcp.pinn.server"),
    "radia_acoustic": lambda: importlib.import_module("radia_mcp.radia_acoustic.server"),
    "radia_ngsolve": lambda: importlib.import_module("radia_mcp.radia_ngsolve.server"),
    "rna_mec": lambda: importlib.import_module("radia_mcp.rna_mec.server"),
    "streamfunction": lambda: importlib.import_module("radia_mcp.streamfunction.server"),
    "team_benchmark": lambda: importlib.import_module("radia_mcp.team_benchmark.server"),
    "topology_optimization": lambda: importlib.import_module("radia_mcp.topology_optimization.server"),
}


def _server_mcp(module):
    return _SERVER_IMPORTS[module]().mcp


def test_server_import_table_covers_every_pack_domain():
    domains = {m for pack in PACKS.values() for profile in pack["profiles"].values()
               for m in profile}
    assert set(_SERVER_IMPORTS) == domains


@pytest.mark.parametrize("pack", PACKS)
def test_pack_preserves_domain_schemas_prompts_resources(pack):
    async def check():
        app = await CapabilityServer(pack).initialize()
        expected = {}
        prompts = {}
        resources = {}
        for module in modules_for(pack):
            source = _server_mcp(module)
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
        CapabilityServer("radia-mcp.optimize", "typo")


def test_duplicate_tool_is_rejected(monkeypatch):
    monkeypatch.setitem(PACKS, "duplicate", {
        "description": "test", "profiles": {"test": ("bayesian_opt", "bayesian_opt")},
    })
    with pytest.raises(ValueError, match="Duplicate capability"):
        asyncio.run(CapabilityServer("duplicate", "test").initialize())


@pytest.mark.parametrize("pack,profile", [("radia-mcp.analysis", "field")])
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
        app = await CapabilityServer("radia-mcp.optimize", "search").initialize()
        module = modules_for("radia-mcp.optimize", "search")[0]
        source = _server_mcp(module)

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
