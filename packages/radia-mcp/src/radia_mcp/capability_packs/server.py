"""Compose existing FastMCP domains through public SDK APIs in one process.

No child transports are started. Calls retain domain validation, annotations,
structured/image results, prompts, resources, hot reload and call logging.
Only stateless server lifecycles are admitted by the explicit pack manifest;
external session owners (Cubit, MATLAB, Mathematica) stay standalone.
"""

from __future__ import annotations

import argparse
import asyncio
import importlib
import json
import sys

from mcp import types
from mcp.server.lowlevel import Server
from mcp.server.stdio import stdio_server

from . import PACKS, modules_for
from .. import __version__
from ..common.status import build_status_payload, _runtime_provenance
from ..common.mcp_contract import audit_tool_definitions, SCHEMA


class CapabilityServer:
    """One immutable discovery surface per startup profile."""

    def __init__(self, pack: str, profile: str = "all"):
        self.pack = pack
        self.profile = profile
        self.modules = modules_for(pack, profile)
        self.tools = {}
        self.prompts = {}
        self.resources = {}
        self.registration_provenance = _runtime_provenance("radia_mcp.capability_packs")
        self.server = Server(
            f"mcp-server-{pack}",
            version=__version__,
            instructions=(
                f"{PACKS[pack]['description']}. Active profile: {profile}. "
                "Call capability_pack_status first. Domain status tools describe their "
                "own capabilities. Preserve production/education and solver ownership; "
                "a successful tool call is not evidence of numerical agreement."
            ),
        )
        self.server.list_tools()(self.list_tools)
        self.server.call_tool()(self.call_tool)
        self.server.list_prompts()(self.list_prompts)
        self.server.get_prompt()(self.get_prompt)
        self.server.list_resources()(self.list_resources)
        self.server.read_resource()(self.read_resource)

    async def initialize(self):
        """Load selected domains, rejecting any ambiguous public ownership."""
        if self.tools:
            raise RuntimeError("Capability server is already initialized")
        instructions = []
        for module in self.modules:
            source = importlib.import_module(f"radia_mcp.{module}.server").mcp
            if await source.list_resource_templates():
                raise ValueError(f"{module}: resource templates require an explicit routing adapter")
            for registry, entries, key in (
                (self.tools, await source.list_tools(), "name"),
                (self.prompts, await source.list_prompts(), "name"),
                (self.resources, await source.list_resources(), "uri"),
            ):
                for entry in entries:
                    name = str(getattr(entry, key))
                    if name == "capability_pack_status" or name in registry:
                        raise ValueError(f"Duplicate capability name {name!r} in {module}")
                    registry[name] = (entry, source, module)
            if source.instructions:
                instructions.append(f"[{module}] {source.instructions}")
        self.server.instructions += "\n" + "\n".join(instructions)
        return self

    def status(self):
        provenance = _runtime_provenance("radia_mcp.capability_packs")
        registered = self.registration_provenance.get("module_file_sha256")
        provenance.update(module_sha256_at_registration=registered,
                          source_changed_since_registration=registered != provenance.get("module_file_sha256"))
        definitions = {name: entry for name, (entry, _, _) in self.tools.items()}
        definitions["capability_pack_status"] = self.status_tool()
        payload = build_status_payload(
            f"mcp-server-{self.pack}", PACKS[self.pack]["description"],
            "radia_mcp.capability_packs", mcp_tools=list(definitions),
            runtime_contract=audit_tool_definitions(definitions),
            runtime_provenance=provenance,
        )
        payload.update({
            "pack": self.pack, "profile": self.profile,
            "modules": list(self.modules),
            "profiles": ["all", *[p for p in PACKS[self.pack]["profiles"] if p != "all"]],
            "owners": {name: owner for name, (_, _, owner) in self.tools.items()},
            "transport": "stdio", "child_mcp_processes": False,
            "version": __version__, "python_executable": sys.executable,
            "source": __file__,
        })
        return payload

    @staticmethod
    def status_tool():
        return types.Tool(
                name="capability_pack_status",
                title="Capability Pack Status",
                description="Describe the active capability profile and exact tool ownership.",
                inputSchema={"type": "object", "properties": {}, "additionalProperties": False},
                outputSchema={"type": "object"},
                _meta={"caeai.contract": SCHEMA, "caeai.control_plane": "status",
                      "caeai.annotation_source": "explicit"},
                annotations=types.ToolAnnotations(readOnlyHint=True, destructiveHint=False,
                                                  idempotentHint=True, openWorldHint=False),
            )

    async def list_tools(self):
        return [
            self.status_tool(),
            *[entry for entry, _, _ in self.tools.values()],
        ]

    async def call_tool(self, name, arguments):
        if name == "capability_pack_status":
            return self.status()
        if name not in self.tools:
            raise ValueError(f"Unknown tool: {name}")
        return await self.tools[name][1].call_tool(name, arguments)

    async def list_prompts(self):
        return [entry for entry, _, _ in self.prompts.values()]

    async def get_prompt(self, name, arguments):
        return await self.prompts[name][1].get_prompt(name, arguments)

    async def list_resources(self):
        return [entry for entry, _, _ in self.resources.values()]

    async def read_resource(self, uri):
        return await self.resources[str(uri)][1].read_resource(uri)


def main(pack: str):
    """Start a fixed profile; changing the surface requires a process restart."""
    parser = argparse.ArgumentParser(description=PACKS[pack]["description"])
    parser.add_argument("--profile", default="all", choices=list(dict.fromkeys(["all", *PACKS[pack]["profiles"]])))
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()

    async def run():
        app = await CapabilityServer(pack, args.profile).initialize()
        if args.selftest:
            print(json.dumps(app.status(), ensure_ascii=True))
            print("PASSED")
            return
        async with stdio_server() as (read, write):
            await app.server.run(read, write, app.server.create_initialization_options())

    asyncio.run(run())
