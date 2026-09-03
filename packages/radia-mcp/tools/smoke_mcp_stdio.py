"""Exercise a radia-mcp server through the real stdio MCP transport.

This is intentionally a protocol probe, not an import test.  It verifies the
same initialize -> tools/list -> status tools/call path used by an MCP client,
including annotations, structured status output, and loaded-source provenance.
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
from pathlib import Path
import sys
from typing import Any

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from radia_mcp.meta.catalog import CATALOG


def _status_payload(result: Any) -> dict[str, Any]:
    structured = getattr(result, "structuredContent", None)
    if isinstance(structured, dict):
        return structured
    for item in getattr(result, "content", ()):
        text = getattr(item, "text", None)
        if not text:
            continue
        try:
            parsed = json.loads(text)
        except (TypeError, ValueError):
            continue
        if isinstance(parsed, dict):
            return parsed
    raise AssertionError("status tools/call returned no structured JSON object")


async def _probe_server(short_name: str) -> dict[str, Any]:
    info = CATALOG[short_name]
    entry_point = str(info["entry_point"])
    subpackage = str(info["subpackage"])
    conventional_status_name = (
        entry_point.removeprefix("mcp-server-").replace("-", "_")
        + "_status"
    )
    environment = os.environ.copy()
    environment["RADIA_MCP_CALL_LOG"] = "0"
    parameters = StdioServerParameters(
        command=sys.executable,
        args=[
            "-c",
            f"from {subpackage}.server import main; main()",
        ],
        env=environment,
    )
    async with stdio_client(parameters) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            initialized = await session.initialize()
            listed = await session.list_tools()
            by_name = {tool.name: tool for tool in listed.tools}
            status_tools = [
                tool.name
                for tool in listed.tools
                if (tool.meta or {}).get("caeai.control_plane") == "status"
            ]
            if len(status_tools) == 1:
                status_name = status_tools[0]
            elif conventional_status_name in by_name:
                status_name = conventional_status_name
            else:
                raise AssertionError(
                    f"{short_name}: could not identify one status tool; "
                    f"marked={status_tools}, conventional="
                    f"{conventional_status_name}"
                )
            incomplete = []
            for name, tool in by_name.items():
                annotations = tool.annotations
                if annotations is None or any(
                    value is None
                    for value in (
                        annotations.readOnlyHint,
                        annotations.destructiveHint,
                        annotations.idempotentHint,
                        annotations.openWorldHint,
                    )
                ):
                    incomplete.append(name)
            if incomplete:
                raise AssertionError(
                    f"{short_name}: incomplete tool annotations: {incomplete}"
                )
            called = await session.call_tool(status_name, {})
            if called.isError:
                raise AssertionError(f"{short_name}: status tools/call failed")
            payload = _status_payload(called)

    if payload.get("schema") != "radia-mcp.server-status.v2":
        raise AssertionError(f"{short_name}: unexpected status schema")
    contract = payload.get("runtime_contract", {})
    if not contract.get("complete"):
        raise AssertionError(f"{short_name}: incomplete runtime contract")
    provenance = payload.get("runtime_provenance", {})
    module_file = provenance.get("module_file")
    if not module_file:
        raise AssertionError(f"{short_name}: status omitted module_file")
    if not provenance.get("module_sha256_at_registration"):
        raise AssertionError(
            f"{short_name}: status omitted registration-time module hash"
        )
    if not provenance.get("module_file_sha256"):
        raise AssertionError(f"{short_name}: status omitted current module hash")
    if provenance.get("source_changed_since_registration"):
        raise AssertionError(
            f"{short_name}: server.py changed after this process registered"
        )
    if contract.get("n_tools") != len(by_name):
        raise AssertionError(
            f"{short_name}: status/list tool-count mismatch "
            f"({contract.get('n_tools')} != {len(by_name)})"
        )
    distribution = provenance.get("distribution", {})
    distribution_version = distribution.get("version")
    if (
        distribution_version
        and distribution_version != "unknown"
        and initialized.serverInfo.version != distribution_version
    ):
        raise AssertionError(
            f"{short_name}: server version {initialized.serverInfo.version} "
            f"does not match distribution {distribution_version}"
        )
    return {
        "server": short_name,
        "server_name": initialized.serverInfo.name,
        "server_version": initialized.serverInfo.version,
        "protocol_version": initialized.protocolVersion,
        "n_tools": len(by_name),
        "status_tool": status_name,
        "module_file": module_file,
        "module_sha256_at_registration": provenance.get(
            "module_sha256_at_registration"
        ),
        "module_file_sha256": provenance.get("module_file_sha256"),
        "source_changed_since_registration": provenance.get(
            "source_changed_since_registration"
        ),
        "distribution": distribution,
        "mcp_sdk_version": provenance.get("mcp_sdk_version"),
        "structured_status": bool(
            getattr(by_name[status_name], "outputSchema", None)
        ),
    }


def probe_server(short_name: str, timeout: float = 60.0) -> dict[str, Any]:
    """Synchronously probe one catalog server through stdio."""
    if short_name not in CATALOG:
        raise KeyError(f"unknown server {short_name!r}")
    return asyncio.run(asyncio.wait_for(_probe_server(short_name), timeout))


def _is_below(path: str, root: str) -> bool:
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except (OSError, ValueError):
        return False


def main() -> int:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--server", choices=sorted(CATALOG))
    group.add_argument("--all", action="store_true")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument(
        "--expect-module-root",
        help="Fail unless each loaded server module is below this directory.",
    )
    args = parser.parse_args()

    selected = sorted(CATALOG) if args.all else [args.server]
    results = []
    for name in selected:
        result = probe_server(name, timeout=args.timeout)
        if args.expect_module_root and not _is_below(
            result["module_file"], args.expect_module_root
        ):
            raise AssertionError(
                f"{name}: loaded {result['module_file']}, expected below "
                f"{args.expect_module_root}"
            )
        results.append(result)
        print(
            f"OK {name}: {result['n_tools']} tools, "
            f"MCP {result['protocol_version']}",
            file=sys.stderr,
        )
    print(json.dumps(results, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
