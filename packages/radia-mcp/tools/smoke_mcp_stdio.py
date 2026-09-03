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
from mcp.types import Implementation

from radia_mcp.meta.catalog import CATALOG


_PROBE_CLIENT = {
    "name": "cae-lab-radia-contract-probe",
    "title": "CAE-AI Lab Runtime Contract Probe",
    "version": "1.0",
}


def _probe_artifact() -> dict[str, Any]:
    return {
        "schema": "cae-ai-lab.solver-artifact-identity.v1",
        "artifact_id": "protocol-probe-001",
        "created_at_utc": "2026-09-03T00:00:00Z",
        "producer": {"name": "protocol-probe", "version": "1.0"},
        "solver": {"name": "analytical-fixture", "version": "1.0"},
        "status": "complete",
        "coordinate_system": "cartesian-right-handed",
        "unit_system": "SI",
        "input_sha256": "a" * 64,
        "result_sha256": "b" * 64,
        "total_compute_seconds": 0.01,
        "timing_breakdown_s": {"verify": 0.01},
    }


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
    local_src = str(Path(__file__).resolve().parents[1] / "src")
    inherited_pythonpath = environment.get("PYTHONPATH", "")
    environment["PYTHONPATH"] = os.pathsep.join(
        part for part in (local_src, inherited_pythonpath) if part
    )
    parameters = StdioServerParameters(
        command=sys.executable,
        args=[
            "-c",
            f"from {subpackage}.server import main; main()",
        ],
        env=environment,
    )
    artifact_gate = None
    async with stdio_client(parameters) as (read_stream, write_stream):
        async with ClientSession(
            read_stream,
            write_stream,
            client_info=Implementation(**_PROBE_CLIENT),
        ) as session:
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
            called_again = await session.call_tool(status_name, {})
            if called_again.isError:
                raise AssertionError(
                    f"{short_name}: second status tools/call failed"
                )
            payload_again = _status_payload(called_again)
            if short_name == "meta":
                gate_name = "radia_mcp_validate_solver_artifact"
                if gate_name not in by_name:
                    raise AssertionError(f"{short_name}: artifact gate missing")
                accepted = await session.call_tool(
                    gate_name, {"artifact": _probe_artifact()}
                )
                stale = _probe_artifact()
                stale["result_sha256"] = "stale"
                rejected = await session.call_tool(
                    gate_name, {"artifact": stale}
                )
                accepted_payload = _status_payload(accepted)
                rejected_payload = _status_payload(rejected)
                if not accepted_payload.get("accepted"):
                    raise AssertionError("meta: complete artifact was rejected")
                if rejected_payload.get("accepted"):
                    raise AssertionError("meta: stale artifact was accepted")
                artifact_gate = {
                    "positive": accepted_payload.get("status"),
                    "negative": rejected_payload.get("status"),
                }

    if payload.get("schema") != "radia-mcp.server-status.v2":
        raise AssertionError(f"{short_name}: unexpected status schema")
    contract = payload.get("runtime_contract", {})
    if not contract.get("complete"):
        raise AssertionError(f"{short_name}: incomplete runtime contract")
    maturity = contract.get("capability_maturity", {})
    if sum(int(value) for value in maturity.values()) != len(by_name):
        raise AssertionError(f"{short_name}: incomplete maturity metadata")
    if contract.get("invalid_maturity"):
        raise AssertionError(f"{short_name}: invalid maturity metadata")
    if contract.get("artifact_identity_contract", {}).get(
        "canonical_container"
    ) != ".hdf5":
        raise AssertionError(f"{short_name}: artifact contract missing")
    connection = payload_again.get("client_connection", {})
    if connection.get("connection_count") != 1:
        raise AssertionError(
            f"{short_name}: client identity was not captured exactly once"
        )
    if connection.get("latest") != _PROBE_CLIENT:
        raise AssertionError(
            f"{short_name}: client identity differs from initialize payload"
        )
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
        "client_connection": connection,
        "capability_maturity": maturity,
        "artifact_gate": artifact_gate,
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
