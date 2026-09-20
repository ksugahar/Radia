"""Audit heterogeneous MCP servers through their real stdio transport.

The JSON configuration is intentionally external so this public utility does
not encode private product names or machine paths. Each server is initialized
with one explicit client identity, listed, and asked for its passive status
twice. No solver or application tool is called.
"""

from __future__ import annotations

import argparse
import asyncio
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import time
from typing import Any

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.types import Implementation


SCHEMA = "cae-ai-lab.mcp-fleet-audit.v1"
MATURITY_LEVELS = [
    "knowledge_only",
    "artifact_gate",
    "solver_ready",
    "live_verified",
    "cross_validated",
]
CLIENT_INFO = {
    "name": "cae-lab-fleet-contract-probe",
    "title": "CAE-AI Lab MCP Fleet Contract Probe",
    "version": "1.0",
}


def _git_state(cwd: str | None) -> dict[str, Any] | None:
    if not cwd:
        return None

    def run(*args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", "-C", cwd, *args],
            capture_output=True,
            check=False,
            encoding="utf-8",
            errors="replace",
            text=True,
            timeout=10,
        )

    head = run("rev-parse", "HEAD")
    if head.returncode != 0:
        return None
    branch = run("branch", "--show-current")
    tracked = run("status", "--porcelain", "--untracked-files=no")
    untracked = run("status", "--porcelain", "--untracked-files=all")
    untracked_count = sum(
        1 for line in untracked.stdout.splitlines() if line.startswith("?? ")
    )
    return {
        "commit": head.stdout.strip(),
        "branch": branch.stdout.strip() if branch.returncode == 0 else "",
        "tracked_clean": tracked.returncode == 0 and not tracked.stdout.strip(),
        "untracked_count": untracked_count,
    }


def _payload(result: Any) -> dict[str, Any]:
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
    raise AssertionError("tools/call returned no JSON object")


def _contract_sections(contract: dict[str, Any]) -> tuple[dict, dict]:
    maturity = contract.get("capability_maturity_contract")
    if not isinstance(maturity, dict):
        maturity = contract.get("capability_maturity", {})
    artifact = contract.get("artifact_identity_contract")
    if not isinstance(artifact, dict):
        artifact = contract.get("artifact_identity", {})
    return maturity, artifact


async def _probe(entry: dict[str, Any]) -> dict[str, Any]:
    server_id = str(entry["id"])
    mode = str(entry.get("mode", "lab_contract"))
    status_tool = str(entry.get("status_tool", ""))
    source_revision = _git_state(str(entry.get("cwd", "")))
    env = os.environ.copy()
    env.update({str(k): str(v) for k, v in entry.get("env", {}).items()})
    params = StdioServerParameters(
        command=str(entry["command"]),
        args=[str(value) for value in entry.get("args", [])],
        cwd=str(entry.get("cwd")) if entry.get("cwd") else None,
        env=env,
    )
    started = time.perf_counter()
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(
            read_stream,
            write_stream,
            client_info=Implementation(**CLIENT_INFO),
        ) as session:
            initialized = await session.initialize()
            listed = await session.list_tools()
            tools = {tool.name: tool for tool in listed.tools}
            if mode == "upstream_reference":
                checks = {
                    "instructions_present": bool(initialized.instructions),
                    "tools_present": bool(tools),
                    "titles_complete": all(tool.title for tool in tools.values()),
                    "annotations_complete": all(
                        tool.annotations is not None for tool in tools.values()
                    ),
                }
                return {
                    "id": server_id,
                    "mode": mode,
                    "ok": all(checks.values()),
                    "checks": checks,
                    "server_name": initialized.serverInfo.name,
                    "server_version": initialized.serverInfo.version,
                    "protocol_version": initialized.protocolVersion,
                    "n_tools": len(tools),
                    "elapsed_seconds": round(time.perf_counter() - started, 6),
                    "source_revision": source_revision,
                    "note": "Upstream reference is not required to expose the lab contract.",
                }
            if status_tool not in tools:
                raise AssertionError(f"{server_id}: status tool {status_tool!r} missing")
            missing_titles = [name for name, tool in tools.items() if not tool.title]
            missing_annotations = [
                name
                for name, tool in tools.items()
                if tool.annotations is None
                or any(
                    getattr(tool.annotations, hint) is None
                    for hint in (
                        "readOnlyHint",
                        "destructiveHint",
                        "idempotentHint",
                        "openWorldHint",
                    )
                )
            ]
            invalid_maturity = [
                name
                for name, tool in tools.items()
                if (tool.meta or {}).get("caeai.maturity") not in MATURITY_LEVELS
            ]
            first = await session.call_tool(status_tool, {})
            second = await session.call_tool(status_tool, {})
            if first.isError or second.isError:
                raise AssertionError(f"{server_id}: passive status call failed")
            status = _payload(second)

    contract = status.get("runtime_contract", {})
    if not isinstance(contract, dict):
        raise AssertionError(f"{server_id}: runtime contract missing")
    maturity, artifact = _contract_sections(contract)
    client = status.get("client_connection", {})
    checks = {
        "instructions_present": bool(initialized.instructions),
        "tools_present": bool(tools),
        "titles_complete": not missing_titles,
        "annotations_complete": not missing_annotations,
        "maturity_metadata_complete": not invalid_maturity,
        "maturity_vocabulary_matches": maturity.get("levels") == MATURITY_LEVELS,
        "artifact_schema_matches": artifact.get("schema")
        == "cae-ai-lab.solver-artifact-identity.v1",
        "canonical_container_hdf5": artifact.get("canonical_container") == ".hdf5",
        "client_captured_once": client.get("connection_count") == 1,
        "client_identity_matches": client.get("latest") == CLIENT_INFO,
        "source_tracked_clean": (
            source_revision is None or source_revision["tracked_clean"]
        ),
    }
    return {
        "id": server_id,
        "ok": all(checks.values()),
        "checks": checks,
        "server_name": initialized.serverInfo.name,
        "server_version": initialized.serverInfo.version,
        "protocol_version": initialized.protocolVersion,
        "n_tools": len(tools),
        "status_tool": status_tool,
        "runtime_contract_schema": contract.get("schema"),
        "elapsed_seconds": round(time.perf_counter() - started, 6),
        "missing_titles": missing_titles,
        "missing_annotations": missing_annotations,
        "invalid_maturity": invalid_maturity,
        "client_connection": client,
        "source_revision": source_revision,
    }


async def audit(config: dict[str, Any]) -> dict[str, Any]:
    entries = config.get("servers")
    if not isinstance(entries, list) or not entries:
        raise ValueError("config.servers must be a non-empty list")
    started = time.perf_counter()
    results = []
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("each server entry must be an object")
        timeout = float(entry.get("timeout_seconds", 60.0))
        try:
            result = await asyncio.wait_for(_probe(entry), timeout=timeout)
        except Exception as exc:
            result = {
                "id": str(entry.get("id", "unknown")),
                "ok": False,
                "error": f"{type(exc).__name__}: {exc}",
            }
        results.append(result)
    return {
        "schema": SCHEMA,
        "artifact_id": str(config.get("artifact_id", "mcp-fleet-audit")),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "producer": {"name": "audit_mcp_fleet", "version": "1.0"},
        "status": "passed" if all(row["ok"] for row in results) else "failed",
        "all_passed": all(row["ok"] for row in results),
        "n_servers": len(results),
        "total_compute_seconds": round(time.perf_counter() - started, 6),
        "servers": results,
    }


def write_report(report: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() in {".h5", ".hdf5"}:
        import h5py

        with h5py.File(path, "w") as handle:
            handle.attrs["schema"] = report["schema"]
            handle.attrs["artifact_id"] = report["artifact_id"]
            handle.attrs["created_at_utc"] = report["created_at_utc"]
            handle.attrs["producer"] = report["producer"]["name"]
            handle.attrs["producer_version"] = report["producer"]["version"]
            handle.attrs["unit_system"] = "not_applicable"
            handle.attrs["coordinate_system"] = "not_applicable"
            handle.attrs["total_compute_seconds"] = report[
                "total_compute_seconds"
            ]
            handle.create_dataset(
                "report_json",
                data=json.dumps(report, ensure_ascii=False, sort_keys=True),
                dtype=h5py.string_dtype(encoding="utf-8"),
            )
            servers = handle.create_group("servers")
            for row in report["servers"]:
                key = re.sub(r"[^A-Za-z0-9_.-]+", "_", row["id"])
                group = servers.create_group(key)
                group.attrs["ok"] = bool(row["ok"])
                group.create_dataset(
                    "report_json",
                    data=json.dumps(row, ensure_ascii=False, sort_keys=True),
                    dtype=h5py.string_dtype(encoding="utf-8"),
                )
        return
    path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    report = asyncio.run(audit(config))
    if args.output:
        write_report(report, args.output)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
