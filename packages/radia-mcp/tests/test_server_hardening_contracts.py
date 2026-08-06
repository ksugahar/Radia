"""Contract tests for the shared MCP-server hardening layer
(radia_mcp.common.server_hardening) across the cubit and build123d
servers.

Locks the holes a rename/typo would open silently:
* a name listed in an explicit classification set that matches NO
  registered tool would silently stop protecting anything (the tool it
  meant to mark destructive would default to READONLY);
* server instructions could be dropped without any test noticing;
* the gate-hiding env could stop working.
"""

import json
import os
import subprocess
import sys

import pytest

from radia_mcp.build123d import server as b3d_server
from radia_mcp.cubit import server as cubit_server
from radia_mcp.common.server_hardening import (
    ANN_DESTRUCTIVE,
    ANN_READONLY,
    error_payload,
)


def _tool_names(mcp) -> set:
    return set(mcp._tool_manager._tools)


# ---------------------------------------------------------------------------
# Classification sets must reference only REAL tool names (typo guard)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("server_mod", [cubit_server, b3d_server],
                         ids=["cubit", "build123d"])
def test_classification_sets_match_registered_tools(server_mod):
    names = _tool_names(server_mod.mcp)
    for set_name in ("_DESTRUCTIVE_TOOLS", "_WRITING_TOOLS", "_WEB_TOOLS"):
        listed = getattr(server_mod, set_name)
        ghosts = listed - names
        assert not ghosts, (
            f"{server_mod.__name__}.{set_name} lists tools that do not "
            f"exist (rename/typo -- the intended tool silently degrades "
            f"to READONLY): {sorted(ghosts)}")


@pytest.mark.parametrize("server_mod", [cubit_server, b3d_server],
                         ids=["cubit", "build123d"])
def test_no_unclassified_tools(server_mod):
    assert server_mod._UNCLASSIFIED_TOOLS == []


def test_critical_tools_keep_their_preset():
    """The safety-critical classifications must never drift."""
    cub = cubit_server.mcp._tool_manager
    b3d = b3d_server.mcp._tool_manager
    assert cub.get_tool("cubit_exec").annotations == ANN_DESTRUCTIVE
    assert cub.get_tool("cubit_session_shutdown").annotations \
        == ANN_DESTRUCTIVE
    assert cub.get_tool("cubit_probe").annotations == ANN_READONLY
    assert cub.get_tool("cubit_doctor").annotations == ANN_READONLY
    assert b3d.get_tool("execute_build123d").annotations == ANN_DESTRUCTIVE
    assert b3d.get_tool("build123d_probe").annotations == ANN_READONLY
    assert b3d.get_tool("build123d_doctor").annotations == ANN_READONLY


# ---------------------------------------------------------------------------
# Server instructions are declared and carry the operating doctrine
# ---------------------------------------------------------------------------

def test_server_instructions_declared():
    cub = cubit_server.mcp.instructions or ""
    b3d = b3d_server.mcp.instructions or ""
    assert "cubit_probe" in cub and "check_vol" in cub.replace("-", "_")
    assert '"kind"' in cub or "kind=" in cub
    assert "label" in b3d and "build123d_probe" in b3d
    assert '"kind"' in b3d or "kind=" in b3d


# ---------------------------------------------------------------------------
# Gate hiding works per server (fresh interpreter: classification runs
# at import time, so the env must be set before the module loads)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("module,env", [
    ("radia_mcp.cubit.server", "RADIA_MCP_CUBIT_GATES"),
    ("radia_mcp.build123d.server", "RADIA_MCP_BUILD123D_GATES"),
], ids=["cubit", "build123d"])
def test_gate_env_hides_gate_tools(module, env):
    code = (
        "import asyncio, json, sys\n"
        f"import {module} as s\n"
        "tools = asyncio.run(s.mcp.list_tools())\n"
        "print(json.dumps([t.name for t in tools]))\n"
    )
    env_full = dict(os.environ, **{env: "0"})
    out = subprocess.run([sys.executable, "-c", code], env=env_full,
                         capture_output=True, text=True, timeout=240)
    assert out.returncode == 0, out.stderr[-800:]
    names = json.loads(out.stdout.strip().splitlines()[-1])
    assert names, "no tools listed"
    assert not [n for n in names if n.endswith("_gate")]


# ---------------------------------------------------------------------------
# Shared error contract
# ---------------------------------------------------------------------------

def test_call_log_rotation(tmp_path):
    from radia_mcp.common.server_hardening import rotate_if_large

    log = tmp_path / "calls.jsonl"
    log.write_text("x" * 100, encoding="utf-8")
    assert rotate_if_large(log, cap_bytes=1000) is False
    assert rotate_if_large(log, cap_bytes=50) is True
    assert not log.exists()
    assert (tmp_path / "calls.jsonl.1").read_text(encoding="utf-8") \
        == "x" * 100
    # a second rotation REPLACES the old .1 (bounded at 2x cap on disk)
    log.write_text("y" * 100, encoding="utf-8")
    assert rotate_if_large(log, cap_bytes=50) is True
    assert (tmp_path / "calls.jsonl.1").read_text(encoding="utf-8") \
        == "y" * 100


def test_error_payload_kind_classification():
    p = error_payload("rpc", "RLM license refused",
                      environment_needles=("license",),
                      log="C:/logs")
    assert p["kind"] == "environment" and p["log"] == "C:/logs"
    p2 = error_payload("input", "bad argument",
                       environment_needles=("license",), log="C:/logs")
    assert p2["kind"] == "input" and "log" not in p2
    p3 = error_payload("x", "boom", kind="internal", log="L")
    assert p3["kind"] == "internal" and p3["log"] == "L"
