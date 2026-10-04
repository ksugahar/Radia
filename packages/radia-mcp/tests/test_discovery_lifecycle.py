"""Cold MCP discovery must not acquire a solver or start external work."""

from __future__ import annotations

import asyncio
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

# A fresh interpreter catches import-time work; the journal catches even swallowed
# exceptions. Guards are test-only and do not alter the live server environment.
_GUARD = r"""
import json, importlib, sys, socket
from pathlib import Path
journal = Path(sys.argv[1])
blocked = tuple(json.loads(sys.argv[2]))
def deny(event, args):
    prohibited = event in {"subprocess.Popen", "os.system", "os.startfile", "socket.connect"}
    # Windows asyncio uses the stdlib loopback socketpair for its wakeup pipe.
    # Admit only that exact stdlib frame, never arbitrary localhost connections.
    if event == "socket.connect" and sys._getframe(1).f_code is getattr(socket.socketpair, "__code__", None):
        prohibited = False
    if event == "import":
        name = str(args[0])
        prohibited = any(name == p or name.startswith(p + ".") for p in blocked)
    if prohibited:
        with journal.open("a", encoding="utf-8") as f:
            f.write(json.dumps({"event": event}) + "\n")
        raise RuntimeError("external work attempted during MCP discovery: " + event)
sys.addaudithook(deny)
if sys.argv[3] == "--negative-control":
    try:
        import subprocess
        if sys.argv[4] == "launch":
            subprocess.run([sys.executable, "-c", "raise SystemExit(0)"])
        elif sys.argv[4] == "connect":
            with socket.socket() as connection:
                connection.connect(("127.0.0.1", 9))
        else:
            __import__(blocked[0])
    except RuntimeError:
        pass
else:
    module = sys.argv[3]
    sys.argv = [module]
    importlib.import_module(module).main()
"""

MODULES = [
    "radia_mcp.radia_ngsolve.server",
    "radia_mcp.radia_analysis.server",
    "radia_mcp.matlab.server",
]
BLOCKED = ["ngsolve", "netgen", "matlab.engine"]
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src"


async def _discover(wrapper: Path, journal: Path, module: str) -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(SOURCE), str(ROOT)])
    params = StdioServerParameters(
        command=sys.executable,
        args=[str(wrapper), str(journal), json.dumps(BLOCKED), module],
        cwd=str(ROOT),
        env=env,
    )
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            first = await session.list_tools()
            second = await session.list_tools()
            assert first.tools
            assert [t.name for t in first.tools] == [t.name for t in second.tools]


@pytest.mark.parametrize("module", MODULES)
def test_repeated_discovery_does_not_acquire_solver(tmp_path: Path, module: str) -> None:
    wrapper = tmp_path / "guarded_server.py"
    wrapper.write_text(_GUARD, encoding="utf-8")
    journal = tmp_path / "attempts.jsonl"
    asyncio.run(asyncio.wait_for(_discover(wrapper, journal, module), timeout=60))
    assert not journal.exists(), journal.read_text() if journal.exists() else ""


@pytest.mark.parametrize(
    "attempt,event",
    [
        ("launch", "subprocess.Popen"),
        ("connect", "socket.connect"),
        ("import", "import"),
    ],
)
def test_discovery_guard_detects_swallowed_attempt(
    tmp_path: Path, attempt: str, event: str
) -> None:
    wrapper = tmp_path / "guarded_server.py"
    wrapper.write_text(_GUARD, encoding="utf-8")
    journal = tmp_path / "attempts.jsonl"
    result = subprocess.run(
        [
            sys.executable,
            str(wrapper),
            str(journal),
            json.dumps(["discovery_test_forbidden_module"]),
            "--negative-control",
            attempt,
        ],
        capture_output=True,
        text=True,
        timeout=10,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(journal.read_text()) == {"event": event}
