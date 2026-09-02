"""A running FastMCP server picks up edited source without a restart."""
from __future__ import annotations

import asyncio
import importlib
import json
import os
import sys
import textwrap
import time

from mcp.server.fastmcp import FastMCP

from radia_mcp._shared import hot_reload


_clock = [time.time() + 2]


def _write(path, body: str) -> None:
    path.write_text(textwrap.dedent(body), encoding="utf-8")
    # Python validates a cached .pyc by whole-second source mtime and size, so
    # two writes of equal length inside one second would reuse the stale
    # bytecode. Give every write its own, strictly later, whole second.
    _clock[0] += 3
    os.utime(path, (_clock[0], _clock[0]))


def _call(mcp: FastMCP, name: str):
    result = asyncio.run(mcp.call_tool(name, {}))
    # Depending on the SDK version this is (content, structured), a dict, or
    # a list of content blocks whose first text block is JSON.
    if isinstance(result, tuple):
        return result[1]
    if isinstance(result, dict):
        return result
    return json.loads(result[0].text)


def test_reload_updates_stale_tools_and_registers_new_ones(tmp_path, monkeypatch):
    pkg = tmp_path / "hotpkg"
    pkg.mkdir()
    _write(pkg / "__init__.py", "")
    _write(pkg / "_helper.py", "VALUE = 'one'\n")
    _write(
        pkg / "tools.py",
        """
        from . import _helper

        def hot_ping() -> dict:
            return {"value": _helper.VALUE}
        """,
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    sys.modules.pop("hotpkg", None)
    tools = importlib.import_module("hotpkg.tools")

    mcp = FastMCP("hot-test")
    mcp.add_tool(tools.hot_ping)
    hot_reload.register_reload_tool(mcp, "hot_reload_code", module_prefix="hotpkg")
    assert _call(mcp, "hot_ping") == {"value": "one"}

    # Edit a dependency and add a tool: both land after one reload.
    _write(pkg / "_helper.py", "VALUE = 'two'  # edited while the server ran\n")
    _write(
        pkg / "tools.py",
        """
        from . import _helper

        def hot_ping() -> dict:
            return {"value": _helper.VALUE}

        def hot_pong() -> dict:
            return {"value": _helper.VALUE + "!"}
        """,
    )
    report = hot_reload.reload_and_refresh(mcp, "hotpkg")

    assert set(report["reloaded"]) >= {"hotpkg._helper", "hotpkg.tools"}
    assert report["errors"] == {}
    assert report["updated"] == ["hot_ping"]
    assert report["added"] == ["hot_pong"]
    assert _call(mcp, "hot_ping") == {"value": "two"}
    assert _call(mcp, "hot_pong") == {"value": "two!"}
    assert "hot_reload_code" in {t.name for t in mcp._tool_manager.list_tools()}

    # Nothing changed since: nothing reloaded, nothing re-registered.
    again = hot_reload.reload_and_refresh(mcp, "hotpkg")
    assert again["reloaded"] == [] and again["updated"] == [] and again["added"] == []

    for name in [m for m in sys.modules if m.startswith("hotpkg")]:
        sys.modules.pop(name, None)


def test_reload_tool_declares_tools_list_changed_capability():
    mcp = FastMCP("cap-test")
    hot_reload.register_reload_tool(mcp, "cap_reload_code", module_prefix="radia_mcp")

    options = mcp._mcp_server.create_initialization_options()

    assert options.capabilities.tools is not None
    assert options.capabilities.tools.listChanged is True
