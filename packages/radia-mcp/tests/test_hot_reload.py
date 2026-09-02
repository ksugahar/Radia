"""A running FastMCP server picks up edited source without a restart."""
from __future__ import annotations

import asyncio
import importlib
import json
import os
import sys
import textwrap
import time
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
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
    try:
        return json.loads(result[0].text)
    except json.JSONDecodeError:
        return result[0].text


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
    assert report["removed"] == []
    assert report["added_tools_need_reconnect_for_server_policy"] is True
    assert _call(mcp, "hot_ping") == {"value": "two"}
    assert _call(mcp, "hot_pong") == {"value": "two!"}
    added = mcp._tool_manager._tools["hot_pong"]
    assert added.annotations.destructiveHint is True
    assert added.annotations.readOnlyHint is False
    assert "hot_reload_code" in {t.name for t in mcp._tool_manager.list_tools()}

    # Nothing changed since: nothing reloaded, nothing re-registered.
    again = hot_reload.reload_and_refresh(mcp, "hotpkg")
    assert again["reloaded"] == [] and again["updated"] == [] and again["added"] == []

    # Removing a source function removes its stale registered tool too.
    _write(
        pkg / "tools.py",
        """
        from . import _helper

        def hot_ping() -> dict:
            return {"value": _helper.VALUE}
        """,
    )
    removed = hot_reload.reload_and_refresh(mcp, "hotpkg")
    assert removed["removed"] == ["hot_pong"]
    assert "hot_pong" not in {t.name for t in mcp._tool_manager.list_tools()}

    for name in [m for m in sys.modules if m.startswith("hotpkg")]:
        sys.modules.pop(name, None)


def test_registered_mtime_baseline_handles_a_file_server_clock_skew(
    tmp_path, monkeypatch
):
    pkg = tmp_path / "skewpkg"
    pkg.mkdir()
    _write(pkg / "__init__.py", "")
    _write(pkg / "tools.py", "def skew_ping(): return 'one'\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    tools = importlib.import_module("skewpkg.tools")
    mcp = FastMCP("skew-test")
    mcp.add_tool(tools.skew_ping)
    hot_reload.register_reload_tool(mcp, "skew_reload_code", module_prefix="skewpkg")

    path = pkg / "tools.py"
    path.write_text("def skew_ping(): return 'two'\n", encoding="utf-8")
    remote_past = time.time() - 3600
    os.utime(path, (remote_past, remote_past))

    report = hot_reload.reload_and_refresh(mcp, "skewpkg")

    assert "skewpkg.tools" in report["reloaded"]
    assert _call(mcp, "skew_ping") == "two"

    for name in [m for m in sys.modules if m.startswith("skewpkg")]:
        sys.modules.pop(name, None)


def test_reload_tool_declares_tools_list_changed_capability():
    mcp = FastMCP("cap-test")
    hot_reload.register_reload_tool(mcp, "cap_reload_code", module_prefix="radia_mcp")

    options = mcp._mcp_server.create_initialization_options()

    assert options.capabilities.tools is not None
    assert options.capabilities.tools.listChanged is True


async def _probe_reload_notification_over_stdio() -> dict:
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(package_root / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    notifications: list[str] = []

    async def handle_message(message) -> None:
        root = getattr(message, "root", None)
        if root is not None:
            notifications.append(type(root).__name__)

    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "radia_mcp.grant_writing.server"],
        cwd=str(package_root),
        env=env,
    )
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(
            read_stream,
            write_stream,
            message_handler=handle_message,
        ) as session:
            initialized = await session.initialize()
            called = await session.call_tool("grant_writing_reload_code", {})
            await asyncio.sleep(0.05)
            return {
                "list_changed": initialized.capabilities.tools.listChanged,
                "payload": json.loads(called.content[0].text),
                "notifications": notifications,
            }


def test_reload_tool_notifies_a_real_stdio_client():
    result = asyncio.run(_probe_reload_notification_over_stdio())

    assert result["list_changed"] is True
    assert result["payload"]["client_notified"] is True
    assert "ToolListChangedNotification" in result["notifications"]


def test_reload_failure_restores_previous_module_and_tool(tmp_path, monkeypatch):
    pkg = tmp_path / "rollbackpkg"
    pkg.mkdir()
    _write(pkg / "__init__.py", "")
    _write(pkg / "tools.py", "def rollback_ping(): return 'one'\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    tools = importlib.import_module("rollbackpkg.tools")
    mcp = FastMCP("rollback-test")
    mcp.add_tool(tools.rollback_ping)

    _write(pkg / "tools.py", "def rollback_ping(:\n    return 'broken'\n")
    syntax = hot_reload.reload_and_refresh(mcp, "rollbackpkg")
    assert syntax["errors"]
    assert syntax["rolled_back"] is False
    assert _call(mcp, "rollback_ping") == "one"

    _write(
        pkg / "tools.py",
        """
        def rollback_ping():
            return 'partial'
        raise RuntimeError('import failed')
        """,
    )
    runtime = hot_reload.reload_and_refresh(mcp, "rollbackpkg")
    assert runtime["errors"]
    assert runtime["rolled_back"] is True
    assert _call(mcp, "rollback_ping") == "one"

    _write(pkg / "tools.py", "def rollback_ping(): return 'two'\n")
    fixed = hot_reload.reload_and_refresh(mcp, "rollbackpkg")
    assert fixed["errors"] == {}
    assert _call(mcp, "rollback_ping") == "two"

    for name in [m for m in sys.modules if m.startswith("rollbackpkg")]:
        sys.modules.pop(name, None)


def test_server_module_that_owns_live_mcp_requires_restart(tmp_path, monkeypatch):
    pkg = tmp_path / "ownerpkg"
    pkg.mkdir()
    _write(pkg / "__init__.py", "")
    _write(pkg / "server.py", "def owner_ping(): return 'one'\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    server = importlib.import_module("ownerpkg.server")
    mcp = FastMCP("owner-test")
    server.LIVE_MCP = mcp
    mcp.add_tool(server.owner_ping)

    _write(pkg / "server.py", "def owner_ping(): return 'two'\n")
    report = hot_reload.reload_and_refresh(mcp, "ownerpkg")
    assert "ownerpkg.server" not in report["reloaded"]
    assert report["restart_required"] == ["ownerpkg.server"]
    assert _call(mcp, "owner_ping") == "one"

    sys.modules.pop("ownerpkg.server", None)
    sys.modules.pop("ownerpkg", None)


def test_tool_refresh_failure_rolls_back_code_and_registry(tmp_path, monkeypatch):
    pkg = tmp_path / "registrypkg"
    pkg.mkdir()
    _write(pkg / "__init__.py", "")
    _write(pkg / "tools.py", "def registry_ping(): return 'one'\n")
    monkeypatch.syspath_prepend(str(tmp_path))
    tools = importlib.import_module("registrypkg.tools")
    mcp = FastMCP("registry-test")
    mcp.add_tool(tools.registry_ping)
    original_add_tool = mcp.add_tool

    _write(pkg / "tools.py", "def registry_ping(): return 'two'\n")

    def fail_add(*args, **kwargs):
        raise RuntimeError("registry rejected replacement")

    monkeypatch.setattr(mcp, "add_tool", fail_add)
    report = hot_reload.reload_and_refresh(mcp, "registrypkg")
    assert report["rolled_back"] is True
    assert "tool_registry" in report["errors"]
    assert _call(mcp, "registry_ping") == "one"

    monkeypatch.setattr(mcp, "add_tool", original_add_tool)
    fixed = hot_reload.reload_and_refresh(mcp, "registrypkg")
    assert fixed["errors"] == {}
    assert _call(mcp, "registry_ping") == "two"

    sys.modules.pop("registrypkg.tools", None)
    sys.modules.pop("registrypkg", None)
