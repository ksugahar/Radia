"""Reload changed radia_mcp modules and re-register their tools without a restart.

Every radia-mcp server is an editable install, and every server still had to
be restarted after a code change (2026-09-02). An editable install only tells
Python which file to read at import time; a running server has already
imported its modules and holds the old function objects, FastMCP's tool
registry points at those objects, and the client froze the tool list when it
connected. Three things therefore have to happen inside the running process:

1. reload the modules whose source changed on disk (dependencies first);
2. replace every registered tool whose function object is now stale, and
   register callables that appeared since start-up;
3. tell the client that the tool list changed
   (``notifications/tools/list_changed``), which the client only honours when
   the server declared ``tools.listChanged`` at initialisation.

``register_reload_tool`` does all three for a FastMCP instance. Reload order is
"deeper module first, ``tools`` and packages last", and the changed set is
reloaded twice so a module that imports a sibling reloaded after it still
ends up bound to the new objects. Module-level state (caches, sessions) is
re-created by a reload; a server that must keep such state across reloads
should hold it in a module that does not change.
"""
from __future__ import annotations

import importlib
import os
import sys
import time
from types import ModuleType
from typing import Any

# Imported at module level on purpose: FastMCP evaluates the reload tool's
# annotations in this module's globals, and ``from __future__ import
# annotations`` turns ``Context`` into a string to look up here.
from mcp.server.fastmcp import Context

# Source mtime each module was last (re)loaded at. A module not yet seen is
# compared with the process start instead. Comparing with the recorded value
# rather than with "the last reload time" keeps this correct when a file
# server's clock runs ahead or behind this machine's.
_seen_mtimes: dict[str, float] = {}


def _source_mtime(module: ModuleType) -> float | None:
    path = getattr(module, "__file__", None)
    if not path or not path.endswith(".py"):
        return None
    try:
        return os.path.getmtime(path)
    except OSError:
        return None


def _reload_order(name: str) -> tuple[int, int, str]:
    parts = name.split(".")
    # Packages (``__init__``) and ``tools`` modules import everything else in
    # their tree, so they go last; deeper modules go first.
    last = parts[-1]
    tail_rank = 2 if last in {"tools", "server"} else 0
    return (-len(parts), tail_rank, name)


def reload_changed_modules(prefix: str = "radia_mcp") -> dict[str, Any]:
    """Reload every imported module under ``prefix`` whose file changed.

    The first call reloads what changed since the process started; later
    calls reload what changed since the previous call. Returns the module
    names reloaded and any reload errors (a module that fails to import is
    reported, and the previous version stays in place).
    """
    changed: dict[str, float] = {}
    for name, module in list(sys.modules.items()):
        if module is None or not (name == prefix or name.startswith(prefix + ".")):
            continue
        mtime = _source_mtime(module)
        if mtime is None:
            continue
        seen = _seen_mtimes.get(name)
        if (seen is None and mtime > _PROCESS_START) or (seen is not None and mtime != seen):
            changed[name] = mtime
    order = sorted(changed, key=_reload_order)

    importlib.invalidate_caches()
    reloaded: list[str] = []
    errors: dict[str, str] = {}
    for _pass in range(2):
        for name in order:
            module = sys.modules.get(name)
            if module is None:
                continue
            try:
                importlib.reload(module)
            except Exception as exc:  # pragma: no cover - reported, not raised
                errors[name] = f"{type(exc).__name__}: {exc}"
            else:
                if name not in reloaded:
                    reloaded.append(name)
    for name, mtime in changed.items():
        if name not in errors:
            _seen_mtimes[name] = mtime
    return {"reloaded": reloaded, "errors": errors}


def _common_tool_prefix(names: list[str]) -> str:
    if not names:
        return ""
    prefix = os.path.commonprefix(names)
    # Cut back to the last underscore so ``grant_writing_check_kanji_ratio``
    # and ``grant_writing_check_misuse`` share ``grant_writing_`` rather than
    # ``grant_writing_check_``.
    if len(names) == 1:
        head = names[0].split("_")
        return head[0] + "_" if len(head) > 1 else names[0]
    return prefix[: prefix.rfind("_") + 1] if "_" in prefix else prefix


def refresh_tools(mcp: Any, module_prefix: str = "radia_mcp") -> dict[str, Any]:
    """Re-register stale tools and pick up callables added since start-up."""
    manager = mcp._tool_manager
    updated: list[str] = []
    added: list[str] = []
    skipped: list[str] = []
    names_by_module: dict[str, list[str]] = {}

    for tool in list(manager.list_tools()):
        fn = getattr(tool, "fn", None)
        module_name = getattr(fn, "__module__", "") or ""
        if not module_name.startswith(module_prefix):
            continue
        module = sys.modules.get(module_name)
        if module is None:
            skipped.append(tool.name)
            continue
        names_by_module.setdefault(module_name, []).append(tool.name)
        fresh = getattr(module, fn.__name__, None)
        if fresh is None or fresh is fn or not callable(fresh):
            continue
        mcp.remove_tool(tool.name)
        mcp.add_tool(
            fresh,
            name=tool.name,
            title=getattr(tool, "title", None),
            annotations=getattr(tool, "annotations", None),
        )
        updated.append(tool.name)

    registered = {tool.name for tool in manager.list_tools()}
    for module_name, names in names_by_module.items():
        module = sys.modules.get(module_name)
        prefix = _common_tool_prefix(sorted(names))
        if module is None or not prefix:
            continue
        for attr in dir(module):
            if not attr.startswith(prefix) or attr in registered:
                continue
            candidate = getattr(module, attr)
            if callable(candidate) and getattr(candidate, "__module__", "") == module_name:
                mcp.add_tool(candidate)
                registered.add(attr)
                added.append(attr)
    return {"updated": updated, "added": added, "skipped": skipped}


def reload_and_refresh(mcp: Any, module_prefix: str = "radia_mcp") -> dict[str, Any]:
    report = reload_changed_modules(module_prefix)
    report.update(refresh_tools(mcp, module_prefix))
    return report


def _declare_tool_list_changed(mcp: Any) -> None:
    """Make the server declare ``tools.listChanged`` at initialisation.

    FastMCP builds its initialisation options with default
    ``NotificationOptions`` (every ``*_changed`` False), so a client is
    entitled to ignore the notification. The instance's factory is wrapped so
    every transport declares the capability; nothing else changes.
    """
    from mcp.server.lowlevel.server import NotificationOptions

    low = mcp._mcp_server
    if getattr(low, "_radia_declares_tool_list_changed", False):
        return
    original = low.create_initialization_options

    def create_initialization_options(notification_options=None, experimental_capabilities=None):
        options = notification_options or NotificationOptions()
        options.tools_changed = True
        return original(options, experimental_capabilities)

    low.create_initialization_options = create_initialization_options
    low._radia_declares_tool_list_changed = True


def register_reload_tool(mcp: Any, tool_name: str, module_prefix: str = "radia_mcp") -> None:
    """Register ``tool_name`` on ``mcp``: reload, refresh, notify the client."""
    _declare_tool_list_changed(mcp)

    async def _reload(ctx: Context) -> dict:
        report = reload_and_refresh(mcp, module_prefix)
        notified = False
        try:
            await ctx.session.send_tool_list_changed()
            notified = True
        except Exception as exc:  # pragma: no cover - depends on transport
            report["notify_error"] = f"{type(exc).__name__}: {exc}"
        report["client_notified"] = notified
        report["note"] = (
            "Changed modules were reloaded and their tools re-registered in this "
            "process. Tool schemas on the client refresh only if it honours "
            "notifications/tools/list_changed; otherwise reconnect the server."
        )
        return report

    _reload.__doc__ = (
        "Reload radia_mcp modules whose source changed on disk and re-register "
        "their tools, without restarting this server (editable install). Call "
        "it after editing the package; the report lists reloaded modules, "
        "updated and added tools, and whether the client was told the tool "
        "list changed."
    )
    mcp.add_tool(_reload, name=tool_name)


_PROCESS_START = time.time()
