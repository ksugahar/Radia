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

import asyncio
import importlib
import os
import sys
import time
from pathlib import Path
from types import ModuleType
from typing import Any

# Imported at module level on purpose: FastMCP evaluates the reload tool's
# annotations in this module's globals, and ``from __future__ import
# annotations`` turns ``Context`` into a string to look up here.
from mcp.server.fastmcp import Context

# Source mtime each module had when the reload tool was registered or when it
# was last reloaded. Comparing two observations of the same source file keeps
# this correct even when a network file server's clock differs from the host.
# A module imported lazily after registration falls back to the process start.
_seen_mtimes: dict[str, int] = {}


def _source_mtime(module: ModuleType) -> int | None:
    path = getattr(module, "__file__", None)
    if not path or not path.endswith(".py"):
        return None
    try:
        return os.stat(path).st_mtime_ns
    except OSError:
        return None


def _reload_order(name: str) -> tuple[int, int, str]:
    parts = name.split(".")
    # Packages (``__init__``) and ``tools`` modules import everything else in
    # their tree, so they go last; deeper modules go first.
    last = parts[-1]
    tail_rank = 2 if last in {"tools", "server"} else 0
    return (-len(parts), tail_rank, name)


_RELOAD_METADATA = frozenset(
    {
        "__name__",
        "__loader__",
        "__package__",
        "__spec__",
        "__path__",
        "__file__",
        "__cached__",
        "__builtins__",
    }
)


def _clean_module_namespace(module: ModuleType) -> None:
    """Clear stale definitions while retaining import-system metadata."""
    namespace = module.__dict__
    metadata = {key: namespace[key] for key in _RELOAD_METADATA if key in namespace}
    namespace.clear()
    namespace.update(metadata)


def _restore_modules(snapshots: dict[str, dict[str, Any]]) -> None:
    for name, snapshot in snapshots.items():
        module = sys.modules.get(name)
        if module is None:
            continue
        module.__dict__.clear()
        module.__dict__.update(snapshot)


def _modules_owning_object(obj: Any, prefix: str) -> set[str]:
    """Find modules that expose ``obj`` as a global (normally server.py)."""
    owners: set[str] = set()
    for name, module in list(sys.modules.items()):
        if module is None or not (name == prefix or name.startswith(prefix + ".")):
            continue
        try:
            if any(value is obj for value in module.__dict__.values()):
                owners.add(name)
        except RuntimeError:
            continue
    return owners


def _prime_mtimes(prefix: str) -> None:
    """Record the startup baseline for modules already imported by a server."""
    for name, module in list(sys.modules.items()):
        if module is None or not (name == prefix or name.startswith(prefix + ".")):
            continue
        mtime = _source_mtime(module)
        if mtime is not None:
            _seen_mtimes.setdefault(name, mtime)


def reload_changed_modules(
    prefix: str = "radia_mcp",
    *,
    protected_modules: set[str] | None = None,
) -> dict[str, Any]:
    """Reload every imported module under ``prefix`` whose file changed.

    The first call reloads what changed since the process started; later
    calls reload what changed since the previous call. Returns the module
    names reloaded and any reload errors (a module that fails to import is
    reported and every touched module is restored to its previous namespace).

    ``protected_modules`` are reported as requiring a process restart. This is
    used for a server module that owns the live FastMCP instance: re-executing
    that module would create a second, orphaned server instead of updating the
    connected one. The reload implementation module itself is protected for the
    same reason.
    """
    protected = set(protected_modules or ())
    protected.add(__name__)
    changed: dict[str, int] = {}
    restart_required: list[str] = []
    for name, module in list(sys.modules.items()):
        if module is None or not (name == prefix or name.startswith(prefix + ".")):
            continue
        mtime = _source_mtime(module)
        if mtime is None:
            continue
        seen = _seen_mtimes.get(name)
        if (seen is None and mtime > _PROCESS_START_NS) or (seen is not None and mtime != seen):
            if name in protected:
                restart_required.append(name)
            else:
                changed[name] = mtime
    order = sorted(changed, key=_reload_order)

    importlib.invalidate_caches()
    errors: dict[str, str] = {}
    compiled: dict[str, Any] = {}
    for name in order:
        module = sys.modules.get(name)
        path = getattr(module, "__file__", None) if module is not None else None
        if not path:
            continue
        try:
            source = Path(path).read_bytes()
            compiled[name] = compile(source, path, "exec")
        except Exception as exc:
            errors[name] = f"{type(exc).__name__}: {exc}"

    if errors:
        return {
            "reloaded": [],
            "errors": errors,
            "rolled_back": False,
            "restart_required": sorted(restart_required),
        }

    snapshots = {
        name: dict(sys.modules[name].__dict__)
        for name in order
        if name in sys.modules
    }
    try:
        # A clean namespace makes deleted definitions disappear. Two passes
        # refresh aliases imported from siblings that were reloaded later in
        # the first pass. Executing the already-compiled source also avoids the
        # same-second/same-size stale-pyc trap of importlib.reload().
        for _pass in range(2):
            for name in order:
                module = sys.modules.get(name)
                if module is None or name not in compiled:
                    continue
                _clean_module_namespace(module)
                exec(compiled[name], module.__dict__)
    except Exception as exc:  # pragma: no cover - exact import failure varies
        errors[name] = f"{type(exc).__name__}: {exc}"
        _restore_modules(snapshots)
        return {
            "reloaded": [],
            "errors": errors,
            "rolled_back": True,
            "restart_required": sorted(restart_required),
        }

    for name, mtime in changed.items():
        _seen_mtimes[name] = mtime
    return {
        "reloaded": list(order),
        "errors": {},
        "rolled_back": False,
        "restart_required": sorted(restart_required),
    }


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


def refresh_tools(
    mcp: Any,
    module_prefix: str = "radia_mcp",
    *,
    reloaded_modules: set[str] | None = None,
) -> dict[str, Any]:
    """Re-register stale tools and pick up callables added since start-up."""
    manager = mcp._tool_manager
    updated: list[str] = []
    added: list[str] = []
    removed: list[str] = []
    skipped: list[str] = []
    names_by_module: dict[str, list[str]] = {}
    selected_modules = set(reloaded_modules or ())
    registry = getattr(manager, "_tools", None)
    if not isinstance(registry, dict):
        raise RuntimeError("unsupported FastMCP ToolManager registry")
    registry_snapshot = dict(registry)

    try:
        for tool in list(manager.list_tools()):
            fn = getattr(tool, "fn", None)
            module_name = getattr(fn, "__module__", "") or ""
            if (
                not module_name.startswith(module_prefix)
                or module_name not in selected_modules
            ):
                continue
            module = sys.modules.get(module_name)
            if module is None:
                skipped.append(tool.name)
                continue
            names_by_module.setdefault(module_name, []).append(tool.name)
            # Closures such as the shared status/reload tools cannot be looked
            # up as module globals. Their global namespace is refreshed in
            # place when the owning module reloads, so leave them registered.
            if "<locals>" in getattr(fn, "__qualname__", ""):
                skipped.append(tool.name)
                continue
            fresh = getattr(module, fn.__name__, None)
            if fresh is None or not callable(fresh):
                mcp.remove_tool(tool.name)
                removed.append(tool.name)
                continue
            if fresh is fn:
                continue
            mcp.remove_tool(tool.name)
            mcp.add_tool(
                fresh,
                name=tool.name,
                title=getattr(tool, "title", None),
                annotations=getattr(tool, "annotations", None),
                icons=getattr(tool, "icons", None),
                meta=getattr(tool, "meta", None),
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
                    # The server-specific annotation/contract pass does not run
                    # again until reconnect. A newly discovered callable is
                    # therefore deliberately conservative for this session.
                    from mcp.types import ToolAnnotations

                    mcp.add_tool(
                        candidate,
                        annotations=ToolAnnotations(
                            readOnlyHint=False,
                            destructiveHint=True,
                            idempotentHint=False,
                            openWorldHint=True,
                        ),
                    )
                    registered.add(attr)
                    added.append(attr)
    except Exception:
        registry.clear()
        registry.update(registry_snapshot)
        raise
    return {
        "updated": updated,
        "added": added,
        "removed": removed,
        "skipped": skipped,
        "added_tools_need_reconnect_for_server_policy": bool(added),
    }


def reload_and_refresh(mcp: Any, module_prefix: str = "radia_mcp") -> dict[str, Any]:
    module_snapshots = {
        name: dict(module.__dict__)
        for name, module in list(sys.modules.items())
        if module is not None
        and (name == module_prefix or name.startswith(module_prefix + "."))
    }
    seen_snapshot = dict(_seen_mtimes)
    protected = _modules_owning_object(mcp, module_prefix)
    report = reload_changed_modules(module_prefix, protected_modules=protected)
    if report["errors"]:
        report.update(
            {
                "updated": [],
                "added": [],
                "removed": [],
                "skipped": [],
                "added_tools_need_reconnect_for_server_policy": False,
            }
        )
        return report
    try:
        report.update(
            refresh_tools(
                mcp,
                module_prefix,
                reloaded_modules=set(report["reloaded"]),
            )
        )
    except Exception as exc:  # restore code and tools as one transaction
        _restore_modules(
            {
                name: module_snapshots[name]
                for name in report["reloaded"]
                if name in module_snapshots
            }
        )
        _seen_mtimes.clear()
        _seen_mtimes.update(seen_snapshot)
        report.update(
            {
                "reloaded": [],
                "errors": {"tool_registry": f"{type(exc).__name__}: {exc}"},
                "rolled_back": True,
                "updated": [],
                "added": [],
                "removed": [],
                "skipped": [],
                "added_tools_need_reconnect_for_server_policy": False,
            }
        )
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
    _prime_mtimes(module_prefix)
    _declare_tool_list_changed(mcp)
    reload_lock = asyncio.Lock()

    async def _reload(ctx: Context) -> dict:
        async with reload_lock:
            report = reload_and_refresh(mcp, module_prefix)
        notified = False
        try:
            await ctx.session.send_tool_list_changed()
            notified = True
        except Exception as exc:  # pragma: no cover - depends on transport
            report["notify_error"] = f"{type(exc).__name__}: {exc}"
        report["client_notified"] = notified
        report["note"] = (
            "Changed implementation modules were reloaded transactionally and "
            "their tools refreshed in this process. A changed server module, "
            "the reload implementation itself, or a newly added tool that needs "
            "the server-specific policy pass still requires reconnecting. Tool "
            "schemas refresh only if the client honours "
            "notifications/tools/list_changed."
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


_PROCESS_START_NS = time.time_ns()
