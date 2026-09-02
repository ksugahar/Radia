"""
server_hardening.py — shared MCP-server hardening infrastructure.

Single source for the MathWorks MATLAB-MCP-server patterns adopted lab-wide
(2026-08-05, first landed in `radia_mcp.cubit`, promoted here for
`radia_mcp.build123d` and future servers):

* **Annotation presets** — every tool is classified through one of four
  fully-specified `ToolAnnotations` presets (MathWorks annotations.go
  discipline: presets only, all four hint fields always set).
* **Gate hiding** — an env var trims the interactive tool surface by
  removing the CI/scenario ``*_gate`` tools (MathWorks exposes 5 tools,
  not 80).
* **All-calls JSONL log** — one choke point wraps ``ToolManager.call_tool``
  and appends ``{ts, tool, args digest, ms, ok, error?}`` per call
  (MathWorks basetool + slog).
* **Error payloads** — uniform ``{status, stage, kind, error, hint?, log?}``
  dicts where ``kind`` tells the LLM audience how to react:
  ``input`` = fix your arguments and retry; ``environment`` = report to
  the user (license/install/hung), do not retry blindly; ``internal`` =
  server bug, do not retry.

Import from here; do not copy these into individual servers (the cubit
server's originals were replaced by these — Discard-the-PoC policy).
"""

from __future__ import annotations

import json
import os
import time

from mcp.types import ToolAnnotations

from . import failure_log as _fl

__all__ = [
    "ANN_READONLY", "ANN_READONLY_WEB", "ANN_WRITES", "ANN_DESTRUCTIVE",
    "classify_tool_annotations", "hide_gate_tools", "install_call_log",
    "error_payload",
    "PROBE_SOLID_CORE_KEYS", "PROBE_FACE_CORE_KEYS",
]

# ---------------------------------------------------------------------------
# Cross-server probe schema contract (cubit <-> build123d <-> external CAD)
# ---------------------------------------------------------------------------
# `cubit_probe(query="entities")` (mesh side) and
# `build123d_probe(query="entities")` (CAD side) MUST both emit these core
# keys per entity, so an agent can compare per-body numbers across the
# STEP -> mesh handoff directly.  Servers may ADD keys (e.g. build123d
# adds `label` -- the named-solid discipline shared with history-based
# CAD systems such as CST); they must not rename or drop these.
# Locked by tests on both servers.
PROBE_SOLID_CORE_KEYS = frozenset(
    {"id", "centroid", "bbox_min", "bbox_max", "extent", "volume"})
PROBE_FACE_CORE_KEYS = frozenset(
    {"id", "center", "bbox_min", "bbox_max", "extent", "area"})


def _ann(read_only: bool, destructive: bool, idempotent: bool,
         open_world: bool) -> ToolAnnotations:
    return ToolAnnotations(
        readOnlyHint=read_only, destructiveHint=destructive,
        idempotentHint=idempotent, openWorldHint=open_world)


ANN_READONLY = _ann(True, False, True, False)
ANN_READONLY_WEB = _ann(True, False, False, True)   # fetches the web
ANN_WRITES = _ann(False, False, False, False)       # writes files; no live-state harm
ANN_DESTRUCTIVE = _ann(False, True, False, False)   # arbitrary code / mutates live state


DEFAULT_READONLY_NAME_HINTS = (
    "_gate", "_docs", "_guide", "_tips", "_reference", "_inventory",
    "_status", "_lookup", "_ask", "_examples", "lint_", "get_",
    "generate_", "netgen_", "_probe", "_diagnose", "_suggest",
    "_failures", "_checkpoints", "_audit", "_doctor", "_usage", "_api",
    "_knowledge", "_manifest", "_handoff", "_crosscheck", "_contract",
    "_discussions", "_inspect", "inspect_",
)


def classify_tool_annotations(mcp, *, destructive=frozenset(),
                              writes=frozenset(), web=frozenset(),
                              readonly_name_hints=DEFAULT_READONLY_NAME_HINTS,
                              ) -> list[str]:
    """Stamp annotation presets onto every registered tool.

    Explicit membership in ``destructive`` / ``writes`` / ``web`` wins;
    everything else falls back to READONLY.  Returns the names that fell
    back WITHOUT a recognizably read-only name shape — surface these in
    the server's --selftest so new tools get a conscious classification.
    """
    unclassified: list[str] = []
    for name, tool in mcp._tool_manager._tools.items():
        if name in destructive:
            tool.annotations = ANN_DESTRUCTIVE
        elif name in writes:
            tool.annotations = ANN_WRITES
        elif name in web:
            tool.annotations = ANN_READONLY_WEB
        elif name.endswith("_reload_code"):
            # Registered by register_status_tool on every server: reloads
            # changed modules and re-registers tools in this process. It
            # writes nothing to disk and touches no live session.
            tool.annotations = ANN_WRITES
        else:
            tool.annotations = ANN_READONLY
            if not any(h in name for h in readonly_name_hints):
                unclassified.append(name)
    return unclassified


def hide_gate_tools(mcp, env_var: str, suffix: str = "_gate") -> int:
    """Remove ``*_gate`` tools when ``env_var`` is set to ``"0"``.

    Default (env unset or any other value) keeps every tool registered so
    the generated docs/TOOLS.md inventory stays stable.  Returns the
    number of tools removed.
    """
    if os.environ.get(env_var, "1") != "0":
        return 0
    removed = 0
    for name in [n for n in list(mcp._tool_manager._tools)
                 if n.endswith(suffix)]:
        mcp._tool_manager.remove_tool(name)
        removed += 1
    return removed


# One-shot size-capped rotation: when the log exceeds this, it is moved
# to <name>.1 (replacing any previous .1) and a fresh file starts -- at
# most 2x the cap on disk, no unbounded growth.
CALL_LOG_ROTATE_BYTES = 5 * 1024 * 1024


def rotate_if_large(log_path, cap_bytes: int = CALL_LOG_ROTATE_BYTES) -> bool:
    """Rotate ``log_path`` to ``<name>.1`` when it exceeds ``cap_bytes``.

    Returns True when a rotation happened.  Best-effort: a locked file
    (concurrent server instance) just skips this round.
    """
    try:
        if log_path.stat().st_size <= cap_bytes:
            return False
        os.replace(log_path, log_path.with_name(log_path.name + ".1"))
        return True
    except OSError:
        return False


def install_call_log(mcp, log_name: str, env_var: str) -> None:
    """Wrap ``ToolManager.call_tool`` with a JSONL all-calls log.

    Every call appends ``{ts, tool, args, ms, ok, error?}`` to
    ``<state_dir>/logs/<log_name>`` (size-capped via
    :func:`rotate_if_large`).  Argument values are recorded as
    truncated reprs (200 chars) so journal bodies / file contents never
    bloat the log.  Set ``env_var=0`` to disable.  Idempotent enough for
    tests: calling again re-wraps against the current ``call_tool``.
    """
    if os.environ.get(env_var, "1") == "0":
        return
    try:
        log_dir = _fl.state_dir() / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / log_name
    except OSError:
        return

    orig_call_tool = mcp._tool_manager.call_tool

    def _digest(arguments: dict) -> dict:
        out = {}
        for k, v in (arguments or {}).items():
            r = repr(v)
            out[k] = r if len(r) <= 200 else r[:200] + f"...({len(r)} ch)"
        return out

    async def logged_call_tool(name, arguments, context=None,
                               convert_result=False):
        t0 = time.time()
        record = {"ts": round(t0, 3), "tool": name,
                  "args": _digest(arguments)}
        try:
            result = await orig_call_tool(
                name, arguments, context=context,
                convert_result=convert_result)
            record["ok"] = True
            return result
        except Exception as exc:
            record["ok"] = False
            record["error"] = f"{type(exc).__name__}: {exc}"
            raise
        finally:
            record["ms"] = round((time.time() - t0) * 1000, 1)
            try:
                rotate_if_large(log_path)
                with open(log_path, "a", encoding="utf-8") as f:
                    f.write(json.dumps(record, ensure_ascii=False,
                                       default=str) + "\n")
            except OSError:
                pass

    mcp._tool_manager.call_tool = logged_call_tool


def error_payload(stage: str, message: str, *, kind: str | None = None,
                  hint: str | None = None,
                  environment_needles: tuple[str, ...] = (),
                  log: str | None = None) -> dict:
    """Uniform error dict: ``{status, stage, kind, error, hint?, log?}``.

    ``kind`` defaults by needle scan: any of ``environment_needles``
    (lower-cased substrings) in the message -> "environment", else
    "input".  Pass ``kind="internal"`` explicitly for server bugs.
    ``log`` (a path/description of where full diagnostics live) is
    attached for environment/internal errors only.
    """
    if kind is None:
        low = message.lower()
        kind = ("environment"
                if any(n in low for n in environment_needles) else "input")
    payload = {"status": "error", "stage": stage, "kind": kind,
               "error": message}
    if hint:
        payload["hint"] = hint
    if log and kind in ("environment", "internal"):
        payload["log"] = log
    return payload
