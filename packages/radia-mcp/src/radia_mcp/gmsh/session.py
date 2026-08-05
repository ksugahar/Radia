"""Persistent gmsh worker session (matlab-mcp-core-server style).

Mirrors the MathWorks MATLAB MCP server architecture on the gmsh side:
a lazily started, PERSISTENT engine process plus a stateful ``evaluate``
(here ``gmsh_exec``), instead of one subprocess per call.  The lab's
radia_mcp.cubit batch daemon uses the same shape (line-delimited
JSON-RPC over stdin/stdout of a worker subprocess).

Design points:

- The worker imports gmsh ONCE and keeps models/options/views alive
  across calls, so a big .msh can be opened once and interrogated many
  times (the one-shot inspect/validate/render tools stay available for
  stateless gating).
- Responses are single lines prefixed with an RS (0x1E) control char;
  anything else the worker prints (gmsh C++ chatter, user ``print``)
  is captured or routed to a log file, so protocol framing survives
  noisy output.
- A reader thread + queue gives every call a hard timeout on Windows
  pipes.  Timeout or worker death KILLS the session and raises --
  fail fast, no silent restart mid-call.  The next call starts a fresh
  worker explicitly (lazy lifecycle, same as MATLAB MCP reconnect).
- The worker stays HEADLESS: no gmsh.fltk in-session (screenshots are
  the one-shot ``gmsh_render``'s job; an FLTK window must never sit
  inside a persistent server-owned process).
"""

from __future__ import annotations

import json
import os
import queue
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path
from typing import Any

from ._gmsh_subprocess import gmsh_available

_RS = "\x1e"  # response-line prefix (record separator)

_WORKER_SCRIPT = r"""
import contextlib
import io
import json
import os
import sys
import time
import traceback

_RS = "\x1e"
_T0 = time.time()

import gmsh
gmsh.initialize(["-noconfig"])
gmsh.option.setNumber("General.Terminal", 0)

_GLOBALS = {"gmsh": gmsh, "__name__": "__gmsh_session__"}
_N_CALLS = 0


def _status():
    try:
        models = list(gmsh.model.list())
    except Exception:
        models = []
    current = ""
    try:
        current = gmsh.model.getCurrent()
    except Exception:
        pass
    return {
        "ok": True,
        "pid": os.getpid(),
        "gmsh_version": gmsh.option.getString("General.Version"),
        "python_version": sys.version.split()[0],
        "uptime_s": round(time.time() - _T0, 3),
        "n_calls": _N_CALLS,
        "models": models,
        "current_model": current,
        "n_views": len(gmsh.view.getTags()),
    }


def _handle(req):
    op = req.get("op")
    if op == "ping":
        return {"ok": True, "pong": True}
    if op == "status":
        return _status()
    if op == "exec":
        code = req.get("code", "")
        buf = io.StringIO()
        _GLOBALS.pop("result", None)
        with contextlib.redirect_stdout(buf):
            exec(compile(code, "<gmsh_exec>", "exec"), _GLOBALS)
        out = {"ok": True, "stdout": buf.getvalue()[-20000:]}
        if "result" in _GLOBALS:
            value = _GLOBALS["result"]
            try:
                json.dumps(value)
                out["result"] = value
            except (TypeError, ValueError):
                out["result"] = repr(value)
                out["result_repr"] = True
        return out
    if op == "shutdown":
        return {"ok": True, "bye": True}
    return {"ok": False, "error": f"unknown op: {op!r}"}


for _line in sys.stdin:
    _line = _line.strip()
    if not _line:
        continue
    try:
        _req = json.loads(_line)
    except json.JSONDecodeError as _exc:
        _resp = {"ok": False, "error": f"bad request line: {_exc}"}
    else:
        _N_CALLS += 1
        try:
            _resp = _handle(_req)
        except SystemExit:
            raise
        except BaseException:
            _resp = {"ok": False, "error": traceback.format_exc()[-4000:]}
    sys.stdout.write(_RS + json.dumps(_resp) + "\n")
    sys.stdout.flush()
    if _resp.get("bye"):
        break

try:
    gmsh.finalize()
except Exception:
    pass
"""

_SESSION_LOCK = threading.Lock()
_SINGLETON: "GmshSession | None" = None

DEFAULT_CALL_TIMEOUT_S = 120.0
_STARTUP_TIMEOUT_S = 60.0


class GmshSessionError(RuntimeError):
    """Raised when the worker dies, times out, or refuses a call."""


class GmshSession:
    """One persistent gmsh worker subprocess (use ``GmshSession.get()``)."""

    def __init__(self):
        if not gmsh_available():
            raise GmshSessionError(
                "gmsh Python package not installed (pip install gmsh)")
        self._log_dir = Path(tempfile.mkdtemp(prefix="radia_mcp_gmsh_session_"))
        self.log_path = self._log_dir / "worker.log"
        self._log_handle = open(self.log_path, "w", encoding="utf-8")
        self._proc = subprocess.Popen(
            [sys.executable, "-u", "-c", _WORKER_SCRIPT],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=self._log_handle,
            text=True,
            encoding="utf-8",
            bufsize=1,
        )
        self._responses: "queue.Queue[str]" = queue.Queue()
        self._reader = threading.Thread(target=self._read_stdout, daemon=True)
        self._reader.start()
        self._call_lock = threading.Lock()
        self.started_at = time.time()
        self.n_calls = 0
        # Fail fast if gmsh cannot even initialize in the worker.
        try:
            info = self.call("ping", timeout_s=_STARTUP_TIMEOUT_S)
        except GmshSessionError:
            self.kill()
            raise
        if not info.get("pong"):
            self.kill()
            raise GmshSessionError(f"gmsh worker failed to start: {info}")

    # ------------------------------------------------------------------
    @classmethod
    def get(cls) -> "GmshSession":
        """Return the live singleton, starting a fresh worker if needed."""
        global _SINGLETON
        with _SESSION_LOCK:
            if _SINGLETON is not None and _SINGLETON.alive():
                return _SINGLETON
            _SINGLETON = cls()
            return _SINGLETON

    @classmethod
    def peek(cls) -> "GmshSession | None":
        """Return the singleton WITHOUT starting one (may be dead)."""
        return _SINGLETON

    @classmethod
    def drop(cls) -> None:
        global _SINGLETON
        with _SESSION_LOCK:
            _SINGLETON = None

    # ------------------------------------------------------------------
    def alive(self) -> bool:
        return self._proc.poll() is None

    def _read_stdout(self) -> None:
        try:
            for line in self._proc.stdout:
                if line.startswith(_RS):
                    self._responses.put(line[1:])
                else:
                    # gmsh C++ chatter or stray writes: keep for debugging
                    self._log_handle.write(line)
                    self._log_handle.flush()
        except (OSError, ValueError):
            pass

    def _log_tail(self, n: int = 1500) -> str:
        try:
            if not self._log_handle.closed:
                self._log_handle.flush()
        except (OSError, ValueError):
            pass
        try:
            return self.log_path.read_text(
                encoding="utf-8", errors="replace")[-n:]
        except OSError:
            return ""

    def call(self, op: str, timeout_s: float = DEFAULT_CALL_TIMEOUT_S,
             **kwargs: Any) -> dict[str, Any]:
        """Send one JSON-RPC request and wait for its response line."""
        with self._call_lock:
            if not self.alive():
                raise GmshSessionError(
                    f"gmsh worker is not running (exit code "
                    f"{self._proc.returncode}); log tail: {self._log_tail()}")
            request = json.dumps({"op": op, **kwargs})
            try:
                self._proc.stdin.write(request + "\n")
                self._proc.stdin.flush()
            except OSError as exc:
                self.kill()
                raise GmshSessionError(
                    f"cannot write to gmsh worker ({exc}); session killed")
            deadline = time.time() + timeout_s
            while True:
                try:
                    raw = self._responses.get(timeout=0.2)
                    break
                except queue.Empty:
                    if not self.alive():
                        tail = self._log_tail()
                        self._close_pipes()
                        raise GmshSessionError(
                            f"gmsh worker died during op={op!r} (exit code "
                            f"{self._proc.returncode}); log tail: {tail}")
                    if time.time() > deadline:
                        self.kill()
                        raise GmshSessionError(
                            f"gmsh worker did not answer op={op!r} within "
                            f"{timeout_s}s; session killed (a hung exec "
                            f"must not block the server). log tail: "
                            f"{self._log_tail()}")
            self.n_calls += 1
            return json.loads(raw)

    def _close_pipes(self) -> None:
        for stream in (self._proc.stdin, self._proc.stdout):
            try:
                if stream is not None:
                    stream.close()
            except OSError:
                pass
        try:
            self._log_handle.close()
        except OSError:
            pass

    def kill(self) -> None:
        """Hard-stop the worker and reap it so alive() flips immediately."""
        try:
            self._proc.kill()
        except OSError:
            pass
        try:
            self._proc.wait(timeout=5)
        except (subprocess.TimeoutExpired, OSError):
            pass
        self._close_pipes()

    def shutdown(self, timeout_s: float = 10.0) -> dict[str, Any]:
        """Graceful shutdown; falls back to kill only after the timeout."""
        result: dict[str, Any] = {"was_alive": self.alive()}
        if self.alive():
            try:
                self.call("shutdown", timeout_s=timeout_s)
            except GmshSessionError as exc:
                result["shutdown_error"] = str(exc)
            try:
                self._proc.wait(timeout=timeout_s)
            except subprocess.TimeoutExpired:
                self.kill()
                result["killed"] = True
        result["exit_code"] = self._proc.poll()
        self._close_pipes()
        GmshSession.drop()
        return result


# ======================================================================
# Tool-facing helpers
# ======================================================================

def session_exec(code: str,
                 timeout_s: float = DEFAULT_CALL_TIMEOUT_S) -> dict[str, Any]:
    """Execute Python code in the persistent gmsh session.

    The session keeps its globals (and gmsh models/views/options) across
    calls; assign to a variable named ``result`` to return a value.
    """
    session = GmshSession.get()
    out = session.call("exec", timeout_s=timeout_s, code=code)
    out["session_pid"] = session._proc.pid
    return out


def session_status() -> dict[str, Any]:
    """Status of the persistent session WITHOUT starting one."""
    session = GmshSession.peek()
    if session is None or not session.alive():
        return {
            "ok": True,
            "running": False,
            "note": ("no persistent gmsh session; the first gmsh_exec "
                     "call starts one"),
        }
    info = session.call("status", timeout_s=30.0)
    info["running"] = True
    info["client_n_calls"] = session.n_calls
    info["worker_log"] = str(session.log_path)
    return info


def session_shutdown() -> dict[str, Any]:
    """Shut the persistent session down (idempotent)."""
    session = GmshSession.peek()
    if session is None:
        return {"ok": True, "running": False, "note": "no session to stop"}
    result = session.shutdown()
    result["ok"] = True
    result["running"] = False
    return result
