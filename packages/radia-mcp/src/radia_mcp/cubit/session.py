"""
session.py — Python-3.12-side client for the Cubit Viewer daemon.

This module is version-agnostic: stdlib only, works on any Python
>= 3.10.

Two transport modes share the same public API (call/ping/shutdown):

    mode="gui"    (DEFAULT, Plan A — 2026-04-19 file-drop bootstrap)
        ─ Launches `coreform_cubit.exe -nojournal bootstrap.py`
        ─ Bootstrap installs a QTimer inside Cubit's Qt event loop
        ─ Client ↔ Bootstrap communicate via atomically-renamed JSON
          files in a per-session drop directory (no sockets)
        ─ Cubit GUI window is live & user-interactive throughout

    mode="batch"  (legacy, CI/scripting)
        ─ Launches Cubit's bundled Python 3.10 with daemon.py
        ─ Client ↔ Daemon communicate via line-delimited JSON-RPC on
          stdin/stdout of the subprocess
        ─ No GUI, no user interaction; pure headless

Both modes preserve a persistent Cubit session across many client calls
(~2 s cold start amortized) — launching per-call would be prohibitively
slow. A process-wide singleton (`CubitSession.get()`) holds the handle.

Protocol version:
    v1 — stdio JSON-RPC (daemon.py)
    v2 — file drop JSON-RPC (bootstrap.py, this module default)

The ready message echoes `protocol_version` for mutual compatibility.
"""

from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from pathlib import Path
from typing import Any


PROTOCOL_VERSION = 2  # default (gui mode)

# Startup timeouts (step 2 of 2026-04-21 speed fix).
# License checkout: RLM server round-trip can take 30+ s on cold start;
# the warmup helper tries pre-warm via rlm_activate.exe with a shorter
# budget (30 s) so the main wait only covers post-warmup Cubit startup.
# Cubit ready: once the license is warm, Qt + bootstrap init finish in
# 2 – 3 s on this lab (measured on LAB 2026-04-21).
LICENSE_WARMUP_TIMEOUT_S = 30.0
CUBIT_READY_TIMEOUT_S = 15.0
# Attach-time responsiveness probe: the bootstrap QTimer polls every
# 200 ms, so a healthy idle daemon answers a ping in well under 1 s.
ATTACH_PING_TIMEOUT_S = 3.0

_SESSION_LOCK = threading.Lock()
_SINGLETON: "CubitSession | None" = None


def _user_daemon_dir() -> Path:
    """Return the per-user stable path for the Cubit daemon's drop-dir.

    Phase 1 of "MCP daemon survives VSCode restart" (2026-04-22).
    Uses ``%LOCALAPPDATA%`` so each Windows user has their own
    daemon, license cache, and drop-dir location.  Stable across
    Python process restarts, so a fresh MCP server can discover
    and attach to an already-running Cubit.
    """
    base = os.environ.get("LOCALAPPDATA") or os.environ.get("TEMP", ".")
    d = Path(base) / "radia-mcp" / "cubit-session"
    return d


def _is_pid_alive(pid: int) -> bool:
    """Check whether the given PID refers to a running process.

    Works on Windows via ``OpenProcess`` with low query rights; on
    POSIX via ``os.kill(pid, 0)``.
    """
    if pid <= 0:
        return False
    try:
        if sys.platform == "win32":
            import ctypes
            PROCESS_QUERY_LIMITED = 0x1000
            STILL_ACTIVE = 259
            h = ctypes.windll.kernel32.OpenProcess(PROCESS_QUERY_LIMITED, False, pid)
            if not h:
                return False
            try:
                code = ctypes.c_ulong()
                ok = ctypes.windll.kernel32.GetExitCodeProcess(h, ctypes.byref(code))
                return bool(ok) and code.value == STILL_ACTIVE
            finally:
                ctypes.windll.kernel32.CloseHandle(h)
        else:
            os.kill(pid, 0)
            return True
    except (OSError, ProcessLookupError):
        return False


def _assign_kill_on_close_job(pid: int):
    """Put ``pid`` into a Windows Job Object with KILL_ON_JOB_CLOSE.

    Used for session_mode="new" ONLY: the private per-process daemon
    must die with its client (a crashed client otherwise leaves a
    license-holding orphan -- observed 2026-08-05).  The returned job
    HANDLE must stay referenced for the client's lifetime; when this
    process exits (normally or not), the OS closes the handle and kills
    the job's processes.  Shared-daemon modes never call this -- their
    persistence across client restarts is the feature.

    Returns the handle. Failure raises ``OSError``: mode ``new`` must
    not continue after losing the cleanup guarantee it promises.
    """
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE = 0x2000
    JobObjectExtendedLimitInformation = 9
    PROCESS_SET_QUOTA = 0x0100
    PROCESS_TERMINATE = 0x0001

    class IO_COUNTERS(ctypes.Structure):
        _fields_ = [(n, ctypes.c_uint64) for n in (
            "ReadOperationCount", "WriteOperationCount",
            "OtherOperationCount", "ReadTransferCount",
            "WriteTransferCount", "OtherTransferCount")]

    class JOBOBJECT_BASIC_LIMIT_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("PerProcessUserTimeLimit", ctypes.c_int64),
            ("PerJobUserTimeLimit", ctypes.c_int64),
            ("LimitFlags", wintypes.DWORD),
            ("MinimumWorkingSetSize", ctypes.c_size_t),
            ("MaximumWorkingSetSize", ctypes.c_size_t),
            ("ActiveProcessLimit", wintypes.DWORD),
            ("Affinity", ctypes.c_size_t),
            ("PriorityClass", wintypes.DWORD),
            ("SchedulingClass", wintypes.DWORD),
        ]

    class JOBOBJECT_EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
        _fields_ = [
            ("BasicLimitInformation", JOBOBJECT_BASIC_LIMIT_INFORMATION),
            ("IoInfo", IO_COUNTERS),
            ("ProcessMemoryLimit", ctypes.c_size_t),
            ("JobMemoryLimit", ctypes.c_size_t),
            ("PeakProcessMemoryUsed", ctypes.c_size_t),
            ("PeakJobMemoryUsed", ctypes.c_size_t),
        ]

    kernel32.CreateJobObjectW.argtypes = [wintypes.LPVOID, wintypes.LPCWSTR]
    kernel32.CreateJobObjectW.restype = wintypes.HANDLE
    kernel32.SetInformationJobObject.argtypes = [
        wintypes.HANDLE, ctypes.c_int, wintypes.LPVOID, wintypes.DWORD]
    kernel32.SetInformationJobObject.restype = wintypes.BOOL
    kernel32.OpenProcess.argtypes = [
        wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.AssignProcessToJobObject.argtypes = [
        wintypes.HANDLE, wintypes.HANDLE]
    kernel32.AssignProcessToJobObject.restype = wintypes.BOOL
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel32.CloseHandle.restype = wintypes.BOOL

    job = None
    hproc = None
    try:
        job = kernel32.CreateJobObjectW(None, None)
        if not job:
            raise ctypes.WinError(ctypes.get_last_error())
        info = JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
        info.BasicLimitInformation.LimitFlags = \
            JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        if not kernel32.SetInformationJobObject(
                job, JobObjectExtendedLimitInformation,
                ctypes.byref(info), ctypes.sizeof(info)):
            raise ctypes.WinError(ctypes.get_last_error())
        hproc = kernel32.OpenProcess(
            PROCESS_SET_QUOTA | PROCESS_TERMINATE, False, pid)
        if not hproc:
            raise ctypes.WinError(ctypes.get_last_error())
        if not kernel32.AssignProcessToJobObject(job, hproc):
            raise ctypes.WinError(ctypes.get_last_error())
        kernel32.CloseHandle(hproc)
        hproc = None
        return job
    except Exception:
        if hproc:
            kernel32.CloseHandle(hproc)
        if job:
            kernel32.CloseHandle(job)
        raise


def _close_windows_handle(handle) -> None:
    """Close a checked Windows HANDLE or raise ``OSError``."""
    import ctypes
    from ctypes import wintypes

    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel32.CloseHandle.restype = wintypes.BOOL
    if not kernel32.CloseHandle(handle):
        raise ctypes.WinError(ctypes.get_last_error())


def _remove_tree_with_retry(path: Path, timeout_s: float = 2.0) -> bool:
    """Remove a process drop directory after Windows releases log handles."""
    deadline = time.monotonic() + max(0.0, timeout_s)
    while True:
        shutil.rmtree(path, ignore_errors=True)
        if not path.exists():
            return True
        if time.monotonic() >= deadline:
            return False
        time.sleep(0.05)


def _try_attach_existing_daemon(drop_dir: Path) -> dict | None:
    """Return the daemon's ready_info if a live daemon is found, else None.

    Attachment criteria (all must hold):
      * ``pid.lock`` exists under ``drop_dir``
      * the PID inside is a live process
      * ``ready`` marker exists (daemon finished bootstrap, JSON-RPC
        ready to accept calls)

    On attach success, the caller reuses ``drop_dir`` as-is.  On
    failure, stale artifacts are removed so the next spawn starts
    clean.
    """
    pid_file = drop_dir / "pid.lock"
    ready_file = drop_dir / "ready"
    if not pid_file.exists() or not ready_file.exists():
        return None
    try:
        pid = int(pid_file.read_text(encoding="utf-8").strip())
    except (OSError, ValueError):
        return None
    if not _is_pid_alive(pid):
        # Stale lock — clean up for the next spawn
        for f in (pid_file, ready_file):
            try:
                f.unlink()
            except OSError:
                pass
        return None
    try:
        info = json.loads(ready_file.read_text(encoding="utf-8"))
        info["attached"] = True  # flag so callers can log "reused"
        info["pid"] = pid
        return info
    except (OSError, json.JSONDecodeError):
        return None

# Module-level overrides (OCP-inspired: cf. set_port / get_port).
_OVERRIDE_BIN_DIR: Path | None = None


def set_cubit_bin_dir(path: str | Path | None) -> None:
    """Override the Cubit `bin/` directory lookup for subsequent sessions."""
    global _OVERRIDE_BIN_DIR
    _OVERRIDE_BIN_DIR = Path(path) if path is not None else None


def get_cubit_bin_dir() -> Path | None:
    """Return the currently-resolved Cubit `bin/` directory, or None."""
    if _OVERRIDE_BIN_DIR is not None:
        return _OVERRIDE_BIN_DIR
    return find_cubit_install()


def _session_info_path() -> Path:
    return Path.home() / ".cubit_viewer" / "session.json"


# ---------------------------------------------------------------------------
# Cubit install auto-discovery
# ---------------------------------------------------------------------------

def find_cubit_install(explicit: str | None = None) -> Path | None:
    """Locate a Coreform Cubit install's `bin/` directory."""
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit))
    env_bin = os.environ.get("CUBIT_BIN_DIR")
    if env_bin:
        candidates.append(Path(env_bin))

    if sys.platform == "win32":
        try:
            import winreg

            with winreg.OpenKey(
                winreg.HKEY_LOCAL_MACHINE, r"SOFTWARE\Coreform"
            ) as hk:
                i = 0
                while True:
                    try:
                        subkey = winreg.EnumKey(hk, i)
                    except OSError:
                        break
                    try:
                        with winreg.OpenKey(hk, subkey) as sk:
                            install_dir, _ = winreg.QueryValueEx(
                                sk, "InstallDir")
                            candidates.append(Path(install_dir) / "bin")
                    except OSError:
                        pass
                    i += 1
        except OSError:
            pass

    for pattern in (
        r"C:/Program Files/Coreform Cubit */bin",
        r"C:/Program Files (x86)/Coreform Cubit */bin",
        "/opt/Coreform-Cubit-*/bin",
        "/opt/Coreform/Cubit*/bin",
    ):
        for path in sorted(glob.glob(pattern), reverse=True):
            candidates.append(Path(path))

    for c in candidates:
        if (c / "cubit.py").exists() or (c / "_cubit3.pyd").exists() \
                or (c / "_cubit3.so").exists():
            return c.resolve()
    return None


def _cubit_python_exe(bin_dir: Path) -> Path:
    """Locate Cubit's bundled Python interpreter under `bin_dir`."""
    candidates = [
        bin_dir / "python3" / "python.exe",
        bin_dir / "python3" / "python3.exe",
        bin_dir / "python3" / "bin" / "python3",
        bin_dir / "python3" / "bin" / "python",
    ]
    for c in candidates:
        if c.exists():
            return c.resolve()
    raise FileNotFoundError(
        f"Cubit bundled Python not found under {bin_dir}/python3/. "
        "Expected python.exe (Windows) or bin/python3 (Linux)."
    )


def _cubit_gui_exe(bin_dir: Path) -> Path:
    """Locate the Cubit GUI launcher (`coreform_cubit.exe` / `coreform_cubit`)."""
    candidates = [
        bin_dir / "coreform_cubit.exe",
        bin_dir / "coreform_cubit",
    ]
    for c in candidates:
        if c.exists():
            return c.resolve()
    raise FileNotFoundError(
        f"coreform_cubit launcher not found under {bin_dir}.")


# ---------------------------------------------------------------------------
# Session (singleton per mcp-server process)
# ---------------------------------------------------------------------------

class CubitSessionError(Exception):
    """Raised when the Cubit daemon cannot be launched or has died."""


class CubitSession:
    """Manages a long-lived Cubit session via one of two transports.

    Thread-safe via an internal lock around each call.

    Usage:
        session = CubitSession.get()   # mode="gui" by default
        r = session.call("cmd", ["create brick x 10"])
        if r["ok"]:
            ...
    """

    def __init__(self,
                 cubit_bin_dir: Path | None = None,
                 mode: str = "gui"):
        self._bin_dir = cubit_bin_dir or find_cubit_install()
        if self._bin_dir is None:
            raise CubitSessionError(
                "Could not locate Coreform Cubit install. Set CUBIT_BIN_DIR "
                "or install to a standard location."
            )
        if mode not in ("gui", "batch"):
            raise ValueError(f"mode must be 'gui' or 'batch', got {mode!r}")
        self._mode = mode
        self._proc: subprocess.Popen | None = None
        self._next_id = 1
        self._client_id = f"{os.getpid():08x}-{uuid.uuid4().hex[:12]}"
        self._lock = threading.Lock()
        self._ready_info: dict | None = None

        # Ownership tag (MathWorks pattern): True when THIS process
        # spawned the Cubit runner; False when we merely attached to a
        # daemon another process started.  Recovery paths must never
        # kill a live session we do not own.
        self._owned = False
        # Windows Job Object handle for session_mode="new" (kill the
        # private daemon when this client dies); None otherwise.
        self._job_handle = None

        # gui-mode (file-drop) state
        self._drop_dir: Path | None = None
        self._outbox: Path | None = None

        # last license pre-warm result (None until first start)
        self._last_license_warmup: dict = {}

        # Per-MCP-process command history for `cubit_session_journal`
        # (the lab's .jou-first reproducibility policy: every live
        # session must be exportable as a portable journal).  Bounded;
        # records {ts, line, ok} for every op=="cmd" line sent.
        self._command_history: list[dict] = []
        self._command_history_max = 20000

        # batch-mode (stdio) stderr retention
        self._stderr_tail: list[bytes] = []
        self._stderr_tail_max = 200

    def _close_private_job(self) -> None:
        """Release the mode-new Job Object, killing any surviving child."""
        handle, self._job_handle = getattr(self, "_job_handle", None), None
        if handle is None or sys.platform != "win32":
            return
        try:
            _close_windows_handle(handle)
        except OSError:
            pass

    def __del__(self):
        # A mode-new session is private to this object. Losing the object
        # must not lose the only handle that enforces kill-on-close.
        try:
            self._close_private_job()
        except Exception:
            pass

    # ---- lifecycle ----

    def ensure_started(self) -> dict:
        """Spawn (or attach to) the Cubit session if not already running.

        Three live-state paths:
          1. Our own child is alive → reuse (``self._proc.poll() is None``)
          2. GUI mode: an external daemon is alive and we've already
             attached to it this Python process → reuse
             (``self._drop_dir`` set + ``pid.lock`` matches a live PID)
          3. Nothing alive → call _start_gui_bootstrap / _start_stdio_daemon
        """
        if self._proc is not None and self._proc.poll() is None:
            return self._ready_info or {}
        if self._proc is not None:
            self._close_private_job()
            self._proc = None
        if self._mode == "gui" and self._drop_dir is not None and self._ready_info is not None:
            # Phase-1 attached mode: verify the daemon PID is still alive
            pid_file = self._drop_dir / "pid.lock"
            if pid_file.exists():
                try:
                    pid = int(pid_file.read_text(encoding="utf-8").strip())
                    if _is_pid_alive(pid):
                        return self._ready_info
                except (OSError, ValueError):
                    pass
            # daemon died — fall through to re-attach or re-spawn
            self._drop_dir = None
            self._outbox = None
            self._ready_info = None
        if self._mode == "gui":
            return self._start_gui_bootstrap()
        return self._start_stdio_daemon()

    def shutdown(self, timeout_s: float = 3.0) -> dict:
        """Politely ask Cubit to exit, then force-kill after timeout.

        This is the EXPLICIT stop path (`cubit_session_shutdown` tool),
        so it may also kill an attached daemon another process started --
        but it reports what it did, so the caller can tell the user which
        process was stopped (``stopped``: "owned-child" |
        "attached-daemon" | "none"; ``pid`` when known).
        """
        report: dict = {"stopped": "none", "pid": None,
                        "owned": self._owned}
        # Best-effort polite shutdown via live transport
        try:
            self._lock_free_shutdown_op(timeout_s=timeout_s)
        except Exception:
            pass
        # Kill our own child process
        if self._proc is not None:
            report["pid"] = self._proc.pid
            try:
                self._proc.wait(timeout=timeout_s)
            except Exception:
                pass
            if self._proc.poll() is None:
                try:
                    self._proc.kill()
                except Exception:
                    pass
            # Reap after kill: releases the child's handles (so the
            # drop-dir rmtree below can succeed on Windows) and avoids
            # the Popen.__del__ "still running" ResourceWarning.
            try:
                self._proc.wait(timeout=10.0)
            except Exception:
                pass
            report["stopped"] = "owned-child"
        self._close_private_job()
        # Kill an attached daemon (not our child) via pid.lock
        if self._drop_dir is not None:
            pid_file = self._drop_dir / "pid.lock"
            if pid_file.exists():
                try:
                    pid = int(pid_file.read_text(encoding="utf-8").strip())
                    if pid > 0 and _is_pid_alive(pid):
                        if report["stopped"] == "none":
                            report["stopped"] = "attached-daemon"
                            report["pid"] = pid
                            report["note"] = (
                                "stopped a daemon this process did not "
                                "start (explicit shutdown request)")
                        if sys.platform == "win32":
                            import ctypes
                            h = ctypes.windll.kernel32.OpenProcess(
                                0x0001, False, pid)  # PROCESS_TERMINATE
                            if h:
                                ctypes.windll.kernel32.TerminateProcess(h, 1)
                                ctypes.windll.kernel32.CloseHandle(h)
                        else:
                            os.kill(pid, 9)
                except (OSError, ValueError):
                    pass
                try:
                    pid_file.unlink()
                except OSError:
                    pass
            # Also remove ready marker so the next ensure_started spawns fresh
            ready = self._drop_dir / "ready"
            if ready.exists():
                try: ready.unlink()
                except OSError: pass
            # Session-mode "new" uses a private per-process drop dir --
            # remove it entirely (the shared dir is named exactly
            # "cubit-session" and is preserved).  Windows may hold the
            # child's log-file handles for a moment after the kill, so
            # retry briefly instead of leaving debris.
            if self._owned and self._drop_dir.name != "cubit-session":
                if not _remove_tree_with_retry(self._drop_dir):
                    report["cleanup_warning"] = (
                        f"private drop directory remains: {self._drop_dir}")
        self._drop_dir = None
        self._outbox = None
        self._proc = None
        self._ready_info = None
        self._owned = False
        return report

    def is_alive(self) -> bool:
        # Our own child
        if self._proc is not None and self._proc.poll() is None:
            return True
        # Phase-1 attached daemon: check pid.lock against a live PID
        if self._mode == "gui" and self._drop_dir is not None:
            pid_file = self._drop_dir / "pid.lock"
            if pid_file.exists():
                try:
                    pid = int(pid_file.read_text(encoding="utf-8").strip())
                    return _is_pid_alive(pid)
                except (OSError, ValueError):
                    pass
        return False

    def ping(self, timeout_s: float = 5.0) -> bool:
        """Heartbeat probe — returns True if Cubit responds."""
        try:
            r = self.call("ping", timeout_s=timeout_s)
            return bool(r.get("ok") and r.get("result") == "pong")
        except Exception:
            return False

    # ---- public RPC entrypoint ----

    def call(self, op: str, args: list | None = None,
             timeout_s: float = 60.0, _recover: bool = True) -> dict:
        """Send one JSON-RPC request, return the response dict.

        Response shape:
            {"id": int, "ok": bool, "result": ..., "error": str?}

        Error recovery (gui mode, 2026-04-19): if Cubit died mid-call,
        the session is transparently reset and the call retried once.
        `result` will then carry `_recovered=True` so callers know state
        was lost (no geometry history; the user may need to re-import
        or `cubit_restore` from a checkpoint). On second failure, raises
        CubitSessionError.

        Pass `_recover=False` to disable the one-shot retry (used by the
        recovery path itself to avoid infinite loops).
        """
        with self._lock:
            try:
                self.ensure_started()
                req_id = self._next_id
                self._next_id += 1
                req = {
                    "id": req_id,
                    "op": op,
                    "args": args or [],
                    "protocol_version": PROTOCOL_VERSION if self._mode == "gui" else 1,
                }
                if self._mode == "gui":
                    resp = self._call_via_filedrop(req, timeout_s=timeout_s)
                else:
                    resp = self._call_via_stdio(req, timeout_s=timeout_s)
                if op == "cmd":
                    self._record_cmd_history(resp)
                return resp
            except CubitSessionError as e:
                if not _recover:
                    raise
                # Auto-restart: Cubit died or misbehaved. Reset and try once.
                self._force_reset()
        # Outside the lock (avoid deadlock) — retry once with recovery off.
        resp = self.call(op, args=args, timeout_s=timeout_s, _recover=False)
        if isinstance(resp, dict):
            resp.setdefault("_recovered", True)
        return resp

    def _record_cmd_history(self, resp) -> None:
        """Append per-line results of an op=="cmd" response to the
        session journal history (bounded)."""
        if not isinstance(resp, dict):
            return
        per_line = resp.get("result")
        if not isinstance(per_line, list):
            return
        now = time.time()
        for step in per_line:
            if not isinstance(step, dict) or "line" not in step:
                continue
            self._command_history.append({
                "ts": now,
                "line": str(step.get("line")),
                "ok": bool(step.get("ok")),
            })
        if len(self._command_history) > self._command_history_max:
            del self._command_history[:len(self._command_history)
                                      - self._command_history_max]

    def _force_reset(self) -> None:
        """Tear down the current session without waiting for graceful exit.

        Used by the recovery path in `call()`. Does NOT re-enter the lock
        (caller is holding it).

        Ownership rule (MathWorks pattern, 2026-08-05): a recovery path
        may kill only a session THIS process spawned.  An attached daemon
        that is still alive belongs to whoever started it (another VSCode
        window's live geometry, potentially) -- we DETACH from it and
        leave its pid.lock/ready markers intact for its other clients.
        Killing an attached-but-hung daemon requires the explicit
        `cubit_session_shutdown` tool, not an implicit retry path.
        """
        proc = self._proc
        if proc is not None and proc.poll() is None:
            try:
                proc.kill()
            except Exception:
                pass
        if self._drop_dir is not None:
            pid_file = self._drop_dir / "pid.lock"
            attached_alive = False
            if not self._owned and pid_file.exists():
                try:
                    pid = int(pid_file.read_text(encoding="utf-8").strip())
                    attached_alive = pid > 0 and _is_pid_alive(pid)
                except (OSError, ValueError):
                    attached_alive = False
            if attached_alive:
                # Detach only: keep pid.lock + ready for the daemon's
                # other clients; the retry will re-attach (or surface
                # the error if the daemon is truly hung).
                pass
            else:
                # Dead daemon (or our own killed child): clean the
                # markers so the next ensure_started spawns fresh.
                if pid_file.exists():
                    try: pid_file.unlink()
                    except OSError: pass
                ready = self._drop_dir / "ready"
                if ready.exists():
                    try: ready.unlink()
                    except OSError: pass
        # Per-process state: cleared.  The per-user drop-dir on disk
        # persists (it's _user_daemon_dir()) so the next spawn reuses it.
        self._drop_dir = None
        self._outbox = None
        self._proc = None
        self._ready_info = None
        self._owned = False

    # ---- mode: gui (file drop) ----

    def _start_gui_bootstrap(self) -> dict:
        """Launch (or attach to) Cubit + bootstrap, wait for ready.

        Phase-1 daemon persistence (2026-04-22): drop-dir is now per-user
        stable (``%LOCALAPPDATA%\\radia-mcp\\cubit-session``), not a
        per-process tempdir.  Before spawning a new Cubit, check if a
        previously-launched daemon is still alive and attach to it —
        this lets a fresh VSCode / MCP-server process reuse the Cubit
        that an earlier session left running (license stays warm,
        imported STEPs survive, sphere in GUI stays as-is).

        Fallback to spawn when no daemon found, with four layers:
        1. Pre-warm Learn license via rlm_activate cache check.
        2. Spawn Cubit with ``DETACHED_PROCESS`` + ``CREATE_NEW_PROCESS_GROUP``
           so the subprocess survives this MCP-server process exit.
           stdout/stderr → DEVNULL to avoid pipe breakage on parent exit.
        3. Wait for ready marker, split timeouts (15 s / 60 s by whether
           license warmup ran).
        4. Record our PID into ``pid.lock`` so the next process can
           find + attach.
        """
        # --- Session mode (MathWorks matlab-session-mode triad) ---------
        #   auto     (default) attach to a live shared daemon, else spawn
        #   new      always spawn a FRESH daemon in a private drop-dir
        #            (hermetic runs -- never reuses / clobbers the shared
        #            per-user daemon)
        #   existing attach only; fail loud when no live shared daemon
        session_mode = os.environ.get(
            "RADIA_CUBIT_SESSION_MODE", "auto").strip().lower()
        if session_mode not in ("auto", "new", "existing"):
            raise CubitSessionError(
                f"Invalid RADIA_CUBIT_SESSION_MODE={session_mode!r}: "
                "expected auto | new | existing")

        # --- Drop-dir ---------------------------------------------------
        # auto/existing share the per-user stable dir; mode "new" gets a
        # per-process private dir so it cannot fight the shared daemon
        # over pid.lock / request files.
        if session_mode == "new":
            drop = _user_daemon_dir().parent / f"cubit-session-{os.getpid()}"
        else:
            drop = _user_daemon_dir()
        drop.mkdir(parents=True, exist_ok=True)
        outbox = drop / "out"
        outbox.mkdir(exist_ok=True)
        ready_path = drop / "ready"
        pid_file = drop / "pid.lock"

        # --- Attach to existing daemon if alive -------------------------
        existing = (None if session_mode == "new"
                    else _try_attach_existing_daemon(drop))
        if session_mode == "existing" and existing is None:
            raise CubitSessionError(
                f"RADIA_CUBIT_SESSION_MODE=existing but no live shared "
                f"Cubit daemon was found at {drop}. Start one (any "
                "cubit_show/cubit_exec in auto mode), or switch the mode.")
        if existing is not None:
            self._drop_dir = drop
            self._outbox = outbox
            self._proc = None  # not OUR child; we're just a client
            self._owned = False
            # Responsiveness probe (MathWorks ping-before-handing-out-
            # the-client): PID-alive passes for a HUNG or busy Cubit
            # whose Qt loop is not polling the drop dir.  Fail loud NOW
            # with a clear message instead of letting the first real
            # call sit in a 60 s timeout.  Non-destructive: the daemon
            # is left running (it may be mid-operation in another
            # window); _ready_info stays unset so the next call re-pings.
            ping_req = {"id": self._next_id, "op": "ping", "args": [],
                        "protocol_version": PROTOCOL_VERSION}
            self._next_id += 1
            try:
                self._call_via_filedrop(ping_req,
                                        timeout_s=ATTACH_PING_TIMEOUT_S)
            except CubitSessionError:
                pid = existing.get("pid")
                raise CubitSessionError(
                    f"Attached Cubit daemon (pid={pid}) is alive but did "
                    f"not answer a ping within {ATTACH_PING_TIMEOUT_S:.0f} s "
                    "-- it is busy with a long operation or hung. Wait and "
                    "retry, or force-stop it with cubit_session_shutdown "
                    "if it is truly stuck.")
            self._ready_info = existing
            self._last_license_warmup = {"status": "skipped",
                                          "reason": "attached to existing daemon"}
            return existing

        # --- No live daemon: new spawn path ----------------------------
        # Resolve the launcher only on the spawn path -- attaching to an
        # already-running daemon must not require (or fail on) the exe.
        gui_exe = _cubit_gui_exe(self._bin_dir)
        bootstrap_path = Path(__file__).with_name("bootstrap.py")
        if not bootstrap_path.exists():
            raise CubitSessionError(
                f"Bootstrap script missing: {bootstrap_path}. Expected "
                "sibling of session.py.")

        # Step 1: license pre-warm
        try:
            from .license_warmup import warmup_license
            warmup = warmup_license(self._bin_dir,
                                     timeout_s=LICENSE_WARMUP_TIMEOUT_S)
            self._last_license_warmup = warmup
        except Exception as exc:
            self._last_license_warmup = {
                "status": "error", "reason": f"warmup module: {exc}"}

        # Clear stale artifacts from any previous spawn in this dir
        startup_error_path = drop / "startup_error.txt"
        try:
            if ready_path.exists():
                ready_path.unlink()
            if pid_file.exists():
                pid_file.unlink()
            if startup_error_path.exists():
                startup_error_path.unlink()
            for f in outbox.iterdir():
                try: f.unlink()
                except OSError: pass
            for f in drop.glob("*.json"):
                try: f.unlink()
                except OSError: pass
        except OSError:
            pass
        self._drop_dir = drop
        self._outbox = outbox

        env = os.environ.copy()
        env["CUBIT_BIN_DIR"] = str(self._bin_dir)
        env["CUBIT_DROP_DIR"] = str(drop)
        # Cubit execs bootstrap.py as a STRING (frames show <string>),
        # so __file__-based sibling imports (probe_ops) resolve against
        # Cubit's CWD and fail -- found by the 2026-08-05 GUI E2E. The
        # package dir is passed explicitly instead.
        env["RADIA_MCP_CUBIT_PKG_DIR"] = str(Path(__file__).parent)

        # Detach the Cubit subprocess so it outlives THIS MCP server
        # process.  DETACHED_PROCESS + CREATE_NEW_PROCESS_GROUP mean: no
        # orphan child on VSCode exit, Cubit keeps running, and the next
        # MCP server process finds it via pid.lock + attaches.
        # EXCEPTION -- session_mode "new": the private daemon is useless
        # without this client, so it joins a kill-on-close Job Object
        # instead (found by the same GUI E2E: a client that died before
        # shutdown left license-holding orphans behind).
        creationflags = 0
        if sys.platform == "win32":
            creationflags = (
                getattr(subprocess, "DETACHED_PROCESS", 0) |
                getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0))
        # Cubit's own console output goes to log FILES in the drop dir
        # (MathWorks matlab_stdout.log pattern), NOT to pipes (pipes
        # block/break when the parent exits -- the detached child must
        # outlive us) and NOT to DEVNULL (which threw away exactly the
        # output needed when bootstrap dies before writing `ready`).
        stdout_log = drop / "cubit_stdout.log"
        stderr_log = drop / "cubit_stderr.log"
        with open(stdout_log, "wb") as out_f, open(stderr_log, "wb") as err_f:
            self._proc = subprocess.Popen(
                [str(gui_exe), "-nojournal", str(bootstrap_path)],
                # stdin=DEVNULL: prevent Cubit from inheriting the MCP
                # server's stdio pipe (Claude Code JSON-RPC). Without it
                # Cubit hangs at startup (resp=False, declining threads).
                # Confirmed 2026-04-25.
                stdin=subprocess.DEVNULL,
                stdout=out_f,
                stderr=err_f,
                env=env,
                bufsize=0,
                creationflags=creationflags,
                close_fds=True,
            )
        self._owned = True
        if session_mode == "new" and sys.platform == "win32":
            try:
                self._job_handle = _assign_kill_on_close_job(self._proc.pid)
            except OSError as exc:
                try:
                    self._proc.kill()
                    self._proc.wait(timeout=3.0)
                except Exception:
                    pass
                self._proc = None
                self._owned = False
                _remove_tree_with_retry(drop)
                self._drop_dir = None
                self._outbox = None
                raise CubitSessionError(
                    "Could not place the private Cubit process in a "
                    "kill-on-close Windows Job Object; refusing to leave "
                    "a mode-new session without orphan protection: "
                    f"{exc}") from exc
        # Write PID lock immediately so a sibling VSCode racing to spawn
        # can also attach (once ready marker appears).
        try:
            pid_file.write_text(str(self._proc.pid), encoding="utf-8")
        except OSError:
            pass
        # NOTE: no stderr drain needed — DEVNULL means no pipe to drain.

        # --- Step 2: ready wait with shorter timeout after warmup ---
        # After license pre-warm, Cubit should be ready in ~3 s.  If
        # the warmup was skipped (cache fresh, login not needed), the
        # budget still includes the full license checkout.  So the
        # actual deadline depends on whether we ran the warmup.
        warmup_ok = (self._last_license_warmup.get("status") in ("ok", "skipped")
                     and self._last_license_warmup.get("action")
                         != "no_rlm_activate")
        budget = CUBIT_READY_TIMEOUT_S if warmup_ok else 60.0
        deadline = time.time() + budget

        def _startup_diag() -> str:
            """Real in-process failure evidence, best first (MathWorks
            mcp_startup_error.txt pattern + log-artifact pointers)."""
            parts = []
            if startup_error_path.exists():
                try:
                    parts.append("bootstrap error:\n"
                                 + startup_error_path.read_text(
                                       encoding="utf-8", errors="replace"))
                except OSError:
                    pass
            try:
                tail = stderr_log.read_bytes()[-2000:]
                if tail.strip():
                    parts.append("cubit stderr tail:\n"
                                 + tail.decode("utf-8", errors="replace"))
            except OSError:
                pass
            parts.append(f"logs: {stdout_log} / {stderr_log}")
            return "\n".join(parts)

        while time.time() < deadline:
            if startup_error_path.exists():
                raise CubitSessionError(
                    "Cubit bootstrap failed during startup.\n"
                    + _startup_diag())
            if ready_path.exists():
                try:
                    info = json.loads(ready_path.read_text(encoding="utf-8"))
                except json.JSONDecodeError:
                    time.sleep(0.05)
                    continue
                self._ready_info = info
                return info
            if self._proc.poll() is not None:
                raise CubitSessionError(
                    f"Cubit exited during bootstrap "
                    f"(rc={self._proc.returncode}).\n" + _startup_diag())
            time.sleep(0.15)
        raise CubitSessionError(
            f"Cubit did not signal ready within {budget:.0f} s. "
            f"drop={drop}\n" + _startup_diag())

    def _request_stem(self, req_id: int) -> str:
        """Return a request filename unique across attached MCP clients."""
        client_id = getattr(self, "_client_id", None)
        if client_id is None:
            client_id = f"{os.getpid():08x}-{uuid.uuid4().hex[:12]}"
            self._client_id = client_id
        return f"{client_id}-{int(req_id):08d}"

    def _call_via_filedrop(self, req: dict, timeout_s: float) -> dict:
        assert self._drop_dir is not None and self._outbox is not None
        req_id = req["id"]
        stem = self._request_stem(req_id)
        req_path = self._drop_dir / f"{stem}.req"
        tmp_path = self._drop_dir / f"{stem}.req.tmp"
        resp_path = self._outbox / f"{stem}.resp"

        tmp_path.write_text(json.dumps(req), encoding="utf-8")
        os.replace(tmp_path, req_path)

        deadline = time.time() + timeout_s
        while time.time() < deadline:
            if resp_path.exists():
                try:
                    text = resp_path.read_text(encoding="utf-8")
                    resp = json.loads(text)
                except (OSError, json.JSONDecodeError):
                    # Partial write: back off and retry.
                    time.sleep(0.03)
                    continue
                resp_path.unlink(missing_ok=True)
                if resp.get("id") != req_id:
                    raise CubitSessionError(
                        f"Response id mismatch for op={req.get('op')!r}: "
                        f"expected {req_id}, got {resp.get('id')!r}")
                return resp
            if self._proc is not None and self._proc.poll() is not None:
                raise CubitSessionError(
                    f"Cubit exited during call (rc={self._proc.returncode}).")
            time.sleep(0.05)
        raise CubitSessionError(
            f"Response timeout after {timeout_s}s for op={req.get('op')!r}")

    # ---- mode: batch (stdio JSON-RPC) ----

    def _start_stdio_daemon(self) -> dict:
        python_exe = _cubit_python_exe(self._bin_dir)
        daemon_path = Path(__file__).with_name("daemon.py")
        if not daemon_path.exists():
            raise CubitSessionError(
                f"Daemon script missing: {daemon_path}. Expected sibling "
                "of session.py in the mcp-server-cubit package.")

        env = os.environ.copy()
        env["CUBIT_DAEMON_MODE"] = self._mode
        env["CUBIT_BIN_DIR"] = str(self._bin_dir)

        self._proc = subprocess.Popen(
            [str(python_exe), str(daemon_path)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            bufsize=0,
        )
        self._owned = True
        self._start_stderr_drain()

        line = self._proc.stdout.readline()
        if not line:
            stderr = self._proc.stderr.read().decode("utf-8", errors="replace")
            raise CubitSessionError(
                f"Daemon died during startup. stderr:\n{stderr[-2000:]}")
        try:
            ready = json.loads(line.decode("utf-8"))
        except json.JSONDecodeError as e:
            raise CubitSessionError(
                f"Daemon emitted non-JSON ready line: {line!r} ({e})")
        if not ready.get("ready"):
            raise CubitSessionError(
                f"Daemon reported startup failure: {ready}")
        self._ready_info = ready
        return ready

    def _call_via_stdio(self, req: dict, timeout_s: float) -> dict:
        self._send_stdio(req)
        line = self._read_stdio_line(timeout_s=timeout_s)
        resp = json.loads(line.decode("utf-8"))
        if resp.get("id") is not None and resp["id"] != req["id"]:
            raise CubitSessionError(
                f"Response id {resp.get('id')} != request id {req['id']}")
        return resp

    def _send_stdio(self, obj: dict) -> None:
        assert self._proc is not None and self._proc.stdin is not None
        line = (json.dumps(obj) + "\n").encode("utf-8")
        try:
            self._proc.stdin.write(line)
            self._proc.stdin.flush()
        except (BrokenPipeError, OSError) as e:
            raise CubitSessionError(f"Daemon pipe broken: {e}")

    def _read_stdio_line(self, timeout_s: float) -> bytes:
        assert self._proc is not None and self._proc.stdout is not None
        result: list[bytes | Exception] = []

        def _reader():
            try:
                line = self._proc.stdout.readline()
                result.append(line)
            except Exception as e:
                result.append(e)

        t = threading.Thread(target=_reader, daemon=True)
        t.start()
        t.join(timeout=timeout_s)
        if t.is_alive():
            raise CubitSessionError(
                f"Daemon response timed out after {timeout_s}s")
        if not result:
            raise CubitSessionError("Daemon returned no data (pipe closed?)")
        out = result[0]
        if isinstance(out, Exception):
            raise CubitSessionError(f"Read error: {out}")
        if not out:
            stderr_tail = b"".join(self._stderr_tail)
            raise CubitSessionError(
                f"Daemon exited unexpectedly (exit={self._proc.poll()}). "
                f"stderr tail:\n"
                f"{stderr_tail.decode('utf-8', errors='replace')[-1500:]}")
        return out

    # ---- common: stderr drain, shutdown op ----

    def _start_stderr_drain(self) -> None:
        """Drain the subprocess' stderr into a bounded ring buffer.

        If nobody reads stderr, Cubit's heavy progress output fills the
        pipe buffer (~64 KB) and subsequent writes block, freezing
        the whole session. Retained for post-mortem only.
        """
        assert self._proc is not None and self._proc.stderr is not None
        self._stderr_tail = []

        def _drain(pipe, buf, buf_max):
            try:
                for line in iter(pipe.readline, b""):
                    buf.append(line)
                    if len(buf) > buf_max:
                        del buf[0]
            except Exception:
                pass

        t = threading.Thread(
            target=_drain,
            args=(self._proc.stderr, self._stderr_tail, self._stderr_tail_max),
            daemon=True,
        )
        t.start()

    def _lock_free_shutdown_op(self, timeout_s: float) -> None:
        """Fire-and-forget shutdown op — tolerate any error."""
        try:
            req_id = self._next_id
            self._next_id += 1
            req = {"id": req_id, "op": "shutdown", "args": []}
            if self._mode == "gui":
                try:
                    self._call_via_filedrop(req, timeout_s=timeout_s)
                except CubitSessionError:
                    pass
            else:
                self._send_stdio(req)
        except Exception:
            pass

    # ---- singleton access ----

    @classmethod
    def get(cls, mode: str = "gui") -> "CubitSession":
        """Return the process-wide singleton, creating on first call."""
        global _SINGLETON
        with _SESSION_LOCK:
            if _SINGLETON is None:
                _SINGLETON = cls(mode=mode)
            return _SINGLETON

    @classmethod
    def reset(cls) -> dict:
        """Shutdown and drop the singleton. Next `get()` will relaunch.

        Returns the shutdown report ({"stopped": ..., "pid": ...}); an
        empty dict when no singleton existed."""
        global _SINGLETON
        report: dict = {}
        with _SESSION_LOCK:
            if _SINGLETON is not None:
                report = _SINGLETON.shutdown()
                _SINGLETON = None
        return report
