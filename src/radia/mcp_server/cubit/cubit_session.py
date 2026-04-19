"""
cubit_session.py — Python-3.12-side client for the Cubit Viewer daemon.

This module is version-agnostic: stdlib only, works on any Python
>= 3.10.  The only version-sensitive piece is the daemon subprocess
itself (Python 3.10 because `_cubit3.pyd` has a Python-ABI dependency),
and its path is located at runtime via `find_cubit_install()`, not
hardcoded.

Architecture (OCP-inspired, transposed to Cubit):

    mcp-server-cubit  (this process, Python 3.12)
        │  JSON-RPC over stdin/stdout
        ▼
    cubit_daemon.py  (subprocess, Cubit-Python 3.10)
        │  cubit.cmd(...)
        ▼
    Cubit GUI  (persistent FLTK window)

Unlike OCP (stateless per-call connect/send/close), our transport is
a persistent stdio pipe because launching Cubit per call would be
prohibitively slow (~5 s cold vs <1 ms hot). The singleton session
holds the subprocess handle and pipes for the lifetime of the
mcp-server process.

Protocol version: 1. Stored in `PROTOCOL_VERSION` — increment if the
JSON-RPC schema changes in a breaking way; the daemon echoes the
version in its ready message for mutual compatibility check.
"""

from __future__ import annotations

import glob
import json
import os
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any


PROTOCOL_VERSION = 1

_SESSION_LOCK = threading.Lock()
_SINGLETON: "CubitSession | None" = None

# Module-level overrides (OCP-inspired: cf. set_port / get_port).
# Use these instead of mutating env vars at runtime.
_OVERRIDE_BIN_DIR: Path | None = None


def set_cubit_bin_dir(path: str | Path | None) -> None:
    """Override the Cubit `bin/` directory lookup for subsequent sessions.

    Mirrors OCP's `set_port(...)` accessor pattern. Calling with None
    restores auto-discovery. Must be called before `CubitSession.get()`
    to take effect on the singleton.
    """
    global _OVERRIDE_BIN_DIR
    _OVERRIDE_BIN_DIR = Path(path) if path is not None else None


def get_cubit_bin_dir() -> Path | None:
    """Return the currently-resolved Cubit `bin/` directory, or None."""
    if _OVERRIDE_BIN_DIR is not None:
        return _OVERRIDE_BIN_DIR
    return find_cubit_install()


def _session_info_path() -> Path:
    """Location for the per-user daemon diagnostic file (OCP-inspired).

    Mirrors OCP's `set_connection_file()` pattern. Used only for
    debugging / external tooling that wants to know if a daemon is
    running and under which pid.
    """
    return Path.home() / ".cubit_viewer" / "session.json"


# ---------------------------------------------------------------------------
# Cubit install auto-discovery (no hardcoded version in a path)
# ---------------------------------------------------------------------------

def find_cubit_install(explicit: str | None = None) -> Path | None:
    """Locate a Coreform Cubit install's `bin/` directory.

    Resolution order (first hit wins):
        1. `explicit` argument (if provided)
        2. `CUBIT_BIN_DIR` environment variable
        3. Windows registry `HKLM\\SOFTWARE\\Coreform\\*\\InstallDir`
        4. Glob `C:/Program Files/Coreform Cubit */bin/`
        5. Glob `/opt/Coreform-Cubit-*/bin/` (Linux, future-proof)
        6. Return None (caller reports error)

    Returns the absolute path to the `bin/` directory (containing
    `coreform_cubit.exe`, `cubit.py`, and `_cubit3.pyd`), or None.
    """
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit))
    env_bin = os.environ.get("CUBIT_BIN_DIR")
    if env_bin:
        candidates.append(Path(env_bin))

    # Windows registry lookup (best-effort)
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

    # Filesystem globs
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
    """Locate Cubit's bundled Python interpreter under `bin_dir`.

    Standard layout on Windows (Cubit 2025.x): `<bin>/python3/python.exe`.
    Linux variants may differ; we also check `<bin>/python3/bin/python3`.
    """
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


# ---------------------------------------------------------------------------
# Session (singleton per mcp-server process)
# ---------------------------------------------------------------------------

class CubitSessionError(Exception):
    """Raised when the Cubit daemon cannot be launched or has died."""


class CubitSession:
    """Manages a long-lived Cubit daemon subprocess via JSON-RPC.

    Thread-safe via an internal lock around each call (Cubit itself is
    single-threaded; serializing at the client side avoids reentrancy).

    Usage from mcp tools:
        session = CubitSession.get()
        r = session.call("cmd", ["create brick x 10", "mesh volume 1"])
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
        self._lock = threading.Lock()
        self._ready_info: dict | None = None

    # ---- lifecycle ----

    def ensure_started(self) -> dict:
        """Spawn the daemon if not already running. Returns ready info."""
        if self._proc is not None and self._proc.poll() is None:
            return self._ready_info or {}

        python_exe = _cubit_python_exe(self._bin_dir)
        daemon_path = Path(__file__).with_name("cubit_daemon.py")
        if not daemon_path.exists():
            raise CubitSessionError(
                f"Daemon script missing: {daemon_path}. Expected sibling "
                "of cubit_session.py in the mcp-server-cubit package.")

        env = os.environ.copy()
        env["CUBIT_DAEMON_MODE"] = self._mode
        # Help the daemon find its own cubit binding even if CUBIT_BIN_DIR
        # is not preset in the user's shell.
        env["CUBIT_BIN_DIR"] = str(self._bin_dir)

        # Use binary mode for stdin/stdout so we can control newline
        # handling precisely across platforms.
        self._proc = subprocess.Popen(
            [str(python_exe), str(daemon_path)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            bufsize=0,
        )

        # Read daemon's "ready" line.
        line = self._proc.stdout.readline()
        if not line:
            stderr = self._proc.stderr.read().decode("utf-8", errors="replace")
            raise CubitSessionError(
                f"Daemon died during startup. stderr:\n{stderr[-2000:]}"
            )
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

    def shutdown(self, timeout_s: float = 3.0) -> None:
        """Politely ask the daemon to exit, then force-kill after timeout."""
        if self._proc is None:
            return
        try:
            self._send({"id": self._next_id, "op": "shutdown", "args": []})
            self._proc.wait(timeout=timeout_s)
        except Exception:
            pass
        if self._proc.poll() is None:
            self._proc.kill()
        self._proc = None
        self._ready_info = None

    def is_alive(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    def ping(self, timeout_s: float = 2.0) -> bool:
        """Heartbeat probe — returns True if daemon responds."""
        try:
            r = self.call("ping", timeout_s=timeout_s)
            return bool(r.get("ok") and r.get("result") == "pong")
        except Exception:
            return False

    # ---- RPC ----

    def call(self, op: str, args: list | None = None,
             timeout_s: float = 60.0) -> dict:
        """Send one JSON-RPC request, return the response dict.

        Response shape:
            {"id": int, "ok": bool, "result": ..., "error": str?}

        On daemon death, raises CubitSessionError. Caller may retry
        after calling `shutdown()` to reset state (rare).
        """
        with self._lock:
            self.ensure_started()
            req_id = self._next_id
            self._next_id += 1
            req = {
                "id": req_id,
                "op": op,
                "args": args or [],
                "protocol_version": PROTOCOL_VERSION,
            }
            self._send(req)
            line = self._read_line(timeout_s=timeout_s)
            resp = json.loads(line.decode("utf-8"))
            # Basic sanity: id should match
            if resp.get("id") is not None and resp["id"] != req_id:
                raise CubitSessionError(
                    f"Response id {resp.get('id')} != request id {req_id}"
                )
            return resp

    # ---- internal ----

    def _send(self, obj: dict) -> None:
        assert self._proc is not None
        line = (json.dumps(obj) + "\n").encode("utf-8")
        try:
            self._proc.stdin.write(line)
            self._proc.stdin.flush()
        except (BrokenPipeError, OSError) as e:
            raise CubitSessionError(f"Daemon pipe broken: {e}")

    def _read_line(self, timeout_s: float) -> bytes:
        """Read one newline-terminated frame from the daemon's stdout.

        Uses a background thread to implement a timeout since
        `subprocess.Popen.stdout.readline()` doesn't honor timeouts.
        """
        assert self._proc is not None
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
            # Daemon hung; caller should probably shutdown and restart
            raise CubitSessionError(
                f"Daemon response timed out after {timeout_s}s")
        if not result:
            raise CubitSessionError("Daemon returned no data (pipe closed?)")
        out = result[0]
        if isinstance(out, Exception):
            raise CubitSessionError(f"Read error: {out}")
        if not out:
            # EOF: daemon exited
            stderr_tail = b""
            try:
                stderr_tail = self._proc.stderr.read() or b""
            except Exception:
                pass
            raise CubitSessionError(
                f"Daemon exited unexpectedly (exit={self._proc.poll()}). "
                f"stderr tail:\n"
                f"{stderr_tail.decode('utf-8', errors='replace')[-1500:]}")
        return out

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
    def reset(cls) -> None:
        """Shutdown and drop the singleton. Next `get()` will relaunch."""
        global _SINGLETON
        with _SESSION_LOCK:
            if _SINGLETON is not None:
                _SINGLETON.shutdown()
                _SINGLETON = None
