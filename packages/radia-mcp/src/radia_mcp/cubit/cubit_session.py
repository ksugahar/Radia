"""
cubit_session.py — Python-3.12-side client for the Cubit Viewer daemon.

This module is version-agnostic: stdlib only, works on any Python
>= 3.10.

Two transport modes share the same public API (call/ping/shutdown):

    mode="gui"    (DEFAULT, Plan A — 2026-04-19 file-drop bootstrap)
        ─ Launches `coreform_cubit.exe -nojournal cubit_bootstrap.py`
        ─ Bootstrap installs a QTimer inside Cubit's Qt event loop
        ─ Client ↔ Bootstrap communicate via atomically-renamed JSON
          files in a per-session drop directory (no sockets)
        ─ Cubit GUI window is live & user-interactive throughout

    mode="batch"  (legacy, CI/scripting)
        ─ Launches Cubit's bundled Python 3.10 with cubit_daemon.py
        ─ Client ↔ Daemon communicate via line-delimited JSON-RPC on
          stdin/stdout of the subprocess
        ─ No GUI, no user interaction; pure headless

Both modes preserve a persistent Cubit session across many client calls
(~2 s cold start amortized) — launching per-call would be prohibitively
slow. A process-wide singleton (`CubitSession.get()`) holds the handle.

Protocol version:
    v1 — stdio JSON-RPC (cubit_daemon.py)
    v2 — file drop JSON-RPC (cubit_bootstrap.py, this module default)

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
from pathlib import Path
from typing import Any


PROTOCOL_VERSION = 2  # default (gui mode)

_SESSION_LOCK = threading.Lock()
_SINGLETON: "CubitSession | None" = None

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
        self._lock = threading.Lock()
        self._ready_info: dict | None = None

        # gui-mode (file-drop) state
        self._drop_dir: Path | None = None
        self._outbox: Path | None = None

        # batch-mode (stdio) stderr retention
        self._stderr_tail: list[bytes] = []
        self._stderr_tail_max = 200

    # ---- lifecycle ----

    def ensure_started(self) -> dict:
        """Spawn the Cubit session if not already running. Returns ready info."""
        if self._proc is not None and self._proc.poll() is None:
            return self._ready_info or {}
        if self._mode == "gui":
            return self._start_gui_bootstrap()
        return self._start_stdio_daemon()

    def shutdown(self, timeout_s: float = 3.0) -> None:
        """Politely ask Cubit to exit, then force-kill after timeout."""
        if self._proc is None:
            return
        try:
            # Best-effort: issue shutdown op via the live transport.
            self._lock_free_shutdown_op(timeout_s=timeout_s)
            self._proc.wait(timeout=timeout_s)
        except Exception:
            pass
        if self._proc.poll() is None:
            self._proc.kill()
        # Clean up gui drop dir
        if self._mode == "gui" and self._drop_dir is not None:
            try:
                shutil.rmtree(self._drop_dir, ignore_errors=True)
            except Exception:
                pass
            self._drop_dir = None
            self._outbox = None
        self._proc = None
        self._ready_info = None

    def is_alive(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

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
                    return self._call_via_filedrop(req, timeout_s=timeout_s)
                return self._call_via_stdio(req, timeout_s=timeout_s)
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

    def _force_reset(self) -> None:
        """Tear down the current session without waiting for graceful exit.

        Used by the recovery path in `call()`. Does NOT re-enter the lock
        (caller is holding it).
        """
        proc = self._proc
        if proc is not None and proc.poll() is None:
            try:
                proc.kill()
            except Exception:
                pass
        if self._mode == "gui" and self._drop_dir is not None:
            try:
                shutil.rmtree(self._drop_dir, ignore_errors=True)
            except Exception:
                pass
            self._drop_dir = None
            self._outbox = None
        self._proc = None
        self._ready_info = None

    # ---- mode: gui (file drop) ----

    def _start_gui_bootstrap(self) -> dict:
        """Launch coreform_cubit.exe + cubit_bootstrap.py, wait for ready."""
        gui_exe = _cubit_gui_exe(self._bin_dir)
        bootstrap_path = Path(__file__).with_name("cubit_bootstrap.py")
        if not bootstrap_path.exists():
            raise CubitSessionError(
                f"Bootstrap script missing: {bootstrap_path}. Expected "
                "sibling of cubit_session.py.")

        drop = Path(tempfile.mkdtemp(prefix="cubit_drop_"))
        (drop / "out").mkdir(exist_ok=True)
        self._drop_dir = drop
        self._outbox = drop / "out"

        env = os.environ.copy()
        env["CUBIT_BIN_DIR"] = str(self._bin_dir)
        env["CUBIT_DROP_DIR"] = str(drop)

        # stdout/stderr go to pipes we drain, same as batch mode — Cubit
        # GUI still writes banner / warnings to these before the Qt
        # widgets take over.
        self._proc = subprocess.Popen(
            [str(gui_exe), "-nojournal", str(bootstrap_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            bufsize=0,
        )
        self._start_stderr_drain()

        # Wait for ready marker. Cold start ~2s on this lab, but allow
        # generous slack for first-license-checkout.
        ready_path = drop / "ready"
        deadline = time.time() + 60.0
        while time.time() < deadline:
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
                    f"(rc={self._proc.returncode}).")
            time.sleep(0.15)
        raise CubitSessionError(
            f"Cubit did not signal ready within 60 s. drop={drop}")

    def _call_via_filedrop(self, req: dict, timeout_s: float) -> dict:
        assert self._drop_dir is not None and self._outbox is not None
        req_id = req["id"]
        stem = f"{req_id:08d}"
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
        daemon_path = Path(__file__).with_name("cubit_daemon.py")
        if not daemon_path.exists():
            raise CubitSessionError(
                f"Daemon script missing: {daemon_path}. Expected sibling "
                "of cubit_session.py in the mcp-server-cubit package.")

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
    def reset(cls) -> None:
        """Shutdown and drop the singleton. Next `get()` will relaunch."""
        global _SINGLETON
        with _SESSION_LOCK:
            if _SINGLETON is not None:
                _SINGLETON.shutdown()
                _SINGLETON = None
