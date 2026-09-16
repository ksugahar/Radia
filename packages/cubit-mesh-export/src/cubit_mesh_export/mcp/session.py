"""Process-owned headless Cubit session using daemon.py stdio protocol v1.

Calls reuse one child within this MCP process. GUI attachment and persistent
file-drop sessions are not supported; human GUI work exchanges saved artifacts.
"""

from __future__ import annotations

import glob
import json
import math
import os
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from pathlib import Path
from typing import Any

PROTOCOL_VERSION = 1


def _cubit_temp_root() -> Path:
    """Honor the caller's scratch root on every supported platform."""
    configured = os.environ.get("CUBIT_MCP_TEMP")
    if configured is not None:
        return Path(configured)
    return Path("C:/temp" if sys.platform == "win32" else tempfile.gettempdir())

# License checkout may take tens of seconds on a cold start.
CUBIT_READY_TIMEOUT_S = 90.0

_SESSION_LOCK = threading.Lock()
_SINGLETON: "CubitSession | None" = None


def _close_proc_streams(proc) -> None:
    """Close a reaped child's pipe handles (batch mode uses three PIPEs;
    leaving them to GC raises unclosed-FileIO ResourceWarnings, which
    pytest promotes to errors)."""
    for stream in (proc.stdin, proc.stdout, proc.stderr):
        try:
            if stream is not None:
                stream.close()
        except Exception:
            pass


def _assign_kill_on_close_job(pid: int):
    """Put ``pid`` into a Windows Job Object with KILL_ON_JOB_CLOSE.

    The process-owned daemon must die with its client. The returned job
    HANDLE must stay referenced for the client's lifetime; when this
    process exits (normally or not), the OS closes the handle and kills
    the job's processes.

    Returns the handle. Failure raises ``OSError``: the client must
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
    env_install = os.environ.get("CUBIT_INSTALL_DIR")
    if env_install:
        candidates.append(Path(env_install) / "bin")

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


def run_headless_journal(
    commands: list[str],
    *,
    timeout_s: float = 300.0,
    working_directory: str | Path | None = None,
    command_plugin_directory: str | Path | None = None,
) -> dict:
    """Run plugin-aware Cubit commands in a disposable headless process.

    The bundled-Python daemon is sufficient for native ``cubit.cmd`` calls,
    but it does not load Cubit's C++ command plugins.  Commands such as
    ``export netgen`` therefore need the real console launcher.  On Windows
    ``coreform_cubit.com`` is deliberately preferred over the GUI-stub
    ``.exe`` so callers can wait for completion and capture diagnostics.

    ``command_plugin_directory`` selects a rebuilt plugin with Cubit's official
    ``-commandplugindir`` switch.  This is the safe plugin-under-test route when
    another user has the deployed plugin open: it does not replace or unload
    the system binary, and it suppresses the user's Cubit init file.
    """
    if not commands:
        raise ValueError("commands must not be empty")
    if any(not isinstance(line, str) or not line.strip() or "\n" in line
           or "\r" in line for line in commands):
        raise ValueError("commands must be nonempty single-line strings")
    try:
        timeout_s = float(timeout_s)
    except (TypeError, ValueError) as exc:
        raise ValueError("timeout_s must be a finite positive number") from exc
    if not math.isfinite(timeout_s) or timeout_s <= 0.0:
        raise ValueError("timeout_s must be a finite positive number")

    bin_dir = get_cubit_bin_dir()
    if bin_dir is None:
        return {
            "status": "error", "stage": "start", "kind": "environment",
            "gui_started": False,
            "error": "Could not locate Coreform Cubit install",
        }
    console = bin_dir / "coreform_cubit.com"
    if not console.exists():
        return {
            "status": "error", "stage": "start", "kind": "environment",
            "gui_started": False,
            "error": (
                f"Headless Cubit console not found: {console}. "
                "Refusing to fall back to the GUI launcher."
            ),
        }

    temp_root = _cubit_temp_root()
    temp_root.mkdir(parents=True, exist_ok=True)
    cwd = Path(working_directory) if working_directory else temp_root
    cwd.mkdir(parents=True, exist_ok=True)
    plugin_dir = None
    if command_plugin_directory is not None:
        plugin_dir = Path(command_plugin_directory)
        if not plugin_dir.is_dir():
            return {
                "status": "error", "stage": "preflight", "kind": "input",
                "gui_started": False,
                "error": f"command plugin directory not found: {plugin_dir}",
            }

    with tempfile.TemporaryDirectory(
            prefix="cubit_mesh_export_headless_", dir=str(temp_root)) as scratch:
        driver = Path(scratch) / "driver.jou"
        driver.write_text("\n".join([*commands, "exit 0", ""]),
                          encoding="utf-8")
        argv = [str(console), "-nographics", "-batch", "-nojournal"]
        if plugin_dir is not None:
            argv.extend(["-noinitfile", "-commandplugindir", str(plugin_dir)])
        argv.append(str(driver))
        try:
            proc = subprocess.run(
                argv, cwd=str(cwd), capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=timeout_s,
                # Cubit reads from an inherited console/stdin during startup.
                # Under an MCP stdio server that handle is the JSON-RPC pipe,
                # which leaves an otherwise headless batch launch waiting
                # indefinitely.  Give the disposable process no input source.
                stdin=subprocess.DEVNULL,
            )
        except subprocess.TimeoutExpired as exc:
            return {
                "status": "error", "stage": "timeout", "kind": "timeout",
                "gui_started": False,
                "error": f"Cubit headless journal exceeded {timeout_s}s",
                "timeout_s": timeout_s,
                "stdout_tail": (exc.stdout or "")[-4000:],
                "stderr_tail": (exc.stderr or "")[-4000:],
            }
        except OSError as exc:
            return {
                "status": "error", "stage": "start", "kind": "environment",
                "gui_started": False,
                "error": f"{type(exc).__name__}: {exc}",
            }

    return {
        # Cubit 2025.12 can return nonzero after successful exports because
        # of classified startup/teardown diagnostics.  Artifact freshness is
        # the caller's completion gate; retain the exit code as evidence.
        "status": "completed",
        "exit_code": int(proc.returncode),
        "console": str(console),
        "headless_flags": [
            "-nographics", "-batch", "-nojournal",
            *(["-noinitfile", "-commandplugindir"] if plugin_dir else []),
        ],
        "command_plugin_directory": str(plugin_dir) if plugin_dir else None,
        "user_init_loaded": plugin_dir is None,
        "gui_started": False,
        "command_count": len(commands),
        "stdout_tail": (proc.stdout or "")[-8000:],
        "stderr_tail": (proc.stderr or "")[-8000:],
    }


# ---------------------------------------------------------------------------
# Session (singleton per mcp-server process)
# ---------------------------------------------------------------------------

class CubitSessionError(Exception):
    """Raised when the Cubit daemon cannot be launched or has died."""


class CubitSession:
    """Manages one process-owned headless Cubit session.

    Thread-safe via an internal lock around each call.

    Usage:
        session = CubitSession.get()   # mode="batch" by default
        r = session.call("cmd", ["create brick x 10"])
        if r["ok"]:
            ...
    """

    def __init__(self,
                 cubit_bin_dir: Path | None = None,
                 mode: str = "batch"):
        if mode != "batch":
            raise ValueError(f"Only headless mode='batch' is supported, got {mode!r}")
        self._bin_dir = cubit_bin_dir or get_cubit_bin_dir()
        if self._bin_dir is None:
            raise CubitSessionError(
                "Could not locate Coreform Cubit install. Set CUBIT_BIN_DIR "
                "or install to a standard location."
            )
        self._mode = mode
        self._proc: subprocess.Popen | None = None
        self._next_id = 1
        self._client_id = f"{os.getpid():08x}-{uuid.uuid4().hex[:12]}"
        self._lock = threading.RLock()
        self._ready_info: dict | None = None

        self._owned = False
        self._job_handle = None

        # Per-MCP-process command history for `cubit_session_journal`.
        # A live session can be captured explicitly as a portable recipe;
        # export itself never requires or creates that journal. Bounded;
        # records {ts, line, ok} for every op=="cmd" line sent.
        self._command_history: list[dict] = []
        self._command_history_max = 20000

        # Cubit's own command record is the provenance source of truth.  The
        # in-memory response history above is retained only for diagnostics
        # (not for reconstructing a journal).  A new file is started after
        # each daemon generation so recovery never splices unlike sessions.
        self._native_journal_path: Path | None = None
        self._native_journal_paths: list[Path] = []
        self._native_journal_error: str | None = None

        # batch-mode (stdio) stderr retention
        self._stderr_tail: list[bytes] = []
        self._stderr_tail_max = 200

    def _close_private_job(self) -> None:
        """Release the Job Object that contains this client's child."""
        handle, self._job_handle = getattr(self, "_job_handle", None), None
        if handle is None or sys.platform != "win32":
            return
        try:
            _close_windows_handle(handle)
        except OSError:
            pass

    def __del__(self):
        # The headless child is private to this object. Losing the object
        # must not lose the only handle that enforces kill-on-close.
        try:
            self._close_private_job()
        except Exception:
            pass

    # ---- lifecycle ----

    def ensure_started(self) -> dict:
        """Start or reuse this client's headless child; never attach by PID."""
        with self._lock:
            return self._ensure_started_locked()

    def _ensure_started_locked(self) -> dict:
        if self._proc is not None and self._proc.poll() is None:
            return self._ready_info or {}
        if self._proc is not None:
            self._force_reset()
            raise CubitSessionError(
                "Previous daemon exited; geometry state was lost. "
                "Restore a checkpoint explicitly before continuing.")
        try:
            return self._start_stdio_daemon()
        except (OSError, ValueError) as exc:
            self._force_reset()
            raise CubitSessionError(f"Daemon startup failed: {exc}") from exc
        except Exception:
            self._force_reset()
            raise

    def shutdown(self, timeout_s: float = 3.0) -> dict:
        """Stop only the child spawned by this client."""
        with self._lock:
            proc = self._proc
            report = {"stopped": "owned-child" if proc else "none",
                      "pid": proc.pid if proc else None, "owned": self._owned}
            if proc is not None:
                self._lock_free_shutdown_op(timeout_s)
                try:
                    proc.wait(timeout=timeout_s)
                except subprocess.TimeoutExpired:
                    pass
            self._force_reset()
            return report

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
             timeout_s: float = 60.0) -> dict:
        """Serialize an RPC. A failed transport is never replayed automatically.

        A command may have run before its response was lost. Resetting and
        retrying would hide both the unknown outcome and the lost geometry.
        """
        with self._lock:
            try:
                self.ensure_started()
                if op == "cmd":
                    self._start_native_journal_locked(timeout_s=timeout_s)
                req_id = self._next_id
                self._next_id += 1
                resp = self._call_via_stdio({
                    "id": req_id, "op": op, "args": args or [],
                    "protocol_version": PROTOCOL_VERSION,
                }, timeout_s=timeout_s)
                resp["execution_mode"] = "batch"
                resp["gui_started"] = False
                if op == "cmd":
                    self._record_cmd_history(resp)
                return resp
            except CubitSessionError as exc:
                self._force_reset()
                raise CubitSessionError(
                    f"{exc}. Session discarded; command outcome may be unknown. "
                    "No automatic replay. Restore a checkpoint explicitly."
                ) from exc

    def _start_native_journal_locked(self, timeout_s: float) -> None:
        """Start Cubit's native ``record \"file\"`` stream once per daemon.

        Caller holds ``self._lock``.  Cubit 2025.12 does not accept the older
        ``record journal ... overwrite`` spelling; the destination must be a
        fresh path.  Failing closed prevents an apparently reproducible AI
        session whose journal was actually reconstructed from RPC responses.
        """
        if getattr(self, "_native_journal_path", None) is not None:
            return
        paths = getattr(self, "_native_journal_paths", None)
        if paths is None:
            paths = self._native_journal_paths = []
        temp_root = _cubit_temp_root()
        journal_dir = temp_root / "cubit-mesh-export" / "journals"
        journal_dir.mkdir(parents=True, exist_ok=True)
        generation = len(paths) + 1
        path = journal_dir / (
            f"ai-{self._client_id}-generation-{generation:03d}.jou")
        if path.exists():
            raise CubitSessionError(
                f"refusing to overwrite native Cubit journal: {path}")
        command = f'record "{str(path).replace(chr(92), "/")}"'
        req_id = self._next_id
        self._next_id += 1
        request = {
            "id": req_id,
            "op": "cmd",
            "args": [command],
            "protocol_version": 1,
        }
        response = self._call_via_stdio(request, timeout_s=timeout_s)
        per_line = response.get("result") if isinstance(response, dict) else None
        ok = bool(response.get("ok")) if isinstance(response, dict) else False
        ok = ok and isinstance(per_line, list) and bool(per_line)
        ok = ok and bool(per_line[0].get("ok"))
        if not ok:
            self._native_journal_error = (
                f"Cubit rejected {command!r}: {response!r}")
            raise CubitSessionError(self._native_journal_error)
        self._native_journal_path = path
        paths.append(path)
        self._native_journal_error = None

    def native_journal_snapshot(self) -> dict:
        """Return Cubit-recorded journal generations without starting Cubit."""
        paths = list(getattr(self, "_native_journal_paths", []))
        chunks: list[str] = []
        readable: list[str] = []
        errors: list[str] = []
        for path in paths:
            try:
                chunks.append(path.read_text(encoding="utf-8-sig"))
                readable.append(str(path))
            except OSError as exc:
                errors.append(f"{path}: {exc}")
        return {
            "journal": "\n".join(chunk.rstrip("\n") for chunk in chunks)
                       + ("\n" if chunks else ""),
            "paths": readable,
            "generation_count": len(paths),
            "errors": errors,
            "recording_error": getattr(self, "_native_journal_error", None),
        }

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
        """Reap our child and close its pipes; caller owns the session lock."""
        self._native_journal_path = None
        proc = self._proc
        if proc is not None:
            if proc.poll() is None:
                proc.kill()
            proc.wait(timeout=10.0)
            _close_proc_streams(proc)
        self._close_private_job()
        self._proc = None
        self._ready_info = None
        self._owned = False


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
        if sys.platform == "win32":
            self._job_handle = _assign_kill_on_close_job(self._proc.pid)
        self._start_stderr_drain()

        line = self._read_stdio_line(timeout_s=CUBIT_READY_TIMEOUT_S)
        try:
            ready = json.loads(line.decode("utf-8"))
        except (ValueError, UnicodeError) as e:
            raise CubitSessionError(
                f"Daemon emitted non-JSON ready line: {line!r} ({e})")
        if not isinstance(ready, dict) or not ready.get("ready") or ready.get("protocol_version") != PROTOCOL_VERSION:
            raise CubitSessionError(
                f"Daemon reported startup failure: {ready}")
        self._ready_info = ready
        return ready

    def _call_via_stdio(self, req: dict, timeout_s: float) -> dict:
        self._send_stdio(req)
        line = self._read_stdio_line(timeout_s=timeout_s)
        try:
            resp = json.loads(line.decode("utf-8"))
        except (ValueError, UnicodeError) as exc:
            raise CubitSessionError(f"Invalid daemon response: {exc}") from exc
        if not isinstance(resp, dict) or resp.get("id") != req["id"]:
            raise CubitSessionError(
                f"Invalid response for request {req['id']}: {resp!r}")
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

        pipe = self._proc.stdout

        def _reader():
            try:
                line = pipe.readline()
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
        """Send shutdown without waiting for a response; tolerate a closed pipe."""
        try:
            req_id = self._next_id
            self._next_id += 1
            self._send_stdio({"id": req_id, "op": "shutdown", "args": [],
                              "protocol_version": PROTOCOL_VERSION})
        except (CubitSessionError, AssertionError):
            pass

    # ---- singleton access ----

    @classmethod
    def get(cls, mode: str = "batch") -> "CubitSession":
        """Return a mode-consistent process-wide singleton.

        A caller requesting headless execution must never inherit a GUI
        singleton created elsewhere in the process.
        """
        global _SINGLETON
        with _SESSION_LOCK:
            if _SINGLETON is None:
                _SINGLETON = cls(mode=mode)
            elif _SINGLETON._mode != mode:
                raise CubitSessionError(
                    f"Cubit singleton already uses mode={_SINGLETON._mode!r}; "
                    f"refusing requested mode={mode!r}. Use "
                    "cubit_session_shutdown to reset the MCP-owned session "
                    "before another headless session."
                )
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
