"""Wrap long-running calc_*.py scripts with AsyncRunner + JSON stdout parsing.

Bridges the lab's existing calc_*.py panel scripts (which run via
`subprocess.Popen` and emit JSON on stdout, per Layer 4 architecture
in CLAUDE.md) with the AsyncRunner pattern.

Usage in an MCP tool:

    from radia_mcp.common.subprocess_runner import run_calc_async, \
        get_run_status, cancel_run, list_runs

    @mcp.tool()
    def motor_solve_start(vol_path: str, n_steps: int = 1000) -> dict:
        return run_calc_async(
            script="calc_motor_transient.py",
            args=["--vol", vol_path, "--n-steps", str(n_steps)],
            label=f"motor_transient({vol_path})",
        )

    @mcp.tool()
    def motor_solve_status(run_id: str) -> dict:
        return get_run_status(run_id)

    @mcp.tool()
    def motor_solve_cancel(run_id: str) -> dict:
        return cancel_run(run_id)

Process is killed on cancel (kill_tree=True on Windows to also kill
NGSolve subprocesses).
"""

from __future__ import annotations
import json
import os
import subprocess
import sys
import threading
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Optional


# Process tracking: {run_id: RunHandle}
_active_runs: dict[str, "RunHandle"] = {}
_runs_lock = threading.Lock()


class RunHandle:
    """A single subprocess + its capture / status."""

    def __init__(self, run_id: str, script: str, args: list[str],
                 label: Optional[str]):
        self.run_id = run_id
        self.script = script
        self.args = args
        self.label = label or f"{Path(script).name}({' '.join(args)})"
        self.start_time = datetime.now()
        self.end_time: Optional[datetime] = None
        self.process: Optional[subprocess.Popen] = None
        self.stdout_lines: list[str] = []
        self.stderr_lines: list[str] = []
        self.parsed_json: Optional[Any] = None
        self.exit_code: Optional[int] = None
        self.status: str = "starting"   # starting / running / completed / failed / cancelled
        self._reader_threads: list[threading.Thread] = []
        self._lock = threading.Lock()

    def _reader(self, stream, sink: list[str]) -> None:
        """Background thread: drain stream into sink line by line."""
        try:
            for line in iter(stream.readline, ""):
                if not line:
                    break
                with self._lock:
                    sink.append(line.rstrip("\n"))
        finally:
            try:
                stream.close()
            except Exception:
                pass

    def start(self, cwd: Optional[str | Path] = None) -> None:
        """Spawn the subprocess + reader threads."""
        cmd = [sys.executable, self.script] + list(self.args)
        env = os.environ.copy()
        env.setdefault("PYTHONUNBUFFERED", "1")
        try:
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                stdin=subprocess.DEVNULL,
                text=True,
                bufsize=1,           # line-buffered
                cwd=str(cwd) if cwd else None,
                env=env,
                # CREATE_NEW_PROCESS_GROUP on Windows -> CTRL+BREAK works
                creationflags=(subprocess.CREATE_NEW_PROCESS_GROUP
                               if os.name == "nt" else 0),
            )
        except FileNotFoundError as e:
            self.status = "failed"
            self.stderr_lines.append(f"Script not found: {e}")
            self.end_time = datetime.now()
            self.exit_code = -1
            return
        self.status = "running"
        # Reader threads
        for stream, sink in [(self.process.stdout, self.stdout_lines),
                              (self.process.stderr, self.stderr_lines)]:
            t = threading.Thread(target=self._reader,
                                 args=(stream, sink), daemon=True)
            t.start()
            self._reader_threads.append(t)

        # Status-update thread
        def _watch():
            self.process.wait()
            for t in self._reader_threads:
                t.join(timeout=2.0)
            with self._lock:
                self.end_time = datetime.now()
                self.exit_code = self.process.returncode
                if self.status == "running":   # not externally cancelled
                    self.status = "completed" if self.exit_code == 0 else "failed"
                # Try to parse JSON from the LAST line of stdout
                if self.stdout_lines:
                    for line in reversed(self.stdout_lines):
                        s = line.strip()
                        if not s:
                            continue
                        if s.startswith("{") or s.startswith("["):
                            try:
                                self.parsed_json = json.loads(s)
                                break
                            except json.JSONDecodeError:
                                pass
                        # Stop after the first non-empty line if it's not JSON
                        break

        threading.Thread(target=_watch, daemon=True).start()

    def cancel(self) -> bool:
        """Kill the process (and its tree on Windows)."""
        if self.process is None:
            return False
        if self.process.poll() is not None:
            return False    # already done
        try:
            if os.name == "nt":
                # taskkill /T kills the whole process tree
                subprocess.run(
                    ["taskkill", "/F", "/T", "/PID", str(self.process.pid)],
                    capture_output=True, check=False)
            else:
                self.process.terminate()
            with self._lock:
                self.status = "cancelled"
            return True
        except Exception:
            return False

    def status_dict(self) -> dict:
        with self._lock:
            elapsed = ((self.end_time or datetime.now()) -
                        self.start_time).total_seconds()
            return {
                "run_id": self.run_id,
                "label": self.label,
                "status": self.status,
                "start_time": self.start_time.isoformat(),
                "end_time": (self.end_time.isoformat()
                              if self.end_time else None),
                "elapsed_seconds": round(elapsed, 1),
                "exit_code": self.exit_code,
                "stdout_lines": len(self.stdout_lines),
                "stderr_lines": len(self.stderr_lines),
                "last_stdout": (self.stdout_lines[-1]
                                  if self.stdout_lines else None),
                "last_stderr": (self.stderr_lines[-1]
                                  if self.stderr_lines else None),
                "result_json": self.parsed_json,
            }


def run_calc_async(
    script: str,
    args: Optional[list[str]] = None,
    label: Optional[str] = None,
    cwd: Optional[str | Path] = None,
) -> dict:
    """Start a calc_*.py subprocess in the background.

    Args:
        script: path to the script (absolute, or relative to cwd)
        args: command-line arguments (list of strings)
        label: user-friendly name (default: script(args))
        cwd: working directory

    Returns:
        {"run_id": "...", "label": "...", "status": "running" / "failed"}
    """
    run_id = str(uuid.uuid4())[:8]
    handle = RunHandle(run_id, script, args or [], label)
    handle.start(cwd=cwd)
    with _runs_lock:
        _active_runs[run_id] = handle
    return handle.status_dict()


def get_run_status(run_id: str) -> dict:
    """Get current status of a running / finished subprocess."""
    with _runs_lock:
        handle = _active_runs.get(run_id)
    if handle is None:
        return {"error": f"run_id '{run_id}' not found",
                "known_runs": list_runs()}
    return handle.status_dict()


def cancel_run(run_id: str) -> dict:
    """Kill the subprocess (tree on Windows)."""
    with _runs_lock:
        handle = _active_runs.get(run_id)
    if handle is None:
        return {"error": f"run_id '{run_id}' not found"}
    cancelled = handle.cancel()
    return {"run_id": run_id, "cancelled": cancelled,
            "status": handle.status_dict()}


def list_runs() -> list[dict]:
    """Summary of all known runs (running + completed)."""
    with _runs_lock:
        return [{"run_id": h.run_id, "label": h.label, "status": h.status}
                for h in _active_runs.values()]


def get_run_output(run_id: str, tail: int = 50) -> dict:
    """Get the captured stdout / stderr of a run.

    Args:
        run_id: id from run_calc_async
        tail: return at most this many lines from each stream (default 50)
    """
    with _runs_lock:
        handle = _active_runs.get(run_id)
    if handle is None:
        return {"error": f"run_id '{run_id}' not found"}
    with handle._lock:
        return {
            "run_id": run_id,
            "label": handle.label,
            "status": handle.status,
            "exit_code": handle.exit_code,
            "stdout": handle.stdout_lines[-tail:],
            "stderr": handle.stderr_lines[-tail:],
            "result_json": handle.parsed_json,
        }


def clear_finished_runs() -> dict:
    """Forget completed / failed / cancelled runs (free memory)."""
    with _runs_lock:
        keep = {rid: h for rid, h in _active_runs.items()
                if h.status in ("starting", "running")}
        removed = len(_active_runs) - len(keep)
        _active_runs.clear()
        _active_runs.update(keep)
    return {"removed": removed, "remaining": len(keep)}
