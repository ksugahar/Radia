"""Threading-based async runner with progress + cancel.

Adapted from wjc9011/COMSOL_Multiphysics_MCP src/async_handler/solver.py
(which wraps long COMSOL solves). Generalized here for arbitrary
long-running commands (Cubit mesh + NGSolve solve + Radia C++ Solve
+ Optuna sweep + MCMC chain, etc.).

Why bother:
  - calc_motor_transient.py / calc_fem_kelvin.py / calc_motor_lamination.py
    routinely run for minutes to hours. A synchronous MCP call would
    block the entire conversation.
  - With AsyncRunner: kick off → poll progress every few seconds →
    cancel if needed.

Usage:
    from radia_mcp.common import AsyncRunner

    runner = AsyncRunner()

    # Start a long task. Returns immediately.
    runner.start(
        target=run_motor_transient,   # any callable
        args=(model_path,),
        kwargs={"n_steps": 1000},
        label="motor_transient_run",
    )

    # Poll status. Called from a separate MCP tool.
    status = runner.get_status()
    # → {"status": "running", "progress": 0.42, "elapsed_seconds": 87, ...}

    # Cancel (cooperative; the target must check `runner.cancel_requested`).
    runner.cancel()

    # Wait for completion (blocking).
    result = runner.wait(timeout=600.0)
    # → {"success": True, "result": ..., "elapsed_seconds": 132}
"""

from __future__ import annotations
import threading
import traceback
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Any, Callable, Optional


class RunStatus(Enum):
    IDLE = "idle"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class RunProgress:
    status: RunStatus = RunStatus.IDLE
    progress: float = 0.0          # 0.0 to 1.0
    message: str = ""
    label: Optional[str] = None    # user-friendly task name
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    result: Any = None             # filled on COMPLETED
    error: Optional[str] = None    # filled on FAILED
    traceback: Optional[str] = None

    def to_dict(self) -> dict:
        now = datetime.now()
        elapsed = ((self.end_time or now) - self.start_time).total_seconds() \
            if self.start_time else 0.0
        return {
            "status": self.status.value,
            "progress": self.progress,
            "message": self.message,
            "label": self.label,
            "start_time": self.start_time.isoformat() if self.start_time else None,
            "end_time": self.end_time.isoformat() if self.end_time else None,
            "elapsed_seconds": elapsed,
            "error": self.error,
        }


class AsyncRunner:
    """Run a single long-running callable in a background thread.

    Single-task model: only one task at a time per AsyncRunner instance.
    Use multiple AsyncRunner instances for parallel tasks.
    """

    def __init__(self):
        self._thread: Optional[threading.Thread] = None
        self._progress: RunProgress = RunProgress()
        self._cancel_flag: bool = False
        self._lock = threading.Lock()

    @property
    def is_running(self) -> bool:
        with self._lock:
            return self._progress.status == RunStatus.RUNNING

    @property
    def cancel_requested(self) -> bool:
        """True if the user has called cancel().

        Long-running targets SHOULD check this periodically and exit
        cleanly. The runner cannot force-kill the thread.
        """
        with self._lock:
            return self._cancel_flag

    def start(
        self,
        target: Callable[..., Any],
        args: tuple = (),
        kwargs: Optional[dict] = None,
        label: Optional[str] = None,
        progress_callback: Optional[Callable[[float, str], None]] = None,
    ) -> bool:
        """Start the target in a background thread.

        Args:
            target: callable to run
            args, kwargs: passed to target
            label: user-friendly task name (shown in status)
            progress_callback: optional `lambda progress, msg: ...`
                               for the TARGET to call (target receives
                               `_progress_cb=progress_callback` if it
                               accepts that kwarg).

        Returns:
            True if started, False if another task is already running.
        """
        kwargs = dict(kwargs or {})
        with self._lock:
            if self._progress.status == RunStatus.RUNNING:
                return False
            self._cancel_flag = False
            self._progress = RunProgress(
                status=RunStatus.RUNNING,
                progress=0.0,
                message="Starting...",
                label=label,
                start_time=datetime.now(),
            )

        def _update_progress(p: float, msg: str = "") -> None:
            """Closure passed to the target as _progress_cb."""
            with self._lock:
                self._progress.progress = max(0.0, min(1.0, p))
                if msg:
                    self._progress.message = msg
            if progress_callback:
                progress_callback(p, msg)

        def _thread_main():
            try:
                # If target accepts _progress_cb, inject it
                import inspect
                sig = inspect.signature(target)
                if "_progress_cb" in sig.parameters:
                    kwargs["_progress_cb"] = _update_progress
                if "_cancel_check" in sig.parameters:
                    kwargs["_cancel_check"] = lambda: self.cancel_requested

                result = target(*args, **kwargs)

                with self._lock:
                    if self._cancel_flag:
                        self._progress.status = RunStatus.CANCELLED
                        self._progress.message = "Cancelled."
                    else:
                        self._progress.status = RunStatus.COMPLETED
                        self._progress.progress = 1.0
                        self._progress.message = "Completed successfully."
                        self._progress.result = result
                    self._progress.end_time = datetime.now()
            except Exception as e:
                with self._lock:
                    self._progress.status = RunStatus.FAILED
                    self._progress.error = str(e)
                    self._progress.traceback = traceback.format_exc()
                    self._progress.message = f"Failed: {e}"
                    self._progress.end_time = datetime.now()

        self._thread = threading.Thread(target=_thread_main, daemon=True)
        self._thread.start()
        return True

    def cancel(self) -> bool:
        """Request cancellation (cooperative).

        Sets `cancel_requested = True`. The running target must check
        this and exit. If the target doesn't cooperate, the thread
        runs to completion.

        Returns:
            True if running task was signalled, False if no task.
        """
        with self._lock:
            if self._progress.status != RunStatus.RUNNING:
                return False
            self._cancel_flag = True
            self._progress.message = "Cancellation requested..."
            return True

    def wait(self, timeout: Optional[float] = None) -> dict:
        """Block until the running task finishes (or timeout).

        Returns the final status dict (with `result` filled if
        COMPLETED).
        """
        if self._thread is not None:
            self._thread.join(timeout=timeout)
        status = self.get_status()
        with self._lock:
            if self._progress.status == RunStatus.COMPLETED:
                status["result"] = self._progress.result
        return status

    def get_status(self) -> dict:
        """Current status dict (safe to call from any thread)."""
        with self._lock:
            return self._progress.to_dict()

    def reset(self) -> None:
        """Forget previous run. Errors if currently running."""
        with self._lock:
            if self._progress.status == RunStatus.RUNNING:
                raise RuntimeError("Cannot reset while task is running")
            self._progress = RunProgress()
            self._cancel_flag = False
