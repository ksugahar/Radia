"""Run the curved ``.vol`` worker in an isolated external Python process.

The MCP server never imports Cubit or Netgen to export.  It starts
:mod:`radia_mcp.cubit.vol_export` with an explicit interpreter and the
isolated flags ``-I -B``, exchanges paths and a JSON status file, and fails
loudly when Cubit or the callback-geometry Netgen build is missing.  The
retired native ``export netgen`` command plugin is never used.

Interpreter selection: the ``python`` argument, else ``$RADIA_CUBIT_WORKER_PYTHON``,
else the server's own interpreter.  The interpreter must have ``radia_mcp`` and
a Netgen build that provides ``Mesh.GeometricRefit`` importable under ``-I``.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

from radia_mcp.cubit import periodic as _periodic

WORKER_MODULE = "radia_mcp.cubit.vol_export"
ISOLATED_FLAGS = ("-I", "-B")

_PROBE = (
    "import json, importlib.util as u\n"
    "spec = u.find_spec('radia_mcp.cubit.vol_export')\n"
    "out = {'worker': spec.origin if spec else None}\n"
    "try:\n"
    "    import netgen, netgen.meshing as m\n"
    "    out['netgen'] = getattr(netgen, '__version__', None)\n"
    "    out['geometric_refit'] = hasattr(m.Mesh, 'GeometricRefit')\n"
    "except Exception as e:\n"
    "    out['netgen_error'] = f'{type(e).__name__}: {e}'\n"
    "print(json.dumps(out))\n"
)


class WorkerRuntimeError(RuntimeError):
    """The worker interpreter, Cubit or the Netgen build is not usable."""


def _creationflags() -> int:
    return getattr(subprocess, "BELOW_NORMAL_PRIORITY_CLASS", 0)


def resolve_python(python: str | os.PathLike[str] | None = None) -> Path:
    """Return the explicit absolute worker interpreter; never search PATH."""
    chosen = python or os.environ.get("RADIA_CUBIT_WORKER_PYTHON") or sys.executable
    path = Path(chosen)
    if not path.is_absolute() or not path.is_file():
        raise WorkerRuntimeError(f"Worker Python must be an existing absolute path: {chosen!r}")
    return path


def resolve_cubit_bin(cubit_bin: str | os.PathLike[str] | None = None) -> Path:
    """Return a Cubit bin directory holding the Cubit Python package."""
    chosen = cubit_bin or os.environ.get("RADIA_CUBIT_BIN")
    if not chosen:
        from radia_mcp.cubit.session import find_cubit_install
        chosen = find_cubit_install()
    if not chosen:
        raise WorkerRuntimeError("Coreform Cubit not found; pass cubit_bin or set RADIA_CUBIT_BIN")
    path = Path(chosen)
    if not (path / "cubit").is_dir():
        raise WorkerRuntimeError(f"Cubit Python package not found under {path}")
    return path


def probe_worker(python: str | os.PathLike[str] | None = None, *, timeout_s: float = 120.0) -> dict[str, Any]:
    """Check, in the isolated worker interpreter, that the worker and GeometricRefit import."""
    exe = resolve_python(python)
    run = subprocess.run([str(exe), *ISOLATED_FLAGS, "-c", _PROBE], stdin=subprocess.DEVNULL,
                         capture_output=True, timeout=timeout_s, creationflags=_creationflags())
    if run.returncode != 0:
        raise WorkerRuntimeError(f"Worker probe failed (rc {run.returncode}): "
                                 f"{run.stderr.decode(errors='replace')[-1000:]}")
    report = json.loads(run.stdout.decode().strip().splitlines()[-1])
    report["python"] = str(exe)
    if not report.get("worker"):
        raise WorkerRuntimeError(f"{WORKER_MODULE} is not importable by {exe} under -I")
    if not report.get("geometric_refit"):
        raise WorkerRuntimeError(
            f"The Netgen build used by {exe} lacks Mesh.GeometricRefit "
            f"({report.get('netgen_error') or report.get('netgen')}); install the "
            "callback-geometry Netgen build in the worker environment")
    return report


def run_export(journal: str | os.PathLike[str], out_path: str | os.PathLike[str], order: int, *,
               cubit_bin: str | os.PathLike[str] | None = None,
               python: str | os.PathLike[str] | None = None,
               periodic: dict[str, Any] | None = None,
               timeout_s: float = 1800.0) -> dict[str, Any]:
    """Export ``journal`` to a curved ``.vol`` through the isolated worker.

    Returns the worker's status JSON with the probe and command line added.
    Raises :class:`WorkerRuntimeError` for a missing runtime; a failed export
    is returned as ``status == "error"`` with the worker's own message.
    """
    if int(order) != order or order < 1:
        raise ValueError(f"order must be a positive integer, got {order!r}")
    if not (isinstance(timeout_s, (int, float)) and math.isfinite(timeout_s) and timeout_s > 0):
        raise ValueError("timeout_s must be finite and positive")
    journal_path = Path(journal).resolve()
    if not journal_path.is_file():
        raise FileNotFoundError(f"journal not found: {journal_path}")
    out = Path(out_path).resolve()
    if out.suffix != ".vol":
        raise ValueError(f"output must end in .vol: {out}")
    if out.exists():
        raise FileExistsError(f"refusing to overwrite an existing .vol: {out}")
    bin_dir = resolve_cubit_bin(cubit_bin)
    probe = probe_worker(python)
    status_path = Path(str(out) + ".status.json")
    command = [probe["python"], *ISOLATED_FLAGS, "-m", WORKER_MODULE, "--cubit-bin", str(bin_dir),
               "--journal", str(journal_path), "--out", str(out), "--order", str(int(order)),
               "--result", str(status_path)]
    for name, value in (periodic or {}).items():
        if name not in _periodic.FIELDS:
            raise ValueError(f"unknown periodic field: {name}")
        if value is not None and value != "":
            command += ["--" + name.replace("_", "-"), str(value)]
    run = subprocess.run(command, stdin=subprocess.DEVNULL, capture_output=True, cwd=str(out.parent),
                         timeout=timeout_s, creationflags=_creationflags())
    if status_path.is_file():
        status = json.loads(status_path.read_text(encoding="utf-8"))
    else:
        status = {"status": "error", "error_type": "WorkerCrash",
                  "error": f"worker exited with {run.returncode} and wrote no status file"}
    status.update({"returncode": run.returncode, "command": command, "worker_probe": probe,
                   "status_json": str(status_path),
                   "stderr_tail": run.stderr.decode(errors="replace")[-2000:]})
    if run.returncode != 0 and status.get("status") == "ok":
        status["status"] = "error"
    return status
