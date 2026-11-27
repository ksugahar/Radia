"""Run the curved ``.vol`` export worker in an isolated external Python process.

The MCP server never imports Cubit or Netgen to export.  It starts
:mod:`radia_mcp.cubit.vol_export` with an explicit interpreter and the
isolated flags ``-I -B``, exchanges paths and a JSON status file, and fails
loudly when Cubit, the worker interpreter or the native plugin is missing.

Interpreter: the ``python`` argument, else ``$RADIA_CUBIT_WORKER_PYTHON``, else
the server's own interpreter; it must import the fork ``ngsolve``/``netgen``
the plugin was built against.  The worker modules are stdlib-only; the runner
puts the server's own ``radia_mcp`` source root first on the isolated worker's
``sys.path`` (``-I`` ignores ``PYTHONPATH``), so worker and server code match.  Plugin: the
``plugin_dir`` argument, else ``$RADIA_CUBIT_PLUGIN_DIR``; it must hold
``radia_cubit_export.dll`` (built from ``packages/radia-mcp/native/cubit_vol_export``).
"""
from __future__ import annotations

import json
import logging
import math
import os
import hashlib
import shutil
import subprocess
import uuid
import sys
from pathlib import Path
from typing import Any

from radia_mcp.cubit import kelvin as _kelvin
from radia_mcp.cubit.vol_export import PLUGIN_DLL, periodic_arguments

WORKER_MODULE = "radia_mcp.cubit.vol_export"
ISOLATED_FLAGS = ("-I", "-B")

_LOG = logging.getLogger(__name__)

SOURCE_ROOT = str(Path(__file__).resolve().parents[2])
_BOOT = ("import sys, runpy; sys.path.insert(0, %r); sys.argv = [%r] + sys.argv[1:]; "
         "runpy.run_module(%r, run_name='__main__', alter_sys=True)")

_PROBE = (
    "import json, sys\n"
    "sys.path.insert(0, %r)\n"
    "import importlib.util as u\n"
    "spec = u.find_spec('radia_mcp.cubit.vol_export')\n"
    "out = {'worker': spec.origin if spec else None}\n"
    "try:\n"
    "    import ngsolve, netgen\n"
    "    out['netgen'] = getattr(netgen, '__version__', None)\n"
    "    out['ngsolve'] = getattr(ngsolve, '__version__', None)\n"
    "except Exception as e:\n"
    "    out['netgen_error'] = f'{type(e).__name__}: {e}'\n"
    "print(json.dumps(out))\n"
) % SOURCE_ROOT


class WorkerRuntimeError(RuntimeError):
    """The worker interpreter, Cubit or the native plugin is not usable."""


class WorkerProtocolError(RuntimeError):
    """The worker wrote a status file that is not a JSON object."""


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


def resolve_plugin_dir(plugin_dir: str | os.PathLike[str] | None = None) -> Path:
    """Return the directory holding radia_cubit_export.dll; fail loudly when absent."""
    chosen = plugin_dir or os.environ.get("RADIA_CUBIT_PLUGIN_DIR")
    if not chosen:
        raise WorkerRuntimeError("Native plugin directory not set; pass plugin_dir or set "
                                 "RADIA_CUBIT_PLUGIN_DIR (build it with native/cubit_vol_export/build.py)")
    path = Path(chosen).resolve()
    if not (path / PLUGIN_DLL).is_file():
        raise WorkerRuntimeError(f"{PLUGIN_DLL} not found in {path}")
    return path


def probe_worker(python: str | os.PathLike[str] | None = None, *, timeout_s: float = 120.0) -> dict[str, Any]:
    """Check, in the isolated worker interpreter, that the worker and ngsolve/netgen import."""
    exe = resolve_python(python)
    run = subprocess.run([str(exe), *ISOLATED_FLAGS, "-c", _PROBE], stdin=subprocess.DEVNULL,
                         capture_output=True, timeout=timeout_s, creationflags=_creationflags())
    if run.returncode != 0:
        raise WorkerRuntimeError(f"Worker probe failed (rc {run.returncode}): "
                                 f"{run.stderr.decode(errors='replace')[-1000:]}")
    report = json.loads(run.stdout.decode().strip().splitlines()[-1])
    report["python"] = str(exe)
    if not report.get("worker"):
        raise WorkerRuntimeError(f"{WORKER_MODULE} is not importable from {SOURCE_ROOT} by {exe}")
    if report.get("netgen_error") or not report.get("netgen"):
        raise WorkerRuntimeError(f"ngsolve/netgen are not importable by {exe}: {report.get('netgen_error')}")
    return report


def run_export(journal: str | os.PathLike[str], out_path: str | os.PathLike[str], order: int, *,
               cubit_bin: str | os.PathLike[str] | None = None,
               python: str | os.PathLike[str] | None = None,
               periodic: dict[str, Any] | None = None,
               timeout_s: float = 1800.0,
               plugin_dir: str | os.PathLike[str] | None = None,
               kelvin: dict[str, Any] | None = None) -> dict[str, Any]:
    """Export ``journal`` to a curved ``.vol`` through the isolated worker.

    Returns the worker's status JSON (companion record and sha256 included)
    with the probe and command line added.  Raises
    :class:`WorkerRuntimeError` for a missing runtime; a failed export is
    returned as ``status == "error"`` with the worker's own message.
    """
    if int(order) != order or order < 1:
        raise ValueError(f"order must be a positive integer, got {order!r}")
    if not (isinstance(timeout_s, (int, float)) and math.isfinite(timeout_s) and timeout_s > 0):
        raise ValueError("timeout_s must be finite and positive")
    periodic_arguments(periodic)  # all nine fields or none; raises otherwise
    kelvin = dict(kelvin or {})
    _kelvin.validate_options(kelvin)
    journal_path = Path(journal).resolve()
    if not journal_path.is_file():
        raise FileNotFoundError(f"journal not found: {journal_path}")
    out = Path(out_path).resolve()
    if out.suffix != ".vol":
        raise ValueError(f"output must end in .vol: {out}")
    companion_path = Path(str(out) + ".json")
    for existing in (out, companion_path):  # early refusal; _claim re-checks under the lock
        if existing.exists():
            raise FileExistsError(f"refusing to overwrite an existing export output: {existing}")
    bin_dir = resolve_cubit_bin(cubit_bin)
    plugin = resolve_plugin_dir(plugin_dir)
    probe = probe_worker(python)
    lock = _claim(out)
    try:
        result = _export_locked(journal_path, out, companion_path, order, bin_dir, plugin, probe,
                                periodic, kelvin, timeout_s)
    except BaseException as exc:
        release_error = _release(lock)
        if release_error:
            exc.add_note(release_error)
        raise
    # The outputs are already published and verified; a lock that cannot be
    # removed is reported with the result, never by discarding that result.
    release_error = _release(lock)
    if release_error:
        result["lock_release_error"] = release_error
    return result


def _release(lock: Path) -> str | None:
    """Remove the claim lock; return a description of the failure, or None."""
    try:
        lock.unlink()
    except FileNotFoundError:
        return None  # already gone: nothing of ours is left to release
    except OSError as error:
        return (f"could not remove export lock {lock} ({type(error).__name__}: {error}); "
                "remove it by hand before exporting to this path again")
    return None


def _claim(out: Path) -> Path:
    """Claim ``out`` with an O_EXCL lock file held through publication; refuse existing outputs."""
    lock = Path(str(out) + ".lock")
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise FileExistsError(f"another export holds {lock}") from None
    # From here the lock is ours: any failure releases it and keeps the
    # primary error, with a failed release attached as a note naming the lock.
    try:
        os.close(descriptor)
        for existing in (out, Path(str(out) + ".json")):
            if existing.exists():
                raise FileExistsError(f"refusing to overwrite an existing export output: {existing}")
    except BaseException as exc:
        release_error = _release(lock)
        if release_error:
            exc.add_note(release_error)
        raise
    return lock


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _export_locked(journal_path, out, companion_path, order, bin_dir, plugin, probe,
                   periodic, kelvin, timeout_s) -> dict[str, Any]:
    token = uuid.uuid4().hex
    stage = out.parent / f".{out.name}.{token}.stage"
    stage.mkdir()
    status: dict[str, Any] | None = None
    try:  # every step after the stage exists runs under the stage-cleanup guard
        staged = stage / out.name
        staged_companion = Path(str(staged) + ".json")
        status_path = stage / "status.json"
        boot = _BOOT % (SOURCE_ROOT, WORKER_MODULE, WORKER_MODULE)
        command = [probe["python"], *ISOLATED_FLAGS, "-c", boot, "--cubit-bin", str(bin_dir),
                   "--plugin-dir", str(plugin), "--journal", str(journal_path), "--out", str(staged),
                   "--order", str(int(order)), "--result", str(status_path),
                   "--token", token]
        if periodic:
            command += ["--periodic", json.dumps(periodic)]
        if kelvin.get("add_kelvin"):
            command += ["--kelvin", json.dumps(kelvin)]
        # cwd stays the destination's parent so relative paths in the journal
        # resolve as before; every staged output is passed as an absolute path.
        run = subprocess.run(command, stdin=subprocess.DEVNULL, capture_output=True, cwd=str(out.parent),
                             timeout=timeout_s, creationflags=_creationflags())
        if status_path.is_file():
            status = _read_status(status_path, run.returncode)
        else:
            status = {"status": "error", "error_type": "WorkerCrash",
                      "error": f"worker exited with {run.returncode} and wrote no status file"}
        status.update({"returncode": run.returncode, "command": command, "worker_probe": probe,
                       "stdout_tail": run.stdout.decode(errors="replace")[-2000:],
                       "stderr_tail": run.stderr.decode(errors="replace")[-2000:]})
        if status.get("status") == "ok":
            problem = _verify_outputs(status, token, staged, staged_companion, run.returncode)
            if not problem:
                problem = _publish(status, token, ((staged, out), (staged_companion, companion_path)))
                if not problem:
                    status.update({"vol": str(out), "companion_json": str(companion_path)})
            if problem:
                status.update({"status": "error", "error_type": "WorkerOutputMismatch", "error": problem})
    except BaseException as exc:
        _cleanup_stage(stage, exc, status)
        raise
    _cleanup_stage(stage, None, status)
    return status


def _read_status(status_path: Path, returncode: int) -> dict[str, Any]:
    """Parse the worker status file; anything but a JSON object is a protocol error."""
    try:
        status = json.loads(status_path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise WorkerProtocolError(f"worker (rc {returncode}) wrote malformed status JSON to "
                                  f"{status_path}: {error}") from error
    if not isinstance(status, dict):
        raise WorkerProtocolError(f"worker (rc {returncode}) status in {status_path} is a "
                                  f"{type(status).__name__}, not a JSON object")
    if status.get("status") not in ("ok", "error"):
        raise WorkerProtocolError(f"worker (rc {returncode}) status in {status_path} has "
                                  f"status {status.get('status')!r}; expected 'ok' or 'error'")
    return status


def _cleanup_stage(stage: Path, exc: BaseException | None, status: dict[str, Any] | None) -> None:
    """Remove the staging directory best-effort; a failure never replaces ``exc``.

    A cleanup failure is logged; it is attached to ``exc`` as a note when an
    error is propagating, else recorded in ``status["stage_cleanup_error"]``.
    """
    try:
        shutil.rmtree(stage)
    except OSError as cleanup:
        message = f"staging cleanup failed for {stage}: {type(cleanup).__name__}: {cleanup}"
        _LOG.error(message)
        if exc is not None:
            exc.add_note(message)
        elif status is not None:
            status["stage_cleanup_error"] = message


def _publish(status: dict[str, Any], token: str, pairs: tuple[tuple[Path, Path], ...]) -> str | None:
    """Hard-link each staged file to its destination, all or nothing.

    Hard links fail when the destination exists, so publication never replaces
    bytes another run returned.  On any failure after the first link, only the
    destinations this run linked (same file as the staged source, confirmed
    before unlinking) are removed, and the original error is re-raised.
    Returns a verification problem (after rolling back) or None.
    """
    published: list[tuple[Path, Path]] = []
    try:
        for staged, dest in pairs:
            os.link(staged, dest)
            published.append((staged, dest))
        problem = _verify_outputs(status, token, pairs[0][1], pairs[1][1], 0)
    except BaseException as exc:
        try:
            _rollback(published)
        except OSError as rollback:
            exc.add_note(f"rollback of published outputs failed: {type(rollback).__name__}: {rollback}")
        raise
    if problem:
        try:
            _rollback(published)
        except OSError as rollback:
            problem += f"; rollback: {rollback}"
    return problem


def _rollback(published: list[tuple[Path, Path]]) -> None:
    """Unlink destinations this run published, only when they are still the staged file.

    A destination replaced by someone else is left untouched, and the rest of
    this run's files are still removed; the refusals are raised at the end.
    """
    problems = []
    for staged, dest in reversed(published):
        try:
            if dest.exists() and os.path.samefile(staged, dest):
                dest.unlink()
            elif dest.exists():
                problems.append(f"{dest} is no longer this run's output (left in place)")
        except OSError as error:
            problems.append(f"{dest}: {type(error).__name__}: {error}")
    if problems:
        raise OSError("rollback incomplete: " + "; ".join(problems))


def _verify_outputs(status: dict[str, Any], token: str, out: Path, companion: Path,
                    returncode: int) -> str | None:
    """Return why this run's ok status is not backed by its own outputs, or None."""
    if returncode != 0:
        return f"worker exited with {returncode} but reported ok"
    if status.get("token") != token:
        return f"status token {status.get('token')!r} does not match this run {token!r}"
    if not out.is_file() or not companion.is_file():
        return f"worker reported ok but {out} or {companion} is missing"
    for path, key in ((out, "sha256"), (companion, "companion_sha256")):
        digest = _sha256(path)
        if status.get(key) != digest:
            return f"worker {key} {status.get(key)!r} does not match {path} ({digest})"
    try:
        record = json.loads(companion.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        return f"{companion} is not valid JSON: {error}"
    if record != status.get("companion"):
        return f"{companion} does not match the companion record the worker reported"
    return None
