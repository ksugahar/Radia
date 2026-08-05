"""
bootstrap.py — Runs INSIDE the Cubit GUI process (Plan A).

Launched by `coreform_cubit.exe -nojournal bootstrap.py`. Cubit's
embedded Python executes this script, at which point:

    - Cubit's Qt event loop is already running (GUI is up)
    - PySide6 is available via Cubit's bundled site-packages
      (Cubit 2025.12 ships PySide6 / Qt6; linked against the SAME QtCore
      Cubit uses -> no ABI split.  Earlier Cubit shipped PyQt5/Qt5.)
    - `cubit` module is already loaded / initialized

We install a QTimer (200 ms) that polls a drop directory for JSON-RPC
request files, dispatches them on the Qt main thread, and writes
JSON responses into a sibling `out/` directory. Client (session.py)
drops a `<id>.req` file and polls for `out/<id>.resp`.

This design completely avoids sockets / named pipes / TCP (lab aesthetic
policy: no COMSOL-pattern), runs commands on the Qt main thread (safe),
and survives an arbitrary number of client requests for the lifetime of
the Cubit GUI window.

Protocol (mirrors daemon.py's stdio JSON-RPC):
    request file  `<drop>/<id>.req`:
        {"id": int, "op": str, "args": list}
    response file `<drop>/out/<id>.resp`:
        {"id": int, "ok": bool, "result": any}  OR
        {"id": int, "ok": false, "error": str}

Ready signal: `<drop>/ready` is written once QTimer is installed and
the bootstrap is operational. Client waits on this before first send.

Protocol version is echoed in the ready file for mutual check.

Environment:
    CUBIT_DROP_DIR — absolute path to the drop directory. Required.
                      Created if missing. `out/` is created too.
"""
from __future__ import annotations

import json
import os
import pathlib
import sys
import traceback


def _write_startup_error(text: str) -> None:
    """Persist a startup failure where session.py's ready-poll can see it.

    MathWorks MATLAB MCP pattern (mcp_startup_error.txt): the launching
    client polls for this file alongside the ready marker, so a license
    denial / Qt failure / import error surfaces as the REAL in-process
    error instead of a generic ready-timeout.  Uses only the env var so
    it works even before DROP is resolved.
    """
    drop = os.environ.get("CUBIT_DROP_DIR")
    if not drop:
        return
    try:
        p = pathlib.Path(drop)
        p.mkdir(parents=True, exist_ok=True)
        (p / "startup_error.txt").write_text(text, encoding="utf-8")
    except OSError:
        pass


try:
    from PySide6.QtCore import QTimer

    import cubit  # Cubit-bundled module; present in the GUI's Python
except Exception:
    _write_startup_error(traceback.format_exc())
    raise


PROTOCOL_VERSION = 2  # file-drop variant of daemon.py's protocol 1

_HERE = os.path.dirname(os.path.abspath(__file__))


def _read_drop_dir() -> pathlib.Path:
    drop = os.environ.get("CUBIT_DROP_DIR")
    if not drop:
        # Fail loud: we don't want to silently squat on a hardcoded path.
        raise RuntimeError(
            "CUBIT_DROP_DIR env var not set — cannot locate drop directory. "
            "session.py should set this before launch."
        )
    p = pathlib.Path(drop)
    p.mkdir(parents=True, exist_ok=True)
    (p / "out").mkdir(exist_ok=True)
    return p


try:
    DROP = _read_drop_dir()
except Exception:
    _write_startup_error(traceback.format_exc())
    raise
OUTBOX = DROP / "out"
LOG = DROP / "bootstrap.log"

# QTimer GC reference holder (QTimer is GC-collected if no strong ref is kept).
_holder: dict = {}


def _log(msg: str) -> None:
    try:
        with LOG.open("a", encoding="utf-8") as f:
            f.write(msg + "\n")
    except Exception:
        pass


def _write_response_atomic(stem: str, payload: dict) -> None:
    """Write response JSON atomically so clients never see a half-written file.

    Write to `<stem>.resp.tmp`, then os.rename to `<stem>.resp`. Rename
    is atomic on NTFS and POSIX.
    """
    tmp = OUTBOX / f"{stem}.resp.tmp"
    dst = OUTBOX / f"{stem}.resp"
    tmp.write_text(json.dumps(payload, default=str), encoding="utf-8")
    os.replace(tmp, dst)


def _cubit_error_count() -> int:
    """Return the cumulative Cubit error count, 0 if API unavailable.

    Cubit's `cubit.cmd()` returns rc=1 (truthy) even for commands that
    print ERROR ("No volume with ID 99" etc.) — the command is accepted
    but silently no-ops. `get_error_count()` (available in modern Cubit
    Python API) increments on each printed ERROR, so we sample the
    counter around each cmd to detect these silent failures.
    """
    try:
        return int(cubit.get_error_count())
    except Exception:
        return 0


def _cubit_last_error() -> str:
    try:
        msg = cubit.get_last_error()
        return str(msg) if msg else ""
    except Exception:
        return ""


def _op_cmd(args: list) -> list:
    """Execute a list of cubit commands. Mirror cubit_daemon._op_cmd.

    A command is considered failed if:
      - `cubit.cmd()` returns falsy rc, OR
      - Cubit's error count increased during the command (silent error).
    """
    results = []
    for line in args:
        err_before = _cubit_error_count()
        try:
            rc = cubit.cmd(str(line))
        except Exception:
            results.append({"line": line, "ok": False,
                            "error": traceback.format_exc(),
                            "error_count_delta": 0})
            break
        err_after = _cubit_error_count()
        delta = err_after - err_before
        rc_ok = bool(rc)
        ok = rc_ok and delta == 0
        step = {
            "line": line,
            "ok": ok,
            "rc": int(rc) if not isinstance(rc, bool) else int(rc),
            "error_count_delta": int(delta),
        }
        if delta > 0:
            msg = _cubit_last_error()
            if msg:
                step["error"] = msg
            else:
                step["error"] = f"Cubit emitted {delta} ERROR line(s) during this command"
        results.append(step)
        if not ok:
            break
    return results


def _op_probe(args: list) -> object:
    """Delegate to the SHARED probe implementation (probe_ops.op_probe)
    used by both the GUI file-drop runner (this file) and the batch
    stdio runner (daemon.py) -- single source, no query drift."""
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)
    import probe_ops
    return probe_ops.op_probe(cubit, args)


def _op_snapshot(args: list) -> dict:
    """Hardcopy current view to PNG. Mirror cubit_daemon._op_snapshot."""
    if len(args) < 1:
        return {"error": "snapshot requires at least 1 arg (path)"}
    out_path = str(args[0])
    w = int(args[1]) if len(args) >= 2 else 800
    h = int(args[2]) if len(args) >= 3 else 600
    cmd_str = f'hardcopy "{out_path}" png window {w} {h}'
    rc = cubit.cmd(cmd_str)
    return {"path": out_path, "ok": bool(rc), "width": w, "height": h}


def _dispatch(req: dict) -> dict:
    """Map one JSON-RPC request dict to a response dict."""
    req_id = req.get("id")
    op = req.get("op", "")
    args = req.get("args", [])
    if op == "ping":
        return {"id": req_id, "ok": True, "result": "pong"}
    if op == "shutdown":
        # Best-effort: ask Cubit to close itself. The QTimer will be
        # torn down implicitly when the Qt event loop exits.
        cubit.cmd("exit")
        return {"id": req_id, "ok": True, "result": "bye"}
    if op == "cmd":
        try:
            return {"id": req_id, "ok": True, "result": _op_cmd(args)}
        except Exception:
            return {"id": req_id, "ok": False, "error": traceback.format_exc()}
    if op == "probe":
        try:
            result = _op_probe(args)
            ok = not (isinstance(result, dict) and "error" in result)
            return {"id": req_id, "ok": ok, "result": result}
        except Exception:
            return {"id": req_id, "ok": False, "error": traceback.format_exc()}
    if op == "snapshot":
        try:
            result = _op_snapshot(args)
            return {"id": req_id, "ok": bool(result.get("ok")), "result": result}
        except Exception:
            return {"id": req_id, "ok": False, "error": traceback.format_exc()}
    return {"id": req_id, "ok": False, "error": f"unknown op: {op!r}"}


def _poll() -> None:
    try:
        reqs = sorted(DROP.glob("*.req"))
    except Exception as e:
        _log(f"glob error: {e}")
        return
    for req_path in reqs:
        stem = req_path.stem
        try:
            content = req_path.read_text(encoding="utf-8")
            req = json.loads(content)
        except Exception:
            _write_response_atomic(stem, {
                "ok": False,
                "error": f"malformed request file: {traceback.format_exc()}",
            })
            req_path.unlink(missing_ok=True)
            continue
        # Dispatch on main thread (QTimer callback is already main thread).
        resp = _dispatch(req)
        try:
            _write_response_atomic(stem, resp)
        except Exception:
            _log(f"failed to write response for {stem}: {traceback.format_exc()}")
        # Consume the request AFTER response is fully materialized.
        req_path.unlink(missing_ok=True)


def _install_ready_marker() -> None:
    ready = DROP / "ready"
    ready.write_text(
        json.dumps({
            "ready": True,
            "protocol_version": PROTOCOL_VERSION,
            "pid": os.getpid(),
            "drop": str(DROP),
        }),
        encoding="utf-8",
    )


def _start() -> None:
    t = QTimer()
    t.timeout.connect(_poll)
    t.start(200)  # 200 ms poll interval — good tradeoff of latency vs CPU
    _holder["timer"] = t
    _log(f"QTimer started, polling {DROP} every 200ms, proto v{PROTOCOL_VERSION}")
    _install_ready_marker()


try:
    _start()
except Exception:
    _write_startup_error(traceback.format_exc())
    raise
