"""
daemon.py — JSON-RPC server exposing a persistent Cubit session.

Runs under the selected Cubit's bundled Python (Python 3.12 in Cubit 2026.8).
The compiled `cubit` binding must use that installation's Python ABI.

Role in the lab stack:
    LLM client
        → cubit-mcp (Python 3.12)
            → THIS DAEMON (Cubit Python, subprocess)
                → cubit.cmd(...) in batch/nographics mode

Design follows the OCP CAD Viewer philosophy transposed to Cubit:
a persistent backend, lightweight client calls, no per-action launch
cost. Latency target: <1ms per command after first-call init.

Protocol: line-delimited JSON-RPC over stdin/stdout (stderr is free
for logging and Cubit's own output).

    request:  {"id": int, "op": str, "args": list}
    response: {"id": int, "ok": bool, "result": any}  OR
              {"id": int, "ok": false, "error": str}

Supported `op` values:
    "cmd"       args=[str, ...]    execute cubit commands, return list of
                                    per-command success flags
    "probe"     args=[query]       return lab-common geometry queries
    "shutdown"  args=[]            graceful exit

On startup the daemon emits exactly one "ready" line before entering
the read loop, e.g. {"ready": true, "version": "..."}.
"""

import json
import os
import sys
import traceback

# Redirect Cubit's noisy stdout (banner, "Successfully created ..." etc.)
# to stderr BEFORE importing/initializing cubit. Otherwise those prints
# go to the same fd we use for JSON-RPC responses and corrupt the
# protocol stream. We preserve the original stdout fd for our own use.
_REAL_STDOUT_FD = os.dup(1)
os.dup2(2, 1)    # fd 1 -> fd 2: any C-level print to stdout now hits stderr

# Keep the binding lookup independent of package imports before Cubit startup.
_HERE = os.path.dirname(os.path.abspath(__file__))

def _load_layout():
    """The package's single BUILD_INFO reader, loaded by path (this runs as a script)."""
    import importlib.util
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "layout.py")
    spec = importlib.util.spec_from_file_location("_radia_mcp_cubit_layout", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _discover_cubit_bin():
    if override := os.environ.get("CUBIT_BIN_DIR"):
        cubit_bin = override
    else:
        root = os.environ.get("CUBIT_INSTALL_DIR") or os.path.join(
            os.environ.get("ProgramFiles", r"C:\Program Files"), "Coreform Cubit 2026.8"
        )
        cubit_bin = os.path.join(root, "bin")
    layout = _load_layout()
    root = os.path.dirname(os.path.normpath(cubit_bin))
    if not layout.is_supported_cubit_root(root):
        raise RuntimeError(
            f"The headless MCP runtime requires Coreform Cubit {layout.SUPPORTED_CUBIT_RELEASE}; "
            f"BUILD_INFO under {root} records {layout.cubit_release(root)!r}")
    return cubit_bin



def _write_response(obj):
    """Emit a single line of JSON to the preserved stdout fd and flush."""
    line = json.dumps(obj, default=str) + "\n"
    os.write(_REAL_STDOUT_FD, line.encode("utf-8"))


def _log(msg):
    """Debug log to stderr (mcp-server will inherit but not parse it)."""
    sys.stderr.write(f"[cubit_daemon] {msg}\n")
    sys.stderr.flush()


def _init_cubit():
    """Import cubit and initialize. Returns the imported module."""
    cubit_bin = _discover_cubit_bin()
    if cubit_bin not in sys.path:
        sys.path.insert(0, cubit_bin)
    import cubit  # noqa: E402  (deferred until CLI validation)

    # -nojournal: we don't need automatic .jou journaling.
    # -noecho:    don't re-echo commands to stdout (we'd corrupt JSON-RPC).
    argv = ["cubit", "-nojournal", "-noecho", "-batch", "-nographics"]
    cubit.init(argv)
    return cubit


def _cubit_error_count(cubit_mod):
    try:
        return int(cubit_mod.get_error_count())
    except Exception:
        return 0


def _cubit_last_error(cubit_mod):
    try:
        msg = cubit_mod.get_last_error()
        return str(msg) if msg else ""
    except Exception:
        return ""


def _op_cmd(cubit_mod, args):
    """Execute a list of cubit commands. Returns per-command success flags.

    Cubit's `cmd` returns rc=1 (truthy) even for commands that print
    ERROR ("No volume with ID 99" etc.) — accepted but silent no-op.
    We sample `get_error_count()` around each call so silent failures
    are reported as ok=False instead of slipping through the batch
    dry-run into the persistent headless session.
    """
    results = []
    for line in args:
        err_before = _cubit_error_count(cubit_mod)
        try:
            rc = cubit_mod.cmd(str(line))
        except Exception:
            results.append({"line": line, "ok": False,
                            "error": traceback.format_exc(),
                            "error_count_delta": 0})
            break
        err_after = _cubit_error_count(cubit_mod)
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
            msg = _cubit_last_error(cubit_mod)
            step["error"] = (msg if msg
                             else f"Cubit emitted {delta} ERROR line(s)")
        results.append(step)
        if not ok:
            break
    return results


def _op_probe(cubit_mod, args):
    """Delegate to the Cubit-side probe implementation."""
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)
    import probe_ops
    return probe_ops.op_probe(cubit_mod, args)


def main():
    """Entry point: initialize Cubit, emit ready, enter RPC loop."""
    requested_mode = os.environ.get("CUBIT_DAEMON_MODE", "batch").lower()
    try:
        if requested_mode != "batch":
            raise ValueError("Only headless CUBIT_DAEMON_MODE=batch is supported")
        cubit_mod = _init_cubit()
        version = cubit_mod.get_version()
        _write_response({
            "ready": True,
            "mode": "batch",
            "protocol_version": 1,
            "cubit_version": version,
        })
    except Exception:
        _write_response({"ready": False, "error": traceback.format_exc()})
        return 2

    # Read line-by-line from the raw bytes stream. Using
    # `for line in sys.stdin` with a piped parent blocks on internal
    # buffering and never releases individual lines until the writer
    # closes. `sys.stdin.buffer.readline()` yields lines as they arrive.
    stdin_bin = sys.stdin.buffer
    while True:
        line_bytes = stdin_bin.readline()
        if not line_bytes:
            break   # EOF
        line = line_bytes.decode("utf-8", errors="replace").strip()
        if not line:
            continue
        try:
            req = json.loads(line)
        except Exception:
            _write_response({"ok": False,
                             "error": f"invalid JSON: {line[:120]!r}"})
            continue

        req_id = req.get("id")
        op = req.get("op", "")
        args = req.get("args", [])

        if op == "shutdown":
            _write_response({"id": req_id, "ok": True, "result": "bye"})
            break
        elif op == "ping":
            _write_response({"id": req_id, "ok": True, "result": "pong"})
        elif op == "cmd":
            try:
                result = _op_cmd(cubit_mod, args)
                _write_response({"id": req_id, "ok": True, "result": result})
            except Exception:
                _write_response({"id": req_id, "ok": False,
                                 "error": traceback.format_exc()})
        elif op == "probe":
            try:
                result = _op_probe(cubit_mod, args)
                ok = not (isinstance(result, dict) and "error" in result)
                _write_response({"id": req_id, "ok": ok, "result": result})
            except Exception:
                _write_response({"id": req_id, "ok": False,
                                 "error": traceback.format_exc()})
        else:
            _write_response({"id": req_id, "ok": False,
                             "error": f"unknown op: {op!r}"})

    return 0


if __name__ == "__main__":
    sys.exit(main())
