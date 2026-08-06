"""
daemon.py — JSON-RPC server exposing a persistent Cubit session.

Runs under Cubit's bundled Python 3.10 (`C:/Program Files/Coreform Cubit
2025.12/bin/python3/python.exe`), NOT the system Python. This is because
the `cubit` Python binding is a compiled extension (`_cubit3.pyd`) tied
to the Cubit ABI + Python 3.10.

Role in the lab stack (2026-04-19 Cubit Viewer):
    VSCode Claude Code
        → mcp-server-cubit (Python 3.12)
            → THIS DAEMON (Python 3.10, subprocess)
                → cubit.cmd(...)
                    → Cubit GUI (FLTK, persistent window)

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
    "snapshot"  args=[path, w, h]  hardcopy current view to PNG
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

# Find Cubit binding directory (sibling of bin/python3/).
# Was hardcoded "Cubit 2025.3" until 2026-05-25; now auto-discovers the
# highest-version "Coreform Cubit *" under C:\Program Files. Honors
# CUBIT_BIN_DIR (legacy var) and CUBIT_INSTALL_DIR (preferred) overrides.
_HERE = os.path.dirname(os.path.abspath(__file__))


def _discover_cubit_bin():
    """Mirror of tools/find_cubit.py logic, inlined to avoid a
    sys.path bootstrap dependency in this daemon's import order."""
    if (override := os.environ.get("CUBIT_BIN_DIR")):
        return override
    if (inst := os.environ.get("CUBIT_INSTALL_DIR")):
        return os.path.join(inst, "bin")
    import glob as _g
    cands = _g.glob(r"C:\Program Files\Coreform Cubit *")
    if not cands:
        return r"C:\Program Files\Coreform Cubit 2025.12\bin"  # last-resort
    def _ver(p):
        name = os.path.basename(p).replace("Coreform Cubit", "", 1).strip()
        try:
            return tuple(int(x) for x in name.split("."))
        except ValueError:
            return (0,)
    return os.path.join(max(cands, key=_ver), "bin")


_CUBIT_BIN = _discover_cubit_bin()
if _CUBIT_BIN not in sys.path:
    sys.path.insert(0, _CUBIT_BIN)


def _write_response(obj):
    """Emit a single line of JSON to the preserved stdout fd and flush."""
    line = json.dumps(obj, default=str) + "\n"
    os.write(_REAL_STDOUT_FD, line.encode("utf-8"))


def _log(msg):
    """Debug log to stderr (mcp-server will inherit but not parse it)."""
    sys.stderr.write(f"[cubit_daemon] {msg}\n")
    sys.stderr.flush()


def _init_cubit(gui: bool):
    """Import cubit and initialize. Returns the imported module."""
    import cubit  # noqa: E402  (deferred import; path set above)

    # -nojournal: we don't need automatic .jou journaling.
    # -noecho:    don't re-echo commands to stdout (we'd corrupt JSON-RPC).
    # For GUI:    omit -batch and -nographics.
    # For headless: add -batch -nographics.
    argv = ["cubit", "-nojournal", "-noecho"]
    if not gui:
        argv.extend(["-batch", "-nographics"])
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
    dry-run into the live GUI.
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
    """Delegate to the SHARED probe implementation (probe_ops.op_probe)
    used by both the batch stdio runner (this file) and the GUI
    file-drop runner (bootstrap.py) -- single source, no query drift."""
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)
    import probe_ops
    return probe_ops.op_probe(cubit_mod, args)



def _op_snapshot(cubit_mod, args):
    """Delegate to the SHARED snapshot implementation (probe_ops)."""
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)
    import probe_ops
    return probe_ops.op_snapshot(cubit_mod, args)


def main():
    """Entry point: initialize Cubit, emit ready, enter RPC loop."""
    gui = os.environ.get("CUBIT_DAEMON_MODE", "gui").lower() != "batch"

    try:
        cubit_mod = _init_cubit(gui=gui)
        version = getattr(cubit_mod, "__version__", "unknown")
        _write_response({
            "ready": True,
            "mode": "gui" if gui else "batch",
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
        elif op == "snapshot":
            try:
                result = _op_snapshot(cubit_mod, args)
                ok = result.get("ok", False)
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
