"""
cubit_daemon.py — JSON-RPC server exposing a persistent Cubit session.

Runs under Cubit's bundled Python 3.10 (`C:/Program Files/Coreform Cubit
2025.3/bin/python3/python.exe`), NOT the system Python. This is because
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
_HERE = os.path.dirname(os.path.abspath(__file__))
_CUBIT_BIN = os.environ.get(
    "CUBIT_BIN_DIR",
    r"C:\Program Files\Coreform Cubit 2025.3\bin",
)
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


def _op_cmd(cubit_mod, args):
    """Execute a list of cubit commands. Returns per-command success flags.

    Cubit's `cmd` returns an int (1 on success, 0 on failure) per call.
    We run them sequentially and stop at the first failure so downstream
    state is recoverable.
    """
    results = []
    for line in args:
        try:
            rc = cubit_mod.cmd(str(line))
        except Exception:
            results.append({"line": line, "ok": False,
                            "error": traceback.format_exc()})
            break
        ok = bool(rc)
        results.append({"line": line, "ok": ok, "rc": int(rc) if not isinstance(rc, bool) else int(rc)})
        if not ok:
            break
    return results


def _op_probe(cubit_mod, args):
    """Return lab-common Cubit state queries by short name."""
    query = args[0] if args else ""
    q = str(query).strip().lower()
    try:
        if q in ("volume_count", "volumes"):
            return int(cubit_mod.get_volume_count())
        if q in ("surface_count", "surfaces"):
            return int(cubit_mod.get_surface_count())
        if q in ("vertex_count", "vertices"):
            return int(cubit_mod.get_vertex_count())
        if q in ("curve_count", "curves"):
            return int(cubit_mod.get_curve_count())
        if q in ("node_count", "nodes"):
            return int(cubit_mod.get_node_count())
        if q in ("hex_count", "hexes", "hex"):
            return int(cubit_mod.get_hex_count())
        if q in ("tet_count", "tets", "tet"):
            return int(cubit_mod.get_tet_count())
        if q in ("summary",):
            return {
                "volumes": int(cubit_mod.get_volume_count()),
                "surfaces": int(cubit_mod.get_surface_count()),
                "curves": int(cubit_mod.get_curve_count()),
                "vertices": int(cubit_mod.get_vertex_count()),
                "nodes": int(cubit_mod.get_node_count()),
                "hexes": int(cubit_mod.get_hex_count()),
                "tets": int(cubit_mod.get_tet_count()),
            }
    except Exception:
        return {"error": traceback.format_exc()}
    return {"error": f"Unknown probe query: {query!r}. "
                     "Try: volume_count / surface_count / vertex_count / "
                     "curve_count / node_count / hex_count / tet_count / summary"}


def _op_snapshot(cubit_mod, args):
    """Hardcopy the current Cubit view to a PNG file."""
    if len(args) < 1:
        return {"error": "snapshot requires at least 1 arg (path)"}
    out_path = str(args[0])
    w = int(args[1]) if len(args) >= 2 else 800
    h = int(args[2]) if len(args) >= 3 else 600
    # Cubit hardcopy command:
    #   hardcopy 'file.png' png window [wxh]
    # GUI mode is recommended for this; batch mode may produce empty
    # images depending on Cubit version.
    cmd_str = f'hardcopy "{out_path}" png window {w} {h}'
    rc = cubit_mod.cmd(cmd_str)
    return {"path": out_path, "ok": bool(rc), "width": w, "height": h}


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

    for line in sys.stdin:
        line = line.strip()
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
