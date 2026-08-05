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


_ENTITY_DUMP_CAP = 300


def _round3(seq, ndigits=6):
    return [round(float(v), ndigits) for v in seq]


def _bbox_fields(bb):
    """Cubit's GeomEntity.bounding_box() returns a 6-tuple
    (min_x, min_y, min_z, max_x, max_y, max_z) -- VERIFIED on Coreform
    Cubit 2025.12 (brick x 1 -> [-0.5]*3 + [0.5]*3).  Extent is derived;
    a flat face has zero extent in its normal direction."""
    lo, hi = bb[0:3], bb[3:6]
    return {
        "bbox_min": _round3(lo),
        "bbox_max": _round3(hi),
        "extent": _round3([h - l for l, h in zip(lo, hi)]),
    }


def _op_probe_entities(cubit_mod):
    """Per-entity geometric dump implementing the lab's Probe-Don't-Guess
    policy: every volume/surface with centroid, bbox, extent, and measure,
    so entity-classification predicates are derived from printed values
    rather than a-priori reasoning about the .jou source."""
    vol_ids = [int(v) for v in cubit_mod.parse_cubit_list("volume", "all")]
    surf_ids = [int(s) for s in cubit_mod.parse_cubit_list("surface", "all")]
    result = {
        "volume_count": len(vol_ids),
        "surface_count": len(surf_ids),
        "volumes": [],
        "surfaces": [],
    }
    for vid in vol_ids[:_ENTITY_DUMP_CAP]:
        v = cubit_mod.volume(vid)
        entry = {
            "id": vid,
            "centroid": _round3(v.centroid()),
            "volume": float(v.volume()),
        }
        entry.update(_bbox_fields(v.bounding_box()))
        result["volumes"].append(entry)
    for sid in surf_ids[:_ENTITY_DUMP_CAP]:
        s = cubit_mod.surface(sid)
        cx, cy, cz = s.center_point()   # Surface API: NOT .centroid()
        entry = {
            "id": sid,
            "center": _round3((cx, cy, cz)),
            "area": float(s.area()),
        }
        entry.update(_bbox_fields(s.bounding_box()))
        result["surfaces"].append(entry)
    if len(vol_ids) > _ENTITY_DUMP_CAP or len(surf_ids) > _ENTITY_DUMP_CAP:
        result["truncated"] = (
            f"listing capped at {_ENTITY_DUMP_CAP} entities per kind; "
            "use parse_cubit_list + per-entity calls via cubit_exec for the rest"
        )
    return result


def _op_probe_labels(cubit_mod):
    """Blocks/sidesets with names + convention audit (phantom/mixed-block
    detection before `export netgen`).

    VERIFIED on Coreform Cubit 2025.12: `block N add <other-kind>` on an
    already-typed block returns success but adds NOTHING (silent no-op),
    and `block N add tri ...` on a hex mesh creates NO block at all
    (phantom).  The primary value of this probe is therefore verifying
    ACTUAL block membership after every block/sideset command -- an
    intended-but-absent element set shows up here as missing counts."""
    blocks = []
    for bid in cubit_mod.get_block_id_list():
        bid = int(bid)
        vol_elems = (len(cubit_mod.get_block_hexes(bid))
                     + len(cubit_mod.get_block_tets(bid))
                     + len(cubit_mod.get_block_wedges(bid))
                     + len(cubit_mod.get_block_pyramids(bid)))
        surf_elems = (len(cubit_mod.get_block_tris(bid))
                      + len(cubit_mod.get_block_faces(bid)))
        blocks.append({
            "id": bid,
            "name": str(cubit_mod.get_exodus_entity_name("block", bid)),
            "volumes": [int(v) for v in cubit_mod.get_block_volumes(bid)],
            "surfaces": [int(s) for s in cubit_mod.get_block_surfaces(bid)],
            "volume_elems": int(vol_elems),
            "surface_elems": int(surf_elems),
        })
    sidesets = []
    for sid in cubit_mod.get_sideset_id_list():
        sid = int(sid)
        sidesets.append({
            "id": sid,
            "name": str(cubit_mod.get_exodus_entity_name("sideset", sid)),
            "surfaces": [int(s) for s in cubit_mod.get_sideset_surfaces(sid)],
        })
    # daemon.py runs standalone under Cubit's Python (no package context);
    # import the shared cubit-free audit module by path.
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)
    from label_audit import audit_label_records
    audit = audit_label_records(blocks, sidesets)
    return {"blocks": blocks, "sidesets": sidesets, "audit": audit}


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
        if q in ("entities", "geometry"):
            return _op_probe_entities(cubit_mod)
        if q in ("labels", "blocks", "block_sideset"):
            return _op_probe_labels(cubit_mod)
        if q in ("quality_summary", "quality"):
            try:
                hex_ids = list(cubit_mod.parse_cubit_list("hex", "all"))
                tet_ids = list(cubit_mod.parse_cubit_list("tet", "all"))
            except Exception:
                hex_ids, tet_ids = [], []
            vals = []
            try:
                if hex_ids:
                    vals.extend(float(v) for v in cubit_mod.get_quality_values(
                        "hex", [int(x) for x in hex_ids], "scaled jacobian"))
            except Exception:
                pass
            try:
                if tet_ids:
                    vals.extend(float(v) for v in cubit_mod.get_quality_values(
                        "tet", [int(x) for x in tet_ids], "scaled jacobian"))
            except Exception:
                pass
            if not vals:
                return {"hex_count": len(hex_ids), "tet_count": len(tet_ids),
                        "min": None, "max": None, "mean": None,
                        "below_0.2": 0}
            below = sum(1 for v in vals if v < 0.2)
            return {
                "hex_count": int(len(hex_ids)),
                "tet_count": int(len(tet_ids)),
                "min": round(min(vals), 6),
                "max": round(max(vals), 6),
                "mean": round(sum(vals) / len(vals), 6),
                "below_0.2": int(below),
                "below_0.2_pct": round(100 * below / len(vals), 2),
            }
    except Exception:
        return {"error": traceback.format_exc()}
    return {"error": f"Unknown probe query: {query!r}. "
                     "Try: volume_count / surface_count / vertex_count / "
                     "curve_count / node_count / hex_count / tet_count / "
                     "summary / quality / entities / labels"}


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
