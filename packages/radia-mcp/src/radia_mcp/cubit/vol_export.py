"""Curved Netgen ``.vol`` export worker: Cubit plus the native ``export netgen`` plugin.

Runs in an isolated external CPython (``python -I -B -m
radia_mcp.cubit.vol_export``), never inside the MCP server or Cubit's embedded
Python.  Order of operations:

1. ``import ngsolve`` first, so the fork nglib the plugin was built against is
   the one already loaded in the process;
2. add the Netgen package and Cubit ``bin`` directories to the DLL search path;
3. point ``CUBIT_PLUGIN_DIR`` at the directory holding ``cubit_mesh_export.dll``;
4. ``import cubit`` and ``cubit.init([... '-noinitfile'])``;
5. play the journal line by line;
6. optionally add the Kelvin open-boundary sphere (:mod:`radia_mcp.cubit.kelvin`);
7. ``export netgen "<out>" order <p> [periodic args] overwrite``.

The plugin writes ``<out>`` and the companion ``<out>.json``.  The worker
writes a status JSON (``--result``) with the companion record and the sha256
of the ``.vol``.  Every failure is reported there and the exit code is 1.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import traceback
from pathlib import Path

PLUGIN_DLL = "cubit_mesh_export.dll"
PERIODIC_FIELDS = ("periodic_name", "periodic_master", "periodic_slave", "periodic_unit",
                   "periodic_patches", "periodic_tx", "periodic_ty", "periodic_tz", "periodic_match")
_PERIODIC_TEXT = {"periodic_name", "periodic_master", "periodic_slave", "periodic_unit", "periodic_patches"}


def periodic_arguments(periodic: dict | None) -> str:
    """Render the named periodic request as ``export netgen`` arguments (all nine or none)."""
    if not periodic:
        return ""
    unknown = set(periodic) - set(PERIODIC_FIELDS)
    if unknown:
        raise ValueError(f"unknown periodic field(s): {sorted(unknown)}")
    missing = [name for name in PERIODIC_FIELDS if periodic.get(name) in (None, "")]
    if missing:
        raise ValueError(f"named periodic export needs all nine fields; missing {missing}")
    parts = []
    for name in PERIODIC_FIELDS:
        value = str(periodic[name])
        if any(ch in value for ch in '"\r\n'):
            raise ValueError(f"{name} must not contain quotes or newlines")
        parts.append(f'{name} "{value}"' if name in _PERIODIC_TEXT else f"{name} {value}")
    return " ".join(parts)


def export_command(out: Path, order: int, periodic: dict | None) -> str:
    extra = periodic_arguments(periodic)
    return f'export netgen "{out.as_posix()}" order {int(order)} ' + (extra + " " if extra else "") + "overwrite"


def run(args: argparse.Namespace) -> dict:
    out = Path(args.out).resolve()
    plugin_dir = Path(args.plugin_dir).resolve()
    if not (plugin_dir / PLUGIN_DLL).is_file():
        raise FileNotFoundError(f"{PLUGIN_DLL} not found in plugin_dir {plugin_dir}")
    periodic = json.loads(args.periodic) if args.periodic else None
    kelvin_options = json.loads(args.kelvin) if args.kelvin else {}
    from radia_mcp.cubit import kelvin as _kelvin
    _kelvin.validate_options(kelvin_options)
    command = export_command(out, args.order, periodic)

    import ngsolve  # noqa: F401  -- load the fork nglib before the plugin binds to it
    import netgen
    os.add_dll_directory(os.path.dirname(netgen.__file__))
    cubit_bin = str(Path(args.cubit_bin).resolve())
    sys.path.append(cubit_bin)
    os.add_dll_directory(cubit_bin)
    os.environ["CUBIT_PLUGIN_DIR"] = str(plugin_dir)
    import cubit
    cubit.init(["cubit", "-nographics", "-nojournal", "-batch", "-noinitfile"])
    started = time.perf_counter()
    journal = Path(args.journal).resolve()
    # Cubit plays the journal itself, so continuation lines, comments and
    # APREPRO blocks keep their journal meaning.
    if not cubit.cmd(f'playback "{journal.as_posix()}"'):
        raise RuntimeError(f"Cubit journal playback failed: {journal}")
    kelvin_info = _kelvin.apply_options(kelvin_options) if kelvin_options else None
    cubit.cmd(command)
    if not out.is_file():
        raise RuntimeError(f"export netgen wrote no file: {command}")
    companion_path = Path(str(out) + ".json")
    companion = json.loads(companion_path.read_text(encoding="utf-8")) if companion_path.is_file() else None
    if companion is None:
        raise RuntimeError(f"export netgen wrote no companion JSON: {companion_path}")
    return {"status": "ok", "vol": str(out), "sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
            "companion_json": str(companion_path), "companion": companion,
            "companion_sha256": hashlib.sha256(companion_path.read_bytes()).hexdigest(),
            "export_command": command, "plugin": str(plugin_dir / PLUGIN_DLL),
            "kelvin": kelvin_info, "netgen": getattr(netgen, "__version__", None),
            "elapsed_s": round(time.perf_counter() - started, 3)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Curved .vol export through the native Cubit plugin")
    parser.add_argument("--cubit-bin", required=True)
    parser.add_argument("--plugin-dir", required=True)
    parser.add_argument("--journal", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--order", type=int, required=True)
    parser.add_argument("--result", required=True)
    parser.add_argument("--token", required=True, help="run token echoed into the status file")
    parser.add_argument("--periodic", default="", help="JSON object with the nine periodic_* fields")
    parser.add_argument("--kelvin", default="", help="JSON object with the Kelvin options")
    args = parser.parse_args(argv)
    try:
        status = run(args)
        code = 0
    except Exception as exc:  # reported through the status file, never swallowed
        status = {"status": "error", "error_type": type(exc).__name__, "error": str(exc),
                  "traceback": traceback.format_exc()[-4000:]}
        code = 1
    status["token"] = args.token
    Path(args.result).write_text(json.dumps(status, indent=1, default=str), encoding="utf-8")
    return code


if __name__ == "__main__":
    sys.exit(main())
