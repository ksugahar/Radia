#!/usr/bin/env python3
"""Import one ACIS SAT in a disposable Cubit process and write volume evidence.

The input SAT is read-only.  The report is intentionally geometric only: it
does not classify materials or infer electrical winding data.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from audit_sat import audit_sat


_DLL_HANDLES: list[object] = []
_REPORT_PREFIX = "__RADIA_CUBIT_INVENTORY__"


def _find_cubit_bin(explicit: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    configured = os.environ.get("CUBIT_PATH")
    if configured:
        configured_path = Path(configured)
        candidates.append(
            configured_path.parent if configured_path.is_file()
            else configured_path)
    if os.name == "nt":
        candidates.extend(sorted(
            Path(r"C:\Program Files").glob("Coreform Cubit */bin"),
            reverse=True))
    for candidate in candidates:
        if (candidate / "coreform_cubit.exe").is_file():
            return candidate.resolve()
    raise RuntimeError("Coreform Cubit 2025.12+ was not found; set CUBIT_PATH")


def _load_cubit(cubit_bin: Path):
    sys.path.insert(0, str(cubit_bin))
    os.environ["PATH"] = str(cubit_bin) + os.pathsep + os.environ.get(
        "PATH", "")
    if hasattr(os, "add_dll_directory"):
        for directory in (cubit_bin, cubit_bin / "python3"):
            if directory.is_dir():
                _DLL_HANDLES.append(os.add_dll_directory(str(directory)))
    import cubit

    cubit.init([
        "cubit", "-nojournal", "-batch", "-nographics",
        "-commandplugindir", str(cubit_bin / "plugins"),
    ])
    return cubit


def _as_float_list(value: Any, count: int, name: str) -> list[float]:
    result = [float(item) for item in value]
    if len(result) != count or not np.all(np.isfinite(result)):
        raise RuntimeError(f"Cubit returned an invalid {name}: {value!r}")
    return result


def _inventory(cubit, source_sat: Path) -> dict[str, Any]:
    cubit.cmd("reset")
    cubit.cmd(f'import acis "{source_sat.as_posix()}"')
    volume_ids = tuple(int(value) for value in cubit.get_entities("volume"))
    if not volume_ids:
        raise RuntimeError("Cubit imported zero volumes from SAT")
    volumes = []
    for volume_id in volume_ids:
        name = str(cubit.get_entity_name("volume", volume_id) or "")
        surface_ids = tuple(
            int(value) for value in cubit.get_relatives(
                "volume", volume_id, "surface"))
        volumes.append({
            "id": volume_id,
            "name": name,
            # Cubit 2025.12 exposes ``arrayd10`` here, not a documented
            # six-coordinate min/max tuple.  Preserve the native sequence
            # losslessly rather than assigning a made-up axis convention.
            "cubit_bounding_box_arrayd10": _as_float_list(
                cubit.get_bounding_box("volume", volume_id), 10,
                "volume bounding box"),
            "centre": _as_float_list(
                cubit.get_center_point("volume", volume_id), 3,
                "volume centre"),
            "surface_count": len(surface_ids),
            "surfaces": [
                {
                    "id": surface_id,
                    "centre": _as_float_list(
                        cubit.get_center_point("surface", surface_id), 3,
                        "surface centre"),
                    "vertices": [
                        _as_float_list(
                            cubit.get_center_point("vertex", vertex_id), 3,
                            "surface vertex centre")
                        for vertex_id in tuple(int(value) for value in
                                                cubit.get_relatives(
                                                    "surface", surface_id,
                                                    "vertex"))
                    ],
                }
                for surface_id in surface_ids
            ],
        })
    sat_bodies = audit_sat(source_sat)["bodies"]
    import_order_candidate = len(sat_bodies) == len(volumes)
    if import_order_candidate:
        for volume, body in zip(volumes, sat_bodies):
            volume["sat_import_order_candidate"] = {
                "sat_body_id": body["sat_body_id"],
                "sat_name": body["name"],
                "verification": "required-cubit-geometry-review",
            }
    return {
        "schema": "radia.cubit-sat-volume-inventory/v1",
        "source_sat": str(source_sat.resolve()),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cubit_version": str(cubit.get_version()),
        "volume_count": len(volumes),
        "sat_import_order_candidate_available": import_order_candidate,
        "volumes": volumes,
        "warning": (
            "This import inventory is geometric evidence only. Assign iron, "
            "coil, air-gap, and discard roles in the DesignSpec after review. "
            "SAT import order, when shown, is an unverified candidate only."
        ),
    }


def _child(args: argparse.Namespace) -> int:
    cubit = _load_cubit(_find_cubit_bin(args.cubit_bin))
    report = _inventory(cubit, args.sat.resolve())
    print(_REPORT_PREFIX + json.dumps(report, ensure_ascii=False))
    return 0


def _parent(args: argparse.Namespace) -> int:
    command = [
        sys.executable, str(Path(__file__).resolve()),
        "--cubit-child", str(args.sat.resolve()),
    ]
    if args.cubit_bin is not None:
        command.extend(("--cubit-bin", str(args.cubit_bin.resolve())))
    completed = subprocess.run(command, capture_output=True, text=True,
                               errors="replace", check=False)
    report = None
    for line in completed.stdout.splitlines():
        if line.startswith(_REPORT_PREFIX):
            report = json.loads(line.removeprefix(_REPORT_PREFIX))
    if completed.returncode or report is None:
        tail = "\n".join((completed.stdout + completed.stderr).splitlines()[-60:])
        raise RuntimeError(
            "Cubit SAT inventory failed "
            f"(exit={completed.returncode}):\n{tail}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8")
    print(
        f"Cubit inventory: version={report['cubit_version']} "
        f"volumes={report['volume_count']}")
    return 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sat", type=Path, help="input ACIS SAT file")
    parser.add_argument("--output", type=Path, required=False,
                        help="JSON inventory output path")
    parser.add_argument("--cubit-bin", type=Path)
    parser.add_argument("--cubit-child", action="store_true",
                        help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if not args.sat.is_file():
        parser.error(f"SAT file does not exist: {args.sat}")
    if not args.cubit_child and args.output is None:
        parser.error("--output is required outside the Cubit child process")
    if args.output is not None and args.output.resolve() == args.sat.resolve():
        parser.error("output must not overwrite the SAT file")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    return _child(args) if args.cubit_child else _parent(args)


if __name__ == "__main__":
    raise SystemExit(main())
