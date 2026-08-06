"""Build mapped HEX meshes for the thin magnetic-conductor disk lane.

Run this script with Cubit's bundled Python.  The generated ``.vol`` is a
validation input and is intentionally not tracked; this generator is the
reproducible source of truth.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


_DLL_HANDLES = []


def _find_cubit_bin(explicit: Path | None) -> Path:
    candidates = []
    if explicit is not None:
        candidates.append(explicit)
    configured = os.environ.get("CUBIT_PATH")
    if configured:
        path = Path(configured)
        candidates.append(path.parent if path.is_file() else path)
    if os.name == "nt":
        candidates.extend(
            sorted(
                Path(r"C:\Program Files").glob("Coreform Cubit */bin"),
                reverse=True,
            )
        )
    for candidate in candidates:
        if (candidate / "coreform_cubit.exe").is_file():
            return candidate.resolve()
    raise RuntimeError("Coreform Cubit 2025.12+ was not found; set CUBIT_PATH")


def _load_cubit(cubit_bin: Path):
    sys.path.insert(0, str(cubit_bin))
    os.environ["PATH"] = str(cubit_bin) + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        for directory in (cubit_bin, cubit_bin / "python3"):
            if directory.is_dir():
                _DLL_HANDLES.append(os.add_dll_directory(str(directory)))
    import cubit

    cubit.init(
        [
            "cubit",
            "-nojournal",
            "-batch",
            "-nographics",
            "-commandplugindir",
            str(cubit_bin / "plugins"),
        ]
    )
    return cubit


def build_mesh(
    cubit,
    *,
    size_mm: float,
    output: Path,
    positive_octant: bool,
) -> dict[str, object]:
    if size_mm <= 0.0:
        raise ValueError("size_mm must be positive")
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    commands = [
        "reset",
        "cylinder height 0.5 radius 10",
        "webcut volume all with plane xplane offset 0 imprint merge",
        "webcut volume all with plane yplane offset 0 imprint merge",
    ]
    if positive_octant:
        commands.append("webcut volume all with plane zplane offset 0 imprint merge")
    commands.extend(("imprint all", "merge all"))
    for command in commands:
        cubit.cmd(command)

    if positive_octant:
        for volume_id in cubit.parse_cubit_list("volume", "all"):
            if any(
                coordinate < -1.0e-10
                for coordinate in cubit.get_center_point("volume", volume_id)
            ):
                cubit.cmd(f"delete volume {volume_id}")
        cubit.cmd("compress")

    for command in (
        f"volume all size {size_mm:g}",
        "mesh volume all",
        "block 1 add hex all",
        'block 1 name "conductor"',
        "volume all scale 0.001",
        f'export netgen "{output.as_posix()}" order 1 overwrite',
    ):
        cubit.cmd(command)

    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"Cubit did not create {output}")
    return {
        "output": str(output),
        "size_mm": float(size_mm),
        "positive_octant": bool(positive_octant),
        "hexes": int(cubit.get_hex_count()),
        "nodes": int(cubit.get_node_count()),
        "bytes": int(output.stat().st_size),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--size-mm", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--positive-octant", action="store_true")
    parser.add_argument("--cubit-bin", type=Path)
    args = parser.parse_args()
    cubit_bin = _find_cubit_bin(args.cubit_bin)
    cubit = _load_cubit(cubit_bin)
    print(
        json.dumps(
            build_mesh(
            cubit,
            size_mm=args.size_mm,
            output=args.output,
            positive_octant=args.positive_octant,
            ),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
