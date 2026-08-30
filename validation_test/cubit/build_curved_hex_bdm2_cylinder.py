"""Build the canonical Cubit 2025.12 curved-HEX BDM2 validation mesh.

The generated order-2 Netgen ``.vol`` is temporary input for
``validate_curved_hex_bdm2_cubit.py``.  Cubit/ACIS owns the cylinder geometry
and ``cubit-mesh-export`` owns the high-order node placement.  No mesh artifact
is tracked in Git.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from radia.install_panels import find_cubit_bin

DEFAULT_OUTPUT = Path(r"C:\temp\radia_hex_bdm2\curved_cylinder_hex_q2.vol")


def build(output: Path) -> Path:
    cubit_bin = find_cubit_bin()
    if not cubit_bin:
        raise RuntimeError("Coreform Cubit 2025.12 is required")
    plugin_dir = Path(cubit_bin).parent / "bin" / "plugins"
    if not plugin_dir.is_dir():
        raise RuntimeError(f"Cubit plugin directory does not exist: {plugin_dir}")
    sys.path.append(str(cubit_bin))
    os.environ["CUBIT_PLUGIN_DIR"] = str(plugin_dir)

    import cubit

    cubit.init(["cubit", "-nojournal", "-batch", "-commandplugindir", str(plugin_dir)])
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    for command in (
        "reset",
        "cylinder height 0.04 radius 0.02",
        "volume 1 size 0.02",
        "mesh volume 1",
        "block 1 add hex all",
        'block 1 name "iron"',
    ):
        cubit.cmd(command)
    if int(cubit.get_hex_count()) != 4:
        raise RuntimeError(
            "canonical curved-cylinder mesh must contain exactly 4 HEX cells "
            f"(got {cubit.get_hex_count()})"
        )
    cubit.cmd(f'export netgen "{output}" order 2 overwrite')
    if not output.is_file():
        raise RuntimeError(f"Cubit did not create {output}")
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(build(args.output), flush=True)


if __name__ == "__main__":
    main()
