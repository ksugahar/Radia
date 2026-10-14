"""Write one extended-validation mesh in an isolated Python process.

NGSolve 6.2.2606 on Windows writes geometry-free structured meshes completely
and then may fault while returning from ``Mesh.Save``. The parent generator
therefore validates the resulting .vol by reopening it and accepts the child
exit status only after that structural check succeeds.
"""

from __future__ import annotations

import sys
import os
import threading
import time
from pathlib import Path

import ngsolve as ng

from extended_catalog import build_mesh


def _exit_after_complete_write(path: Path) -> None:
    """Avoid the 6.2.2606 geometry-free Save return-path fault on Windows."""

    previous = -1
    stable = 0
    for _ in range(240):
        time.sleep(0.25)
        size = path.stat().st_size if path.is_file() else 0
        if size > 0 and size == previous:
            stable += 1
            if stable >= 8:
                os._exit(0)
        else:
            stable = 0
        previous = size


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: write_extended_mesh.py MESH_ID OUTPUT_VOL")
    with ng.TaskManager():
        mesh = build_mesh(sys.argv[1])
        path = Path(sys.argv[2])
        path.parent.mkdir(parents=True, exist_ok=True)
        netgen_mesh = mesh.ngmesh if isinstance(mesh, ng.Mesh) else mesh
        watchdog = threading.Thread(
            target=_exit_after_complete_write, args=(path,), daemon=True)
        watchdog.start()
        netgen_mesh.Save(str(path))


if __name__ == "__main__":
    main()
