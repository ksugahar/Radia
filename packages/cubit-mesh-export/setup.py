"""setup.py with a pre-build staleness check.

pyproject.toml drives metadata and package-data; this file exists only
to gate wheel/sdist creation on a freshness invariant:

    No bundled Cubit plugin binary (.ccm / .ccl / .pyd) may be older
    than any source file under src/cubit_plugin/.

If the invariant is violated, the build aborts BEFORE setuptools bundles
the stale file into a wheel. Without this guard, ``pip install`` or
``pip wheel`` would happily package the latest .cpp changes' *non-built*
binaries, which is how 100号機 got a post-6a8d2e5 Python package with a
pre-6a8d2e5 .ccl on 2026-04-14.

Override with ``CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK=1`` only as a
last resort (e.g. emergency release when the build box is offline).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

from setuptools import setup


def _check_binary_freshness():
    # Env-var escape hatch. Use sparingly.
    if os.environ.get("CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK") == "1":
        print("cubit-mesh-export: freshness check SKIPPED via env var.",
              file=sys.stderr)
        return

    here = Path(__file__).resolve().parent
    pkg_dir = here / "src" / "cubit_mesh_export"
    # Walk up the monorepo to find src/cubit_plugin. Two levels up from
    # packages/cubit-mesh-export/.
    repo_root = here.parent.parent
    cpp_dir = repo_root / "src" / "cubit_plugin"

    if not cpp_dir.is_dir():
        # Building from an sdist — the C++ source isn't shipped, so the
        # freshness check doesn't apply. The sdist was already built
        # with the guard upstream.
        return

    # Only .ccm and .ccl are in this freshness gate — these are what Cubit
    # loads directly and are rebuilt by every Build.ps1 run. The .pyd
    # (Cubit-less HO-mesh module) is compiled only when the build machine
    # has pybind11 + full Netgen (non-compact), which is not the CI path;
    # its freshness is tracked separately by its own CMake target.
    bundled = [
        pkg_dir / "radia_cubit.ccm",
        pkg_dir / "radia_cubit.ccl",
    ]
    missing = [p for p in bundled if not p.is_file()]
    if missing:
        sys.stderr.write(
            "\ncubit-mesh-export: FATAL — expected bundled binary missing:\n")
        for p in missing:
            sys.stderr.write(f"  - {p}\n")
        sys.stderr.write(
            "\n  Rebuild the Cubit plugin and propagate the artifacts:\n"
            "    powershell Build.ps1 -Rebuild   # from repo root\n\n")
        sys.exit(1)

    # Find the latest mtime among .cpp/.hpp/.h/.cmake under src/cubit_plugin/.
    latest_src_mtime = 0.0
    latest_src_file = None
    SRC_EXT = {".cpp", ".cc", ".cxx", ".c", ".h", ".hpp", ".hh", ".hxx",
               ".cmake", ".txt"}  # .txt for CMakeLists.txt
    # Skip obviously regenerated subdirs.
    SKIP_DIRS = {"build-pyd", "build-ccm", "build", "compact_netgen"}
    for root, dirs, files in os.walk(cpp_dir):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for f in files:
            if Path(f).suffix.lower() in SRC_EXT:
                mt = (Path(root) / f).stat().st_mtime
                if mt > latest_src_mtime:
                    latest_src_mtime = mt
                    latest_src_file = Path(root) / f

    stale = []
    for b in bundled:
        mt = b.stat().st_mtime
        if mt + 1 < latest_src_mtime:  # 1 s tolerance for fs granularity
            stale.append((b, mt))

    if stale:
        sys.stderr.write(
            "\ncubit-mesh-export: FATAL — bundled Cubit plugin binaries "
            "are older than src/cubit_plugin/ source files.\n")
        sys.stderr.write(
            f"\n  Newest source: {latest_src_file}\n"
            f"    mtime: {latest_src_mtime}\n")
        sys.stderr.write("  Stale bundled binaries:\n")
        for b, mt in stale:
            delta_h = (latest_src_mtime - mt) / 3600.0
            sys.stderr.write(f"    - {b} (mtime {mt}, {delta_h:.1f} h older)\n")
        sys.stderr.write(
            "\n  Rebuild BEFORE packaging:\n"
            "    powershell.exe -ExecutionPolicy Bypass "
            "-File Build.ps1 -Rebuild    # full rebuild (~10 min)\n"
            "  OR: targeted rebuild — see the `release-triple` skill,\n"
            "  Phase 0. Then re-run `pip wheel` / `pip install`.\n"
            "\n  To override (not recommended):\n"
            "    set CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK=1\n\n")
        sys.exit(1)


_check_binary_freshness()
setup()
