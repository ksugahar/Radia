"""setup.py with a pre-build staleness check.

pyproject.toml drives metadata and package-data; this file exists only
to gate wheel/sdist creation on a freshness invariant:

    No bundled Cubit plugin binary (.ccm / .pyd) may be older
    than any source file under src/cubit_plugin/.

If the invariant is violated, the build aborts BEFORE setuptools bundles
the stale file into a wheel. Without this guard, ``pip install`` or
``pip wheel`` would happily package the latest .cpp changes' *non-built*
binaries, which is how 100号機 got a post-6a8d2e5 Python package with a
pre-6a8d2e5 .ccm on 2026-04-14.

Note (radia 4.80.0): the .ccl was removed (Qt5 GUI deleted; PySide6
toolbar at src/radia/panels/radia_export_menu.py replaces it).  Only
.ccm is now in the freshness gate.

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

    # Only .ccm is in this freshness gate — what Cubit loads directly
    # and is rebuilt by every Build.ps1 run. The .pyd (Cubit-less
    # HO-mesh module) is compiled only when the build machine has
    # pybind11 + full Netgen (non-compact), which is not the CI path;
    # its freshness is tracked separately by its own CMake target.
    #
    # Note (radia 4.80.0): the .ccl was removed from this gate (and
    # from the wheel package_data) because the Qt5 GUI .ccl was
    # deleted; PySide6 toolbar at radia/panels/radia_export_menu.py
    # replaces it.
    bundled = [
        pkg_dir / "cubit_mesh_export.ccm",
    ]
    present = [p for p in bundled if p.is_file()]
    missing = [p for p in bundled if not p.is_file()]
    if not present:
        # No binaries shipped at all — assume sdist / source-only build.
        # The wheel will still be installable (Python-only fallback paths)
        # but cubit-plugin-install will fail loud at deploy time.
        sys.stderr.write(
            "cubit-mesh-export: WARN — no bundled plugin binaries found "
            "in package source dir.\n")
        sys.stderr.write(
            "  Building sdist or a Python-only wheel. Wheel will lack "
            "cubit_mesh_export.ccm and cubit-plugin-install will refuse "
            "to deploy from it.\n")
        return
    if missing:
        sys.stderr.write(
            "\ncubit-mesh-export: FATAL — partial set of bundled binaries:\n")
        for p in present:
            sys.stderr.write(f"  [present] {p.name}\n")
        for p in missing:
            sys.stderr.write(f"  [MISSING] {p}\n")
        sys.stderr.write(
            "\n  All-or-nothing: either ship every binary or none. Rebuild "
            "and re-propagate (cubit_mesh_export.ccm is built by Build.ps1; the\n"
            "  .ccl was removed in radia 4.80.0) before retrying.\n\n")
        sys.exit(1)
    bundled = present  # only freshness-check the binaries we will ship

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
            "  OR: targeted rebuild — see the `release-qud` skill,\n"
            "  Phase 0. Then re-run `pip wheel` / `pip install`.\n"
            "\n  To override (not recommended):\n"
            "    set CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK=1\n\n")
        sys.exit(1)


_check_binary_freshness()
setup()
