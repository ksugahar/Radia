"""setup.py with a content-addressed native provenance check.

pyproject.toml drives metadata and package-data; this file exists only
to gate wheel/sdist creation on a provenance invariant:

    The native source-tree digest and both bundled payload digests must match
    native_payloads.json. Filesystem timestamps are never evidence.

If the invariant is violated, the build aborts BEFORE setuptools bundles
the stale file into a wheel. Without this guard, ``pip install`` or
``pip wheel`` would happily package the latest .cpp changes' *non-built*
binaries, which is how 100号機 got a post-6a8d2e5 Python package with a
pre-6a8d2e5 .ccm on 2026-04-14.

Note (radia 4.80.0): the .ccl was removed (Qt5 GUI deleted; PySide6
dialogs and the Claro-owned menu replace it). The .ccm and .pyd are both
mandatory wheel payloads and are both freshness-gated.

Override with ``CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK=1`` only as a
last resort (e.g. emergency release when the build box is offline).
"""

from __future__ import annotations

import os
import sys
import importlib.util
from pathlib import Path

from setuptools import Distribution, setup


class BinaryDistribution(Distribution):
    """Mark the prebuilt .pyd/.ccm payload as a native wheel."""

    def has_ext_modules(self):
        return True


def _load_provenance_module(package_dir: Path):
    path = package_dir / "_native_provenance.py"
    spec = importlib.util.spec_from_file_location("cme_native_provenance", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load native provenance helper: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _check_binary_provenance():
    # Env-var escape hatch. Use sparingly.
    if os.environ.get("CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK") == "1":
        print("cubit-mesh-export: native provenance check SKIPPED via env var.",
              file=sys.stderr)
        return

    here = Path(__file__).resolve().parent
    pkg_dir = here / "src" / "cubit_mesh_export"
    # Walk up the monorepo to find src/cubit_plugin. Two levels up from
    # packages/cubit-mesh-export/.
    repo_root = here.parent.parent
    cpp_dir = repo_root / "src" / "cubit_plugin"

    if not cpp_dir.is_dir():
        # An sdist omits C++ sources, not the mandatory payload contract.
        # Skip only source-tree comparison; still verify schema and binaries.
        repo_root = None

    provenance = _load_provenance_module(pkg_dir)
    errors = provenance.verify_manifest(repo_root, pkg_dir)
    if errors:
        sys.stderr.write(
            "\ncubit-mesh-export: FATAL — native provenance mismatch.\n")
        for error in errors:
            sys.stderr.write(f"  - {error}\n")
        sys.stderr.write(
            "\n  Rebuild BEFORE packaging:\n"
            "    pwsh -File src/cubit_plugin/cubit_build.ps1 -Rebuild\n"
            "  This rebuilds and SHA-verifies both mandatory payloads.\n"
            "  Then re-run `pip wheel` / `pip install`.\n"
            "\n  To override (not recommended):\n"
            "    set CUBIT_MESH_EXPORT_SKIP_FRESHNESS_CHECK=1\n\n")
        sys.exit(1)


_check_binary_provenance()
setup(distclass=BinaryDistribution)
