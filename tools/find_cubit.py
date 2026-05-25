"""Find the latest Coreform Cubit install. Python helper.

Companion to ``tools/find_cubit.ps1`` for PowerShell. Both follow the
same discovery rules:

  1. ``$CUBIT_INSTALL_DIR`` env var if set and exists.
  2. Latest ``C:\\Program Files\\Coreform Cubit *`` directory
     (version-aware sort: 2025.12 > 2025.3).
  3. Returns None + logs a warning if nothing found.

Replaces the prior hardcoded path pattern in
``radia_mcp.cubit.daemon`` and ``radia.install_panels``:

    # OLD (hardcoded -- breaks on Cubit upgrade):
    candidates = [r"C:\\Program Files\\Coreform Cubit 2025.3\\bin", ...]

    # NEW (auto-discovery):
    from tools.find_cubit import find_cubit_bin
    cubit_bin = find_cubit_bin()

Or, since ``tools/`` is not a package, copy this small helper into
the consumer file or vendor as needed.

Standalone CLI:

    python tools/find_cubit.py             # human-readable
    python tools/find_cubit.py --json      # machine-readable JSON

Exit code: 0 on success, 1 if no install found.
"""
from __future__ import annotations

import os
import sys
from functools import lru_cache
from pathlib import Path
from typing import Optional


def _parse_version(name: str) -> tuple:
    """Parse 'Coreform Cubit 2025.12' -> (2025, 12) for version sort."""
    stem = name.replace("Coreform Cubit", "", 1).strip()
    parts = stem.split(".")
    out = []
    for p in parts:
        try:
            out.append(int(p))
        except ValueError:
            out.append(0)
            break
    return tuple(out) if out else (0,)


@lru_cache(maxsize=1)
def find_cubit_install() -> Optional[Path]:
    """Return Path to latest Coreform Cubit install, or None.

    Cached: first call scans Program Files, subsequent calls return
    the cached Path. To force rediscovery (rare; e.g. after a Cubit
    install/upgrade within the same Python process), call
    ``find_cubit_install.cache_clear()``.
    """
    # 1. Explicit env override
    env = os.environ.get("CUBIT_INSTALL_DIR")
    if env:
        p = Path(env)
        if p.is_dir():
            return p

    # 2. Scan C:\Program Files for the highest-version Coreform Cubit
    base = Path(r"C:\Program Files")
    if not base.is_dir():
        return None
    cands = [p for p in base.iterdir()
             if p.is_dir() and p.name.startswith("Coreform Cubit")]
    if not cands:
        return None
    return max(cands, key=lambda p: _parse_version(p.name))


def find_cubit_bin() -> Optional[Path]:
    """Return Path to the bin/ subdir, or None if no install."""
    inst = find_cubit_install()
    return inst / "bin" if inst else None


def find_cubit_cmake() -> Optional[Path]:
    """Return Path to the cmake/ subdir (for CubitConfig.cmake), or None."""
    inst = find_cubit_install()
    return inst / "cmake" if inst else None


def find_cubit_exe() -> Optional[Path]:
    """Return Path to coreform_cubit.exe, or None."""
    bin_dir = find_cubit_bin()
    if not bin_dir:
        return None
    exe = bin_dir / "coreform_cubit.exe"
    return exe if exe.exists() else None


def get_cubit_version() -> Optional[str]:
    """Return the version string '2025.12' (parsed from dir name), or None."""
    inst = find_cubit_install()
    if not inst:
        return None
    return inst.name.replace("Coreform Cubit", "", 1).strip()


def _main(argv):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--json", action="store_true",
                    help="emit JSON instead of human-readable lines")
    args = ap.parse_args(argv[1:])

    inst = find_cubit_install()
    if not inst:
        print("ERROR: no Coreform Cubit install found under C:\\Program Files",
              file=sys.stderr)
        print("  Override: set CUBIT_INSTALL_DIR=<path>", file=sys.stderr)
        return 1

    info = {
        "install_dir": str(inst),
        "bin_dir":     str(find_cubit_bin()),
        "cmake_dir":   str(find_cubit_cmake()),
        "exe":         str(find_cubit_exe() or ""),
        "version":     get_cubit_version(),
    }
    if args.json:
        import json
        print(json.dumps(info, indent=2))
    else:
        for k, v in info.items():
            print(f"{k} = {v}")
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv))
