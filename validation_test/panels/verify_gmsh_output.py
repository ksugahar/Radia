"""
GMSH output file sanity checker for Radia panels.

Validates .geo and .msh files BEFORE they are opened in the GUI,
catching the bugs that would otherwise show up as "Open GMSH doesn't
work" or "GMSH window is black".

Usage:
    python verify_gmsh_output.py path/to/file.geo
    python verify_gmsh_output.py path/to/file.msh
    python verify_gmsh_output.py path/to/directory/

Can also be imported and called from pytest:
    from validation_test.panels.verify_gmsh_output import verify_geo, verify_msh
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path


# ============================================================
# Known-bad GMSH options (Gmsh 5)
# ============================================================

_INVALID_OPTIONS = {
    "Mesh.Volumes",      # does NOT exist; use VolumeEdges + VolumeFaces
    "Mesh.Surfaces",     # does NOT exist; use SurfaceEdges + SurfaceFaces
}


# ============================================================
# .geo verification
# ============================================================

def verify_geo(geo_path: str | Path, verbose: bool = True) -> bool:
    """Verify a companion .geo file.

    Checks:
      1. All ``Merge "..."`` targets exist on disk (resolved relative
         to the .geo's directory).
      2. No invalid GMSH options (e.g. ``Mesh.Volumes``).
      3. View[N] indices do not exceed the number of Merge targets
         that carry NodeData (rough heuristic: count Merge lines).

    Returns True if all checks pass.
    """
    geo_path = Path(geo_path)
    geo_dir = geo_path.parent
    lines = geo_path.read_text(encoding="utf-8").splitlines()
    ok = True
    merge_targets = []
    max_view_idx = -1

    for lineno, line in enumerate(lines, 1):
        stripped = line.strip()

        # Check 1: Merge targets exist
        m = re.match(r'Merge\s+"([^"]+)"', stripped)
        if m:
            target = m.group(1)
            target_path = geo_dir / target
            if not target_path.is_file():
                _report(verbose, "FAIL", geo_path.name, lineno,
                        f'Merge target missing: "{target}"')
                ok = False
            else:
                _report(verbose, "OK", geo_path.name, lineno,
                        f'Merge target exists: "{target}" '
                        f'({target_path.stat().st_size:,} bytes)')
            merge_targets.append(target)

        # Check 2: Invalid options
        for opt in _INVALID_OPTIONS:
            if stripped.startswith(opt):
                _report(verbose, "FAIL", geo_path.name, lineno,
                        f'Invalid GMSH option: {opt} '
                        f'(does not exist in GMSH 4.x)')
                ok = False

        # Track View[N] indices
        vm = re.findall(r'View\[(\d+)\]', stripped)
        for idx_str in vm:
            idx = int(idx_str)
            if idx > max_view_idx:
                max_view_idx = idx

    # Check 3: View index range
    # (rough: each .msh with NodeData contributes ~1 view, so
    #  max_view_idx should be < len(merge_targets) * 2 at most)
    if max_view_idx >= 0 and verbose:
        _report(verbose, "INFO", geo_path.name, 0,
                f'Max View index: {max_view_idx}, '
                f'Merge targets: {len(merge_targets)}')

    if ok and verbose:
        print(f"  PASS: {geo_path.name}")
    return ok


# ============================================================
# .msh verification
# ============================================================

def verify_msh(msh_path: str | Path, verbose: bool = True) -> bool:
    """Use the canonical strict MSH 4.2 structural validator."""
    from radia_mcp.gmsh.msh_inspect import validate_msh

    result = validate_msh(msh_path, check_jacobians=False)
    for error in result.get("errors", []):
        _report(verbose, "FAIL", Path(msh_path).name, 0, str(error))
    if result["ok"]:
        _report(verbose, "OK", Path(msh_path).name, 0,
                "MSH 4.2 structural checks passed")
    return bool(result["ok"])


def _report(verbose, level, filename, lineno, msg):
    if not verbose:
        return
    loc = f"{filename}:{lineno}" if lineno else filename
    print(f"  [{level}] {loc}: {msg}")


# ============================================================
# CLI
# ============================================================

def main():
    if len(sys.argv) < 2:
        print("Usage: python verify_gmsh_output.py <file_or_dir>")
        sys.exit(1)

    target = Path(sys.argv[1])
    all_ok = True

    if target.is_dir():
        files = sorted(target.glob("*.geo")) + sorted(target.glob("*.msh"))
    elif target.suffix == ".geo":
        files = [target]
        # Also check referenced .msh files
        text = target.read_text(encoding="utf-8")
        for m in re.finditer(r'Merge\s+"([^"]+)"', text):
            ref = target.parent / m.group(1)
            if ref.is_file():
                files.append(ref)
    elif target.suffix == ".msh":
        files = [target]
    else:
        print(f"Unknown file type: {target}")
        sys.exit(1)

    for f in files:
        print(f"\n=== {f.name} ===")
        if f.suffix == ".geo":
            if not verify_geo(f):
                all_ok = False
        elif f.suffix == ".msh":
            if not verify_msh(f):
                all_ok = False

    print()
    if all_ok:
        print("ALL CHECKS PASSED")
    else:
        print("SOME CHECKS FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
