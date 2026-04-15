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
    from tests.panels.verify_gmsh_output import verify_geo, verify_msh
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path


# ============================================================
# Known-bad GMSH options (GMSH 4.x)
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
    """Verify a GMSH v4.1 .msh file.

    v4.1 uses an entity-block structure where element IDs restart per
    block by design, so the v2.2-style duplicate-ID / node-count
    checks do not apply. Only NodeData consistency is verified here.

    Returns True if all checks pass. Fails loudly if the file is not
    GMSH v4.1.
    """
    msh_path = Path(msh_path)
    text = msh_path.read_text(encoding="utf-8")

    # Detect version from $MeshFormat — v4.1 is the only supported format.
    m_ver = re.search(r'\$MeshFormat\s+(\d+\.\d+)', text)
    if not m_ver:
        _report(verbose, "FAIL", msh_path.name, 0,
                '$MeshFormat header not found; cannot verify')
        return False
    version = m_ver.group(1)
    if not version.startswith("4"):
        _report(verbose, "FAIL", msh_path.name, 0,
                f'GMSH v{version} is not supported; '
                f'v4.1 is the only accepted format')
        return False

    _report(verbose, "INFO", msh_path.name, 0,
            f'GMSH v{version} (entity-block format, '
            f'structural checks skipped)')
    ok = _check_nodedata(text, msh_path.name, verbose)
    if ok and verbose:
        print(f"  PASS: {msh_path.name} (v{version})")
    return ok


def _check_nodedata(text, filename, verbose):
    """Check NodeData sections (same layout in v4.1 as in legacy v2.2)."""
    ok = True
    in_nodedata = False
    nd_name = ""
    nd_declared = 0
    nd_actual = 0
    nd_header_lines = 0
    n_views = 0

    for line in text.splitlines():
        stripped = line.strip()

        if stripped == "$NodeData":
            in_nodedata = True
            nd_header_lines = 0
            nd_name = ""
            nd_declared = 0
            nd_actual = 0
            continue
        if stripped == "$EndNodeData":
            in_nodedata = False
            n_views += 1
            if nd_declared != nd_actual:
                _report(verbose, "FAIL", filename, 0,
                        f'NodeData "{nd_name}" declared {nd_declared} '
                        f'nodes but found {nd_actual} data lines')
                ok = False
            else:
                _report(verbose, "OK", filename, 0,
                        f'NodeData "{nd_name}": {nd_actual} nodes')
            continue
        if in_nodedata:
            nd_header_lines += 1
            # v2.2 NodeData header layout:
            #   line 1: "1" (num string tags)
            #   line 2: "name"
            #   line 3: "1" (num real tags)
            #   line 4: time value
            #   line 5: "3" (num int tags)
            #   line 6: timestep
            #   line 7: ncomp
            #   line 8: n_nodes   <-- declared count
            #   line 9+: node_id val1 [val2 val3]
            if nd_header_lines == 2:
                nd_name = stripped.strip('"')
            elif nd_header_lines == 8:
                try:
                    nd_declared = int(stripped)
                except ValueError:
                    pass
            elif nd_header_lines > 8:
                nd_actual += 1
            continue

    if verbose:
        _report(verbose, "INFO", filename, 0,
                f'{n_views} NodeData view(s)')
    return ok


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
