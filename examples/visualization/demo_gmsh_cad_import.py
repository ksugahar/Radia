#!/usr/bin/env python
"""CAD import workflow for Radia meshes with optional GMSH display.

Radia policy:
  * generate analysis meshes with Netgen/OCC or Cubit/Coreform
  * store NGSolve inputs as Netgen .vol files
  * use GMSH only to inspect existing CAD, .geo, or .msh visualization files

This example keeps the historical filename, but no longer uses the GMSH Python
API for mesh generation.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def build_vol_from_cad(step_file: str | Path,
                       vol_file: str | Path | None = None,
                       maxh: float = 0.01,
                       curve_order: int = 2) -> dict | None:
    """Import a STEP file with Netgen/OCC and save a .vol mesh.

    The heavy imports stay inside the function so the demo remains runnable on
    machines that only want to read the workflow text.
    """
    step_path = Path(step_file)
    if not step_path.exists():
        print(f"[INFO] STEP file not found: {step_path}")
        return None

    from netgen.occ import OCCGeometry
    from ngsolve import Mesh, TaskManager

    out_path = Path(vol_file) if vol_file else step_path.with_suffix(".vol")
    print("=" * 60)
    print("CAD -> Netgen .vol mesh")
    print("=" * 60)
    print(f"STEP: {step_path}")
    print(f"maxh: {maxh}")

    geo = OCCGeometry(str(step_path))
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)
        ngmesh.Save(str(out_path))
        mesh = Mesh(str(out_path))
        if curve_order > 1:
            mesh.Curve(curve_order)

    print(f"[OK] wrote {out_path}")
    print(f"     elements={mesh.ne}, vertices={mesh.nv}, curve_order={curve_order}")
    return {
        "step": str(step_path),
        "vol": str(out_path),
        "elements": mesh.ne,
        "vertices": mesh.nv,
        "curve_order": curve_order,
    }


def open_for_gmsh_display(path: str | Path) -> None:
    """Open an existing file in standalone GMSH for display only."""
    display_path = Path(path)
    if not display_path.exists():
        print(f"[INFO] display file not found: {display_path}")
        return
    subprocess.run(["gmsh", str(display_path)], check=False)


def workflow_example() -> None:
    print("=" * 60)
    print("Radia CAD import workflow")
    print("=" * 60)
    print()
    print("Recommended analysis path:")
    print("  CAD STEP -> Netgen/OCC or Cubit/Coreform -> .vol -> NGSolve/Radia")
    print()
    print("Recommended visualization path:")
    print("  existing STEP/.msh/.geo -> standalone gmsh viewer")
    print()
    print("Why .vol is preferred for solves:")
    print("  [OK] direct NGSolve Mesh('model.vol') input")
    print("  [OK] material and boundary labels stay explicit")
    print("  [OK] avoids pip-gmsh runtime dependency in public examples")
    print("  [OK] keeps GMSH as visualization/post-processing only")


def comparison_table() -> None:
    print()
    print("=" * 60)
    print("Mesh workflow roles")
    print("=" * 60)
    print("| Tool | Role in Radia examples |")
    print("|------|------------------------|")
    print("| Netgen/OCC | tetrahedral .vol generation for NGSolve |")
    print("| Cubit/Coreform | high-quality tet/hex .vol export |")
    print("| GMSH | display existing STEP/.geo/.msh outputs |")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", default="motor_rotor.step",
                        help="STEP file to mesh with Netgen/OCC")
    parser.add_argument("--vol", default=None,
                        help="optional output .vol path")
    parser.add_argument("--maxh", type=float, default=0.01,
                        help="Netgen maximum element size")
    parser.add_argument("--curve-order", type=int, default=2,
                        help="NGSolve mesh.Curve order for curved display/solve")
    parser.add_argument("--open-gmsh", action="store_true",
                        help="open the STEP file in standalone GMSH for display")
    args = parser.parse_args()

    workflow_example()
    comparison_table()

    result = build_vol_from_cad(args.step, args.vol, args.maxh, args.curve_order)
    if args.open_gmsh:
        open_for_gmsh_display(args.step)

    if result is None:
        print()
        print("To run with real CAD:")
        print("  python demo_gmsh_cad_import.py --step path/to/model.step")


if __name__ == "__main__":
    main()
