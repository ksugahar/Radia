#!/usr/bin/env python
"""PEEC conductor boundary surface from a Netgen .vol mesh.

This directory still contains historical GMSH fixtures, but new Radia examples
should not use the GMSH Python API for analysis mesh generation.  PEEC/BEM
workflows can use the boundary triangles of a Netgen/Cubit .vol mesh as the
surface-current support.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def build_rectangular_torus_boundary(radius: float = 0.05,
                                     width: float = 0.002,
                                     height: float = 0.002,
                                     maxh: float = 0.002,
                                     curve_order: int = 2,
                                     vol_file: str | Path | None = None) -> dict:
    from netgen.occ import WorkPlane, Axes, Axis, Pnt, Dir, OCCGeometry
    from ngsolve import BND, CF, Integrate, Mesh, TaskManager

    wp = WorkPlane(Axes(p=Pnt(radius, 0, 0),
                        n=Dir(0, 1, 0),
                        h=Dir(0, 0, 1)))
    section = wp.Rectangle(width, height).Face()
    conductor = section.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360)
    conductor.name = "conductor_volume"
    for face in conductor.faces:
        face.name = "conductor"
        face.maxh = maxh

    geo = OCCGeometry(conductor)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)
        if vol_file:
            ngmesh.Save(str(vol_file))
        mesh = Mesh(ngmesh)
        if curve_order > 1:
            mesh.Curve(curve_order)
        area = Integrate(CF(1.0), mesh, BND)

    exact = 2.0 * (width + height) * 2.0 * math.pi * radius
    return {
        "volume_elements": int(mesh.ne),
        "boundary_elements": int(mesh.GetNE(BND)),
        "vertices": int(mesh.nv),
        "boundary_area": float(area),
        "analytic_area": float(exact),
        "relative_error": float(abs(area - exact) / exact),
        "vol_file": str(vol_file) if vol_file else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-vol", default=None,
                        help="optional .vol output path")
    parser.add_argument("--maxh-mm", type=float, default=2.0,
                        help="mesh size in mm")
    args = parser.parse_args()

    print()
    print("PEEC conductor boundary surface from .vol")
    print("=" * 60)
    print("Key concept:")
    print("  PEEC surface unknowns live on boundary triangles.")
    print("  The conductor volume can still be meshed as .vol for robust labels.")
    print("  GMSH is optional display only, not the mesh generator.")
    print()

    result = build_rectangular_torus_boundary(
        maxh=args.maxh_mm / 1000.0,
        vol_file=args.write_vol,
    )
    for key in ("volume_elements", "boundary_elements", "vertices"):
        print(f"{key:>18s}: {result[key]}")
    print(f"{'boundary_area':>18s}: {result['boundary_area']:.8e}")
    print(f"{'analytic_area':>18s}: {result['analytic_area']:.8e}")
    print(f"{'relative_error':>18s}: {result['relative_error']:.3e}")
    if result["vol_file"]:
        print(f"{'vol_file':>18s}: {result['vol_file']}")

    if result["relative_error"] > 0.05:
        raise SystemExit("FAIL: boundary area error exceeds 5%")
    print()
    print("Overall: PASS")


if __name__ == "__main__":
    main()
