#!/usr/bin/env python
"""Circular coil boundary-surface demo for PEEC-style workflows.

The filename is historical.  The current Radia policy is:
  * create analysis meshes with Netgen/OCC or Cubit/Coreform
  * save NGSolve inputs as .vol
  * use GMSH only to display existing .msh/.geo outputs

This example builds a rectangular-section toroidal coil with Netgen/OCC,
reports its boundary surface area, and optionally writes a .vol file.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def rectangular_torus_surface_area(mean_radius: float,
                                   wire_width: float,
                                   wire_height: float) -> float:
    """Analytic surface area for a swept rectangular section."""
    section_perimeter = 2.0 * (wire_width + wire_height)
    centerline_length = 2.0 * math.pi * mean_radius
    return section_perimeter * centerline_length


def create_circular_coil_vol(mean_radius: float = 0.05,
                             wire_width: float = 0.004,
                             wire_height: float = 0.004,
                             maxh: float = 0.003,
                             curve_order: int = 2,
                             vol_file: str | Path | None = None) -> dict:
    """Build a rectangular-section toroidal coil and inspect its boundary."""
    from netgen.occ import WorkPlane, Axes, Axis, Pnt, Dir, OCCGeometry
    from ngsolve import BND, CF, Integrate, Mesh, TaskManager

    print("=" * 60)
    print("Circular coil via Netgen/OCC")
    print("=" * 60)
    print(f"Mean radius: {mean_radius * 1000:.1f} mm")
    print(f"Wire width:  {wire_width * 1000:.2f} mm")
    print(f"Wire height: {wire_height * 1000:.2f} mm")
    print(f"maxh:        {maxh * 1000:.2f} mm")

    wp = WorkPlane(Axes(p=Pnt(mean_radius, 0, 0),
                        n=Dir(0, 1, 0),
                        h=Dir(0, 0, 1)))
    section = wp.Rectangle(wire_width, wire_height).Face()
    coil = section.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360)
    coil.name = "coil"
    for face in coil.faces:
        face.name = "conductor"
        face.maxh = maxh

    geo = OCCGeometry(coil)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)
        if vol_file:
            ngmesh.Save(str(vol_file))
        mesh = Mesh(ngmesh)
        if curve_order > 1:
            mesh.Curve(curve_order)
        boundary_area = Integrate(CF(1.0), mesh, BND)

    exact_area = rectangular_torus_surface_area(
        mean_radius, wire_width, wire_height)
    rel_error = abs(boundary_area - exact_area) / exact_area
    boundary_elements = mesh.GetNE(BND)

    print()
    print("Boundary surface summary")
    print("-" * 60)
    print(f"Volume elements:   {mesh.ne}")
    print(f"Boundary elements: {boundary_elements}")
    print(f"Vertices:          {mesh.nv}")
    print(f"Boundary area:     {boundary_area:.8e} m^2")
    print(f"Analytic area:     {exact_area:.8e} m^2")
    print(f"Relative error:    {rel_error:.3e}")
    if vol_file:
        print(f"Saved .vol:        {vol_file}")

    return {
        "vol_file": str(vol_file) if vol_file else None,
        "volume_elements": int(mesh.ne),
        "boundary_elements": int(boundary_elements),
        "vertices": int(mesh.nv),
        "boundary_area": float(boundary_area),
        "analytic_area": float(exact_area),
        "relative_error": float(rel_error),
        "curve_order": int(curve_order),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=float, default=50.0,
                        help="mean radius in mm")
    parser.add_argument("--width", type=float, default=4.0,
                        help="wire width in mm")
    parser.add_argument("--height", type=float, default=4.0,
                        help="wire height in mm")
    parser.add_argument("--maxh", type=float, default=3.0,
                        help="mesh size in mm")
    parser.add_argument("--curve-order", type=int, default=2,
                        help="NGSolve mesh.Curve order")
    parser.add_argument("--write-vol", default=None,
                        help="optional .vol output path")
    args = parser.parse_args()

    result = create_circular_coil_vol(
        mean_radius=args.radius / 1000.0,
        wire_width=args.width / 1000.0,
        wire_height=args.height / 1000.0,
        maxh=args.maxh / 1000.0,
        curve_order=args.curve_order,
        vol_file=args.write_vol,
    )

    if result["relative_error"] > 0.05:
        raise SystemExit("FAIL: boundary area error exceeds 5%")

    print()
    print("Overall: PASS")


if __name__ == "__main__":
    main()
