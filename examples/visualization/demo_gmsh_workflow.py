#!/usr/bin/env python
"""
Netgen -> NGSolve -> Radia -> GMSH display workflow.

This example keeps mesh generation in Netgen/NGSolve and uses GMSH only as a
file-based viewer for exported mesh and post-processing data.

Requirements:
- netgen / ngsolve
- radia
- optional standalone gmsh executable for --open-gmsh
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import radia as rad


BOX_MIN = (-0.05, -0.05, -0.05)
BOX_MAX = (0.05, 0.05, 0.05)


def create_box_mesh_with_netgen(maxh: float):
    """Create the demonstration box mesh with Netgen/OCC."""
    from netgen.occ import Box, OCCGeometry, Pnt
    from ngsolve import Mesh

    pmin = Pnt(*BOX_MIN)
    pmax = Pnt(*BOX_MAX)
    shape = Box(pmin, pmax)
    ngmesh = OCCGeometry(shape).GenerateMesh(maxh=maxh)
    return Mesh(ngmesh)


def validate_box_mesh(mesh) -> dict[str, float]:
    """Validate volume and boundary area against the analytic box values."""
    from ngsolve import BND, CF, Integrate

    dx = BOX_MAX[0] - BOX_MIN[0]
    dy = BOX_MAX[1] - BOX_MIN[1]
    dz = BOX_MAX[2] - BOX_MIN[2]
    exact_volume = dx * dy * dz
    exact_area = 2.0 * (dx * dy + dx * dz + dy * dz)

    volume = float(Integrate(CF(1.0), mesh))
    area = float(Integrate(CF(1.0), mesh, BND))
    return {
        "volume": volume,
        "exact_volume": exact_volume,
        "volume_rel_error": abs(volume - exact_volume) / exact_volume,
        "surface_area": area,
        "exact_surface_area": exact_area,
        "surface_area_rel_error": abs(area - exact_area) / exact_area,
    }


def write_gmsh_msh22(mesh, path: Path) -> dict[str, int]:
    """Write a Netgen/NGSolve tri/tet mesh as GMSH v2.2 ASCII for display."""
    from ngsolve import BND, VOL

    vertices = sorted(list(mesh.vertices), key=lambda vertex: vertex.nr)
    boundary = list(mesh.Elements(BND))
    volume = list(mesh.Elements(VOL))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        fh.write("$Nodes\n")
        fh.write(f"{len(vertices)}\n")
        for vertex in vertices:
            x, y, z = vertex.point
            fh.write(f"{vertex.nr + 1} {x:.16g} {y:.16g} {z:.16g}\n")
        fh.write("$EndNodes\n")

        fh.write("$Elements\n")
        fh.write(f"{len(boundary) + len(volume)}\n")
        element_id = 1
        for element in boundary:
            node_ids = [vertex.nr + 1 for vertex in element.vertices]
            if len(node_ids) != 3:
                raise ValueError("Expected triangular boundary elements")
            fh.write(
                f"{element_id} 2 2 2 2 "
                f"{node_ids[0]} {node_ids[1]} {node_ids[2]}\n"
            )
            element_id += 1
        for element in volume:
            node_ids = [vertex.nr + 1 for vertex in element.vertices]
            if len(node_ids) != 4:
                raise ValueError("Expected tetrahedral volume elements")
            fh.write(
                f"{element_id} 4 2 1 1 "
                f"{node_ids[0]} {node_ids[1]} {node_ids[2]} "
                f"{node_ids[3]}\n"
            )
            element_id += 1
        fh.write("$EndElements\n")

    return {
        "nodes": len(vertices),
        "boundary_triangles": len(boundary),
        "tetrahedra": len(volume),
    }


def compute_radia_field_samples(mesh, sample_count: int) -> list[tuple[tuple[float, float, float], np.ndarray]]:
    """Compute Radia B-field vectors at representative mesh vertices."""
    rad.UtiDelAll()
    magnet = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 954930])

    vertices = sorted(list(mesh.vertices), key=lambda vertex: vertex.nr)
    step = max(1, math.floor(len(vertices) / max(1, sample_count)))
    selected = vertices[::step][:sample_count]

    samples = []
    for vertex in selected:
        point = tuple(float(v) for v in vertex.point)
        field = np.asarray(rad.Fld(magnet, "b", point), dtype=float)
        samples.append((point, field))
    return samples


def write_gmsh_pos_vectors(
    samples: list[tuple[tuple[float, float, float], np.ndarray]],
    path: Path,
) -> None:
    """Write sampled Radia vectors as a GMSH .pos post-processing view."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write('View "Radia B-field samples" {\n')
        for point, field in samples:
            x, y, z = point
            bx, by, bz = field
            fh.write(
                f"VP({x:.16g},{y:.16g},{z:.16g})"
                f"{{{bx:.16g},{by:.16g},{bz:.16g}}};\n"
            )
        fh.write("};\n")


def write_gmsh_display_companion(msh_path: Path, pos_path: Path, geo_path: Path) -> None:
    """Write a .geo companion that opens the exported mesh and vector view."""
    msh_name = msh_path.name.replace("\\", "/")
    pos_name = pos_path.name.replace("\\", "/")
    geo_path.parent.mkdir(parents=True, exist_ok=True)
    with geo_path.open("w", encoding="ascii", newline="\n") as fh:
        fh.write("// Open with: gmsh netgen_box_display.geo\n")
        fh.write("Mesh.NumSubEdges = 4;\n")
        fh.write(f'Merge "{msh_name}";\n')
        fh.write(f'Merge "{pos_name}";\n')
        fh.write("General.BackgroundGradient = 0;\n")
        fh.write("General.Trackball = 0;\n")
        fh.write("View[0].VectorType = 5;\n")
        fh.write("View[0].ShowScale = 0;\n")


def maybe_open_gmsh(display_geo: Path) -> None:
    """Open the display companion in standalone GMSH when explicitly requested."""
    gmsh = shutil.which("gmsh")
    if gmsh is None:
        print("  [SKIP] standalone gmsh executable was not found on PATH")
        return
    subprocess.run([gmsh, str(display_geo)], check=False)


def run_workflow(maxh: float, out_dir: Path, sample_count: int) -> dict[str, object]:
    print("=" * 60)
    print("Netgen/NGSolve mesh with GMSH display files")
    print("=" * 60)

    print("\n[1] Generating mesh with Netgen/OCC...")
    mesh = create_box_mesh_with_netgen(maxh=maxh)
    quality = validate_box_mesh(mesh)
    print(f"    Vertices: {mesh.nv}")
    print(f"    Tetrahedra: {mesh.ne}")
    print(f"    Volume: {quality['volume']:.8e} "
          f"(rel err {quality['volume_rel_error']:.3e})")
    print(f"    Surface area: {quality['surface_area']:.8e} "
          f"(rel err {quality['surface_area_rel_error']:.3e})")

    out_dir.mkdir(parents=True, exist_ok=True)
    vol_path = out_dir / "netgen_box.vol"
    msh_path = out_dir / "netgen_box_for_gmsh.msh"
    pos_path = out_dir / "radia_box_field.pos"
    geo_path = out_dir / "netgen_box_display.geo"

    print("\n[2] Writing mesh files...")
    mesh.ngmesh.Save(str(vol_path))
    gmsh_stats = write_gmsh_msh22(mesh, msh_path)
    print(f"    Netgen .vol: {vol_path}")
    print(f"    GMSH display .msh: {msh_path}")

    print("\n[3] Sampling Radia field and writing GMSH .pos view...")
    samples = compute_radia_field_samples(mesh, sample_count=sample_count)
    write_gmsh_pos_vectors(samples, pos_path)
    write_gmsh_display_companion(msh_path, pos_path, geo_path)
    magnitudes = [float(np.linalg.norm(field)) for _, field in samples]
    print(f"    Samples: {len(samples)}")
    print(f"    |B| range: {min(magnitudes):.6g} T .. {max(magnitudes):.6g} T")
    print(f"    GMSH .pos: {pos_path}")
    print(f"    GMSH display companion: {geo_path}")

    print("\n[4] Viewer commands:")
    print(f"    netgen {vol_path}")
    print(f"    gmsh {geo_path}")
    print("\n[OK] Workflow complete")

    return {
        "quality": quality,
        "gmsh_stats": gmsh_stats,
        "samples": len(samples),
        "out_dir": str(out_dir),
        "vol_path": str(vol_path),
        "msh_path": str(msh_path),
        "pos_path": str(pos_path),
        "geo_path": str(geo_path),
        "b_min": min(magnitudes),
        "b_max": max(magnitudes),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    default_out = Path(tempfile.gettempdir()) / "radia_gmsh_workflow"
    parser.add_argument("--maxh", type=float, default=0.02,
                        help="Netgen maximum mesh size in meters")
    parser.add_argument("--sample-count", type=int, default=8,
                        help="number of mesh vertices used for Radia field samples")
    parser.add_argument("--out-dir", type=Path, default=default_out,
                        help="directory for .vol/.msh/.pos/.geo outputs")
    parser.add_argument("--open-gmsh", action="store_true",
                        help="open the generated display companion in standalone gmsh")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = run_workflow(
        maxh=args.maxh,
        out_dir=args.out_dir,
        sample_count=args.sample_count,
    )
    if args.open_gmsh:
        maybe_open_gmsh(Path(result["geo_path"]))


if __name__ == "__main__":
    main()
