"""Validation-class Netgen .vol surface-closure checks for FEM/BEM intake.

This example is intentionally outside the fast unit-test path. It records the
geometry checks that matter when a Cubit/Coreform or Netgen-style ``.vol`` file
is used as a shared FEM/BEM mesh:

* oriented boundary triangle area vectors should close to zero;
* boundary-triangle signed volume should match the tetrahedral volume in
  absolute value;
* flipping one triangle should be detectable by both checks;
* optional external ``--vol`` input can analyze a freshly exported mesh.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol, read_netgen_tri_tet_vol


POINTS = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
    (0.5, 0.5, 0.5),
]

OUTWARD_SURFACE_TRIANGLES = [
    (1, 3, 2),
    (1, 4, 3),
    (5, 6, 7),
    (5, 7, 8),
    (1, 2, 6),
    (1, 6, 5),
    (4, 7, 3),
    (4, 8, 7),
    (1, 5, 8),
    (1, 8, 4),
    (2, 3, 7),
    (2, 7, 6),
]


def _star_cube_vol(surface_triangles):
    lines = [
        "mesh3d",
        "dimension",
        "3",
        "geomtype",
        "0",
        "facedescriptors",
        "1",
        "1 1 0 1 1",
        "surfaceelements",
        str(len(surface_triangles)),
    ]
    for tri in surface_triangles:
        lines.append(f"1 1 1 0 3 {tri[0]} {tri[1]} {tri[2]}")
    lines.extend(["volumeelements", str(len(surface_triangles))])
    for tri in surface_triangles:
        lines.append(f"1 4 {tri[0]} {tri[1]} {tri[2]} 9")
    lines.extend(["points", str(len(POINTS))])
    lines.extend(f"{x} {y} {z}" for x, y, z in POINTS)
    lines.extend(["pointelements", "0", "materials", "1", "1 air", "bcnames", "1", "1 outer", "endmesh"])
    return "\n".join(lines) + "\n"


def _record(name, mesh):
    closure = mesh.surface_closure_summary()
    return {
        "name": name,
        "mesh_summary": mesh.summary(),
        "bounding_box": mesh.bounding_box(),
        "closure": closure,
        "edge_length_summary": mesh.tetrahedron_edge_length_summary(),
    }


def build_summary(vol_path=None):
    outward = parse_netgen_tri_tet_vol(_star_cube_vol(OUTWARD_SURFACE_TRIANGLES))
    flipped = list(OUTWARD_SURFACE_TRIANGLES)
    flipped[2] = tuple(reversed(flipped[2]))
    one_triangle_flipped = parse_netgen_tri_tet_vol(_star_cube_vol(flipped))

    records = [
        _record("builtin_star_cube_outward", outward),
        _record("builtin_star_cube_one_triangle_flipped", one_triangle_flipped),
    ]
    if vol_path is not None:
        records.append(_record("external_vol", read_netgen_tri_tet_vol(vol_path)))

    clean = records[0]["closure"]
    bad = records[1]["closure"]
    checks = {
        "clean_volume": clean["tetrahedron_total_volume"],
        "clean_surface_area": clean["total_surface_area"],
        "clean_vector_area_norm_over_area": clean["surface_vector_area_norm_over_area"],
        "clean_abs_volume_rel_error": clean["surface_abs_volume_rel_error"],
        "clean_boundary_orientation": clean["boundary_orientation"],
        "flipped_vector_area_norm_over_area": bad["surface_vector_area_norm_over_area"],
        "flipped_abs_volume_rel_error": bad["surface_abs_volume_rel_error"],
    }

    if abs(clean["tetrahedron_total_volume"] - 1.0) > 1e-14:
        raise AssertionError("builtin cube tetra volume drifted")
    if abs(clean["total_surface_area"] - 6.0) > 1e-14:
        raise AssertionError("builtin cube surface area drifted")
    if clean["surface_vector_area_norm_over_area"] > 1e-14:
        raise AssertionError("closed cube surface vector area did not cancel")
    if clean["surface_abs_volume_rel_error"] > 1e-14:
        raise AssertionError("closed cube boundary volume does not match tet volume")
    if bad["surface_vector_area_norm_over_area"] < 0.05:
        raise AssertionError("single flipped triangle was not detected by vector area")
    if bad["surface_abs_volume_rel_error"] < 0.05:
        raise AssertionError("single flipped triangle was not detected by volume mismatch")

    return {
        "problem": "Netgen .vol boundary closure and orientation validation",
        "checks": checks,
        "records": records,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vol", help="Optional external tri/tet Netgen .vol file")
    parser.add_argument(
        "--output",
        default=str(Path(__file__).with_name("validation_vol_surface_closure_summary.json")),
    )
    args = parser.parse_args()

    summary = build_summary(args.vol)
    out = Path(args.output)
    out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary["checks"], indent=2, sort_keys=True))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
