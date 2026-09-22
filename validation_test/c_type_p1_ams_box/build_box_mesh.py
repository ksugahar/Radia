"""Build one finite-box tetrahedral mesh of the ESRF Example-5 C-type dipole.

The iron is the Example-5 C-yoke of ``radia.esrf_examples``; the coil is the
rounded-rectangle solid current the three-engine lane drives, here also
meshed so that a solver which needs a meshed current region receives exactly
the same tetrahedra as the solvers which do not.  Materials: ``iron``,
``coil``, ``air``.  Boundaries: every iron face is ``iron_air_interface``,
every coil face ``coil_air_interface``, the six box faces ``outer``; the box
corner at (-L, -L, -L) is the ``GND`` vertex for scalar-potential point gauges.

Netgen/OCC is the explicit mesh route of this lane (the Kelvin three-engine
lane keeps its Cubit/ACIS authority).  ``--scale`` multiplies every mesh size;
the default sizes give about 195k tetrahedra at 0.75 and 530k at 0.5.  The
mesh is written as Netgen ``.vol`` beside a JSON contract with counts, labels
and the ``.vol`` SHA-256 that the solver lane verifies before running.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE_ROOT = HERE.parents[1] / "src"
if (SOURCE_ROOT / "radia" / "__init__.py").is_file():
    sys.path.insert(0, str(SOURCE_ROOT))

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, Glue, OCCGeometry, Pnt  # noqa: E402

from radia.esrf_examples import build_esrf_coils, build_esrf_occ  # noqa: E402

SCHEMA = "radia.validation.c-type-box-mesh.v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fuse_all(shape):
    solids = list(shape.solids)
    fused = solids[0]
    for solid in solids[1:]:
        fused = fused + solid
    return fused


def build_geometry(*, box_half: float, scale: float, gap_h: float, iron_h: float,
                   coil_h: float, air_h: float, mesh_coil: bool = True):
    shapes = build_esrf_occ(5, include_coils=True)
    iron = fuse_all(shapes["iron"])
    if len(list(iron.solids)) != 1:
        raise RuntimeError(f"iron fuse left {len(list(iron.solids))} solids")
    iron.mat("iron")
    iron.faces.name = "iron_air_interface"
    iron.maxh = iron_h * scale

    # The CoilBuilder solid is a Glue of swept segments that share faces; every
    # piece carries the same material name.  A reduced formulation with an
    # exact source never needs it meshed (``mesh_coil=False``): the coil
    # volume is then plain air and the mesh is what such a solver would use.
    coil_solids = []
    if mesh_coil:
        coil_shape = build_esrf_coils(5)[0].to_occ()
        coil_solids = list(coil_shape.solids)
        for solid in coil_solids:
            solid.mat("coil")
            solid.maxh = coil_h * scale
        coil_shape.faces.name = "coil_air_interface"

    # Gap refinement: an air box spanning the 10 mm gap and the pole faces.
    gap_box = Box(Pnt(-0.030, -0.025, -0.006), Pnt(0.030, 0.025, 0.006))
    gap_box.mat("air")
    gap_box.maxh = gap_h * scale
    gap_air = gap_box - iron

    half = box_half
    outer = Box(Pnt(-half, -half, -half), Pnt(half, half, half))
    outer.faces.name = "outer"
    outer.mat("air")
    outer.maxh = air_h * scale
    air = outer - iron - gap_box
    if mesh_coil:
        air = air - coil_shape
    gauge_named = False
    for vertex in air.vertices:
        p = vertex.p
        if all(abs(p[k] + half) < 1e-9 for k in range(3)):
            vertex.name = "GND"
            gauge_named = True
    if not gauge_named:
        raise RuntimeError("GND corner vertex not found on the air box")
    return OCCGeometry(Glue([iron] + coil_solids + [gap_air, air]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--box-half", type=float, default=0.5)
    parser.add_argument("--scale", type=float, default=1.0,
                        help="global multiplier on every mesh size")
    parser.add_argument("--gap-h", type=float, default=0.0025)
    parser.add_argument("--iron-h", type=float, default=0.008)
    parser.add_argument("--coil-h", type=float, default=0.010)
    parser.add_argument("--air-h", type=float, default=0.12)
    parser.add_argument("--grading", type=float, default=0.3)
    parser.add_argument("--no-coil", action="store_true",
                        help="leave the coil volume as air (reduced formulations with an exact source)")
    options = parser.parse_args()
    if options.scale <= 0.0 or options.box_half <= 0.2:
        raise ValueError("--scale must be positive and --box-half must exceed the magnet")

    started = time.perf_counter()
    geometry = build_geometry(box_half=options.box_half, scale=options.scale,
                              gap_h=options.gap_h, iron_h=options.iron_h,
                              coil_h=options.coil_h, air_h=options.air_h,
                              mesh_coil=not options.no_coil)
    ngmesh = geometry.GenerateMesh(maxh=options.air_h * options.scale,
                                   grading=options.grading)
    build_seconds = time.perf_counter() - started
    mesh = ng.Mesh(ngmesh)
    options.output.parent.mkdir(parents=True, exist_ok=True)
    ngmesh.Save(str(options.output))

    materials = set(mesh.GetMaterials())
    boundaries = set(mesh.GetBoundaries())
    bbboundaries = set(mesh.GetBBBoundaries())
    for required in (("iron", "air") if options.no_coil else ("iron", "coil", "air")):
        if required not in materials:
            raise RuntimeError(f"material {required!r} missing: {sorted(materials)}")
    for required in ("iron_air_interface", "outer"):
        if required not in boundaries:
            raise RuntimeError(f"boundary {required!r} missing: {sorted(boundaries)}")
    if "GND" not in bbboundaries:
        raise RuntimeError(f"GND vertex missing: {sorted(bbboundaries)}")
    counts: dict[str, int] = {}
    for element in mesh.Elements(ng.VOL):
        counts[str(element.mat)] = counts.get(str(element.mat), 0) + 1
    contract = {
        "schema": SCHEMA,
        "geometry": ("ESRF Example 5 C-type dipole; iron + finite air box (coil volume left as air)"
                     if options.no_coil else
                     "ESRF Example 5 C-type dipole; iron + meshed coil + finite air box"),
        "coil_meshed": not options.no_coil,
        "vol": str(options.output.resolve()),
        "vol_sha256": sha256(options.output),
        "box_half_m": options.box_half,
        "scale": options.scale,
        "mesh_sizes_m": {"gap": options.gap_h * options.scale,
                         "iron": options.iron_h * options.scale,
                         "coil": options.coil_h * options.scale,
                         "air": options.air_h * options.scale},
        "grading": options.grading,
        "curve_order": 1,
        "elements": int(mesh.ne),
        "vertices": int(mesh.nv),
        "elements_by_material": counts,
        "materials": sorted(materials),
        "boundaries": sorted(boundaries),
        "bbboundaries": sorted(bbboundaries),
        "build_seconds": build_seconds,
        "netgen": ng.__version__,
    }
    options.output.with_suffix(".json").write_text(
        json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(contract, indent=2))


if __name__ == "__main__":
    main()
