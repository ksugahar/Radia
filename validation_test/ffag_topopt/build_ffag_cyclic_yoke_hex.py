"""Build canonical curved-Q2 HEX FFAG cyclic-yoke meshes with Cubit 2025.12.

The sector and explicit full ring are both ACIS annular sweeps.  Cubit owns
the curved geometry, mapped HEX mesh, boundary sidesets, and order-2 Netgen
export.  The ordinary project Python then adds only Netgen's rotational point
identifications to the named sector cuts; no geometry is reconstructed.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np


def _load_axisymmetric_builder():
    path = Path(__file__).with_name("build_ffag_axisymmetric_c_yoke_hex.py")
    spec = importlib.util.spec_from_file_location("_ffag_cubit_runtime", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _create_annular_sector(cubit, *, start_angle_deg, stop_angle_deg,
                           inner_radius, outer_radius, half_height):
    def arc(radius):
        cubit.cmd(
            "create curve arc radius %.17g center location 0 0 0 "
            "normal 0 0 1 start angle %.17g stop angle %.17g"
            % (radius, start_angle_deg, stop_angle_deg))
        curve = int(cubit.get_last_id("curve"))
        return curve, tuple(int(value) for value in
                            cubit.get_relatives("curve", curve, "vertex"))

    inner_curve, inner_vertices = arc(inner_radius)
    outer_curve, outer_vertices = arc(outer_radius)
    radial_curves = []
    for inner_vertex in inner_vertices:
        point = np.asarray(
            cubit.get_center_point("vertex", inner_vertex), dtype=float)
        direction = point[:2] / np.linalg.norm(point[:2])

        def mismatch(outer_vertex, direction=direction):
            outer = np.asarray(
                cubit.get_center_point("vertex", outer_vertex), dtype=float)
            return np.linalg.norm(outer[:2] / np.linalg.norm(outer[:2])
                                  - direction)

        outer_vertex = min(outer_vertices, key=mismatch)
        cubit.cmd(f"create curve vertex {inner_vertex} {outer_vertex}")
        radial_curves.append(int(cubit.get_last_id("curve")))
    cubit.cmd(
        "create surface curve %d %d %d %d"
        % (inner_curve, outer_curve, radial_curves[0], radial_curves[1]))
    source = int(cubit.get_last_id("surface"))
    cubit.cmd(f"move surface {source} z {-half_height:.17g}")
    cubit.cmd(f"sweep surface {source} vector 0 0 {2.0*half_height:.17g}")
    return int(cubit.get_last_id("volume"))


def _classify_sector_surfaces(cubit, *, fold, inner_radius, outer_radius,
                              half_height):
    theta = 2.0 * math.pi / fold
    middle_radius = 0.5 * (inner_radius + outer_radius)
    scale = max(outer_radius, half_height, 1.0)
    tolerance = 1.0e-8 * scale
    periodic = []
    skin = []
    for surface in cubit.get_entities("surface"):
        center = np.asarray(
            cubit.get_center_point("surface", surface), dtype=float)
        radius = float(np.linalg.norm(center[:2]))
        if (abs(center[2]) <= tolerance
                and abs(radius - middle_radius) <= tolerance):
            angle = math.atan2(center[1], center[0]) % (2.0 * math.pi)
            distance_min = min(angle, 2.0 * math.pi-angle)
            distance_max = abs(angle-theta)
            periodic.append((distance_min, distance_max, int(surface)))
        else:
            skin.append(int(surface))
    if len(periodic) != 2:
        raise RuntimeError(
            "Cubit sector does not expose exactly two azimuthal cut surfaces")
    minimum = min(periodic, key=lambda row: row[0])[2]
    maximum = min(periodic, key=lambda row: row[1])[2]
    if minimum == maximum:
        raise RuntimeError("Cubit cyclic cut-surface classification collapsed")
    return minimum, maximum, tuple(skin)


def _mesh_and_export(cubit, *, output, fold, intervals, inner_radius,
                     outer_radius, half_height, full_ring):
    cubit.cmd("reset")
    theta_deg = 360.0 / fold
    sectors = range(fold) if full_ring else range(1)
    for sector in sectors:
        _create_annular_sector(
            cubit, start_angle_deg=sector*theta_deg,
            stop_angle_deg=(sector+1)*theta_deg,
            inner_radius=inner_radius, outer_radius=outer_radius,
            half_height=half_height)
    if full_ring:
        cubit.cmd("imprint volume all")
        cubit.cmd("merge volume all")
    cubit.cmd("volume all scheme sweep")
    cubit.cmd(f"curve all interval {int(intervals)}")
    cubit.cmd("mesh volume all")
    cubit.cmd("block 1 add hex all")
    cubit.cmd('block 1 name "yoke"')
    if full_ring:
        skin = tuple(
            int(surface) for surface in cubit.get_entities("surface")
            if len(cubit.get_relatives("surface", surface, "volume")) == 1)
        if not skin:
            raise RuntimeError("Cubit full ring has no exterior skin surfaces")
        cubit.cmd("sideset 22 add surface " + " ".join(map(str, skin)))
        cubit.cmd('sideset 22 name "skin"')
        periodic_surfaces = ()
    else:
        minimum, maximum, skin = _classify_sector_surfaces(
            cubit, fold=fold, inner_radius=inner_radius,
            outer_radius=outer_radius, half_height=half_height)
        cubit.cmd(f"sideset 20 add surface {minimum}")
        cubit.cmd('sideset 20 name "periodic_min"')
        cubit.cmd(f"sideset 21 add surface {maximum}")
        cubit.cmd('sideset 21 name "periodic_max"')
        cubit.cmd("sideset 22 add surface " + " ".join(map(str, skin)))
        cubit.cmd('sideset 22 name "skin"')
        periodic_surfaces = (minimum, maximum)
    output.parent.mkdir(parents=True, exist_ok=True)
    cubit.cmd(f'export netgen "{output.as_posix()}" order 2 overwrite')
    cubit.cmd(
        f'save cub5 "{output.with_suffix(".cub5").as_posix()}" overwrite')
    if (cubit.get_hex_count() <= 0 or not output.is_file()
            or output.stat().st_size == 0):
        raise RuntimeError("Cubit did not export a nonempty curved HEX mesh")
    return {
        "output": str(output),
        "full_ring": bool(full_ring),
        "hexes": int(cubit.get_hex_count()),
        "nodes": int(cubit.get_node_count()),
        "periodic_surface_ids": list(periodic_surfaces),
        "bytes": int(output.stat().st_size),
    }


def _child(args):
    runtime = _load_axisymmetric_builder()
    cubit = runtime._load_cubit(runtime._find_cubit_bin(args.cubit_bin))
    sector = _mesh_and_export(
        cubit, output=args.sector_output.resolve(), fold=args.fold,
        intervals=args.intervals, inner_radius=args.inner_radius,
        outer_radius=args.outer_radius, half_height=args.half_height,
        full_ring=False)
    full = _mesh_and_export(
        cubit, output=args.full_output.resolve(), fold=args.fold,
        intervals=args.intervals, inner_radius=args.inner_radius,
        outer_radius=args.outer_radius, half_height=args.half_height,
        full_ring=True)
    print("__RADIA_CUBIT_REPORT__" + json.dumps({
        "cubit_version": "2025.12",
        "sector": sector,
        "full_ring": full,
    }))


def _postprocess(args, child_report):
    import ngsolve as ng
    from cubit_mesh_export.check import check_consistency
    from radia.ffag_topopt import (
        build_ffag_cyclic_density_map,
        identify_ffag_cyclic_sector_vertices,
    )

    sector_path = args.sector_output.resolve()
    full_path = args.full_output.resolve()
    sector = ng.Mesh(str(sector_path))
    identification = identify_ffag_cyclic_sector_vertices(sector, args.fold)
    sector.ngmesh.Save(str(sector_path))
    sector = ng.Mesh(str(sector_path))
    full = ng.Mesh(str(full_path))
    sector_check = check_consistency(
        sector_path, min_curve_order=2, required_materials=("yoke",),
        required_boundaries=("skin", "periodic_min", "periodic_max"))
    full_check = check_consistency(
        full_path, min_curve_order=2, required_materials=("yoke",),
        required_boundaries=("skin",))
    if not sector_check["passed"] or not full_check["passed"]:
        raise RuntimeError("check-vol rejected the Cubit FFAG mesh pair")
    if ({len(element.vertices) for element in sector.Elements(ng.VOL)} != {8}
            or {len(element.vertices) for element in full.Elements(ng.VOL)} != {8}
            or sector.GetCurveOrder() != 2 or full.GetCurveOrder() != 2):
        raise RuntimeError("FFAG Cubit lane requires pure curved-Q2 HEX")
    density_map = build_ffag_cyclic_density_map(sector)
    with ng.TaskManager():
        sector_volume = float(ng.Integrate(ng.CF(1.0), sector))
        full_volume = float(ng.Integrate(ng.CF(1.0), full))
    exact_sector = ((2.0*math.pi/args.fold)
                    * (args.outer_radius**2-args.inner_radius**2)
                    * args.half_height)
    exact_full = (2.0*math.pi
                  * (args.outer_radius**2-args.inner_radius**2)
                  * args.half_height)
    report = {
        "schema": "radia.ffag-cubit-cyclic-yoke-hex/v1",
        "cubit": child_report,
        "geometry": {
            "fold": int(args.fold),
            "inner_radius_m": float(args.inner_radius),
            "outer_radius_m": float(args.outer_radius),
            "half_height_m": float(args.half_height),
            "curve_order": 2,
            "intervals_per_curve": int(args.intervals),
        },
        "sector": {
            "elements": int(sector.ne),
            "vertices": int(sector.nv),
            "volume_m3": sector_volume,
            "exact_volume_m3": exact_sector,
            "relative_volume_error": abs(sector_volume-exact_sector)/exact_sector,
            "periodic_identification": identification,
            "density_variables": int(density_map.variable_count),
            "density_elements": int(density_map.element_count),
            "density_boundary_facet_pairs": int(
                density_map.boundary_pair_count),
            "check_vol_passed": bool(sector_check["passed"]),
        },
        "full_ring": {
            "elements": int(full.ne),
            "vertices": int(full.nv),
            "volume_m3": full_volume,
            "exact_volume_m3": exact_full,
            "relative_volume_error": abs(full_volume-exact_full)/exact_full,
            "check_vol_passed": bool(full_check["passed"]),
        },
    }
    report["checks"] = {
        "cubit_sector_is_curved_q2_pure_hex": True,
        "cubit_full_ring_is_curved_q2_pure_hex": True,
        "periodic_vertices_identified": identification["pair_count"] > 0,
        "periodic_density_variables_are_reduced": (
            density_map.variable_count < density_map.element_count),
        "sector_volume_matches_full_ring_share": abs(
            args.fold*sector_volume-full_volume)/exact_full < 2.0e-5,
    }
    report["passed"] = all(report["checks"].values())
    output = args.report.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise RuntimeError("Cubit FFAG cyclic-yoke build gate failed")


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    root = Path(r"C:\temp\radia_ffag_cyclic")
    parser.add_argument("--sector-output", type=Path,
                        default=root / "ffag_cyclic_sector_q2.vol")
    parser.add_argument("--full-output", type=Path,
                        default=root / "ffag_cyclic_full_q2.vol")
    parser.add_argument("--report", type=Path,
                        default=root / "ffag_cyclic_mesh_build.json")
    parser.add_argument("--fold", type=int, default=12)
    parser.add_argument("--inner-radius", type=float, default=1.0)
    parser.add_argument("--outer-radius", type=float, default=2.0)
    parser.add_argument("--half-height", type=float, default=0.25)
    parser.add_argument("--intervals", type=int, default=2)
    parser.add_argument("--cubit-bin", type=Path)
    parser.add_argument("--cubit-child", action="store_true",
                        help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def main():
    args = parse_args()
    if (args.fold < 3 or args.intervals < 2
            or not 0.0 < args.inner_radius < args.outer_radius
            or args.half_height <= 0.0):
        raise ValueError("invalid FFAG cyclic Cubit mesh settings")
    if args.cubit_child:
        _child(args)
        return 0
    for path in (args.sector_output, args.full_output,
                 args.sector_output.with_suffix(".cub5"),
                 args.full_output.with_suffix(".cub5"), args.report):
        path.resolve().unlink(missing_ok=True)
    command = [sys.executable, str(Path(__file__).resolve()),
               "--cubit-child",
               "--sector-output", str(args.sector_output.resolve()),
               "--full-output", str(args.full_output.resolve()),
               "--report", str(args.report.resolve()),
               "--fold", str(args.fold),
               "--inner-radius", str(args.inner_radius),
               "--outer-radius", str(args.outer_radius),
               "--half-height", str(args.half_height),
               "--intervals", str(args.intervals)]
    if args.cubit_bin is not None:
        command.extend(("--cubit-bin", str(args.cubit_bin.resolve())))
    completed = subprocess.run(
        command, capture_output=True, text=True, errors="replace", check=False)
    report = None
    for line in completed.stdout.splitlines():
        if line.startswith("__RADIA_CUBIT_REPORT__"):
            report = json.loads(line.removeprefix("__RADIA_CUBIT_REPORT__"))
    if completed.returncode or report is None:
        tail = "\n".join((completed.stdout+completed.stderr).splitlines()[-80:])
        raise RuntimeError(
            "Cubit FFAG cyclic mesh build failed "
            f"(exit={completed.returncode}):\n{tail}")
    _postprocess(args, report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
