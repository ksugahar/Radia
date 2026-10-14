"""Build the iron-only HEX superset for the full-field FFAG PoC.

The fixture is an explicit fully axisymmetric annular C-yoke: two annular pole
packs border the unmeshed beam gap and an outer annular return leg joins them.
It contains no azimuthal sector cut or periodic FE trace.  The source
coils and particle trajectories are deliberately not meshed.  Run this file
with the ordinary project Python; it locates and embeds Coreform Cubit.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np


_DLL_HANDLES = []


def _find_cubit_bin(explicit: Path | None) -> Path:
    candidates = []
    if explicit is not None:
        candidates.append(explicit)
    configured = os.environ.get("CUBIT_PATH")
    if configured:
        path = Path(configured)
        candidates.append(path.parent if path.is_file() else path)
    if os.name == "nt":
        candidates.extend(sorted(
            Path(r"C:\Program Files").glob("Coreform Cubit */bin"),
            reverse=True))
    for candidate in candidates:
        if (candidate / "coreform_cubit.exe").is_file():
            return candidate.resolve()
    raise RuntimeError("Coreform Cubit 2025.12+ was not found; set CUBIT_PATH")


def _load_cubit(cubit_bin: Path):
    sys.path.insert(0, str(cubit_bin))
    os.environ["PATH"] = str(cubit_bin) + os.pathsep + os.environ.get(
        "PATH", "")
    if hasattr(os, "add_dll_directory"):
        for directory in (cubit_bin, cubit_bin / "python3"):
            if directory.is_dir():
                _DLL_HANDLES.append(os.add_dll_directory(str(directory)))
    import cubit

    cubit.init([
        "cubit", "-nojournal", "-batch", "-nographics",
        "-commandplugindir", str(cubit_bin / "plugins")])
    return cubit


def _annulus(cubit, *, inner_radius, outer_radius, height, center_z):
    cubit.cmd(f"create cylinder height {height:.12g} radius {outer_radius:.12g}")
    outer = int(cubit.get_last_id("volume"))
    cubit.cmd(f"create cylinder height {height:.12g} radius {inner_radius:.12g}")
    inner = int(cubit.get_last_id("volume"))
    cubit.cmd(f"subtract volume {inner} from volume {outer}")
    result = int(cubit.get_last_id("volume"))
    cubit.cmd(f"move volume {result} z {center_z:.12g}")
    return result


def build_mesh(cubit, *, output: Path, step_output: Path | None,
               mesh_size: float,
               aperture_radius: float, pole_radius: float,
               return_outer_radius: float, half_gap: float,
               pole_thickness: float,
               axial_pre_scale: float) -> dict[str, object]:
    dimensions = [mesh_size, aperture_radius, pole_radius,
                  return_outer_radius, half_gap, pole_thickness,
                  axial_pre_scale]
    if (any(not 0.0 < value for value in dimensions)
            or not aperture_radius < pole_radius < return_outer_radius):
        raise ValueError("invalid positive C-yoke dimensions")
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    step_output = (output.with_suffix(".step") if step_output is None
                   else step_output.resolve())
    step_output.parent.mkdir(parents=True, exist_ok=True)
    scaled_output = output.with_name(output.stem + "_scaled.vol")
    cubit.cmd("reset")

    pole_center = half_gap + 0.5 * pole_thickness
    top = _annulus(
        cubit, inner_radius=aperture_radius, outer_radius=pole_radius,
        height=pole_thickness, center_z=pole_center)
    bottom = _annulus(
        cubit, inner_radius=aperture_radius, outer_radius=pole_radius,
        height=pole_thickness, center_z=-pole_center)
    return_height = 2.0 * (half_gap + pole_thickness)
    return_leg = _annulus(
        cubit, inner_radius=pole_radius,
        outer_radius=return_outer_radius, height=return_height,
        center_z=0.0)
    cubit.cmd(f"unite volume {top} {bottom} {return_leg}")
    source_volumes = [int(cubit.get_last_id("volume"))]
    for command in (
            f'export step "{step_output.as_posix()}" volume all overwrite',
            f"volume all scale z {axial_pre_scale:.12g}",
            (f"sculpt volume all processors 1 size {mesh_size:.12g} "
             "gen_sidesets 2 stair 1"),
            "block 1 add hex all",
            'block 1 name "yoke"',
            # Do not rely only on Sculpt's generated free sideset: the
            # topology-optimization charge basis needs every exterior face.
            # Materializing the single material block skin also makes this
            # contract explicit to the Netgen exporter.
            "skin block 1 make sideset 20",
            'sideset 20 name "yoke_boundary"',
            f'export netgen "{scaled_output.as_posix()}" order 1 overwrite'):
        cubit.cmd(command)
    hexes = int(cubit.get_hex_count())
    sidesets = []
    for sideset_id in cubit.get_sideset_id_list():
        sidesets.append({
            "id": int(sideset_id),
            "name": str(cubit.get_exodus_entity_name(
                "sideset", sideset_id) or ""),
            "quads": len(tuple(cubit.get_sideset_quads(sideset_id))),
            "surfaces": len(tuple(cubit.get_sideset_surfaces(sideset_id))),
        })
    if (hexes <= 0 or not scaled_output.is_file()
            or scaled_output.stat().st_size == 0):
        raise RuntimeError("Cubit did not create a nonempty HEX Netgen mesh")
    return {
        "schema": "radia.ffag-axisymmetric-c-yoke-hex/v1",
        "output": str(output),
        "step_output": str(step_output),
        "mesh_size_m": float(mesh_size),
        "axial_mesh_size_m": float(mesh_size / axial_pre_scale),
        "axial_pre_scale": float(axial_pre_scale),
        "aperture_radius_m": float(aperture_radius),
        "pole_radius_m": float(pole_radius),
        "return_outer_radius_m": float(return_outer_radius),
        "half_gap_m": float(half_gap),
        "pole_thickness_m": float(pole_thickness),
        "source_volume_ids": source_volumes,
        "sculpt_mode": "stair-full",
        "hexes": hexes,
        "nodes": int(cubit.get_node_count()),
        "sidesets": sidesets,
        "bytes": int(scaled_output.stat().st_size),
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--step-output", type=Path)
    parser.add_argument("--mesh-size", type=float, default=2.50)
    parser.add_argument("--aperture-radius", type=float, default=1.50)
    parser.add_argument("--pole-radius", type=float, default=7.50)
    parser.add_argument("--return-outer-radius", type=float, default=10.0)
    parser.add_argument("--half-gap", type=float, default=0.12)
    parser.add_argument("--pole-thickness", type=float, default=0.38)
    parser.add_argument("--axial-pre-scale", type=float, default=75.0)
    parser.add_argument("--cubit-bin", type=Path)
    parser.add_argument("--cubit-child", action="store_true",
                        help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def main() -> int:
    args = parse_args()
    if args.cubit_child:
        cubit = _load_cubit(_find_cubit_bin(args.cubit_bin))
        report = build_mesh(
            cubit, output=args.output, step_output=args.step_output,
            mesh_size=args.mesh_size,
            aperture_radius=args.aperture_radius,
            pole_radius=args.pole_radius,
            return_outer_radius=args.return_outer_radius,
            half_gap=args.half_gap, pole_thickness=args.pole_thickness,
            axial_pre_scale=args.axial_pre_scale)
        print("__RADIA_CUBIT_REPORT__" + json.dumps(report))
        return 0

    output = args.output.resolve()
    scaled_output = output.with_name(output.stem + "_scaled.vol")
    stale_sidecars = (
        output.with_suffix(output.suffix + ".json"),
        scaled_output.with_suffix(scaled_output.suffix + ".json"),
    )
    for path in (output, scaled_output, *stale_sidecars):
        path.unlink(missing_ok=True)
    child = [sys.executable, str(Path(__file__).resolve()),
             "--cubit-child", "--output", str(output),
             "--mesh-size", str(args.mesh_size),
             "--aperture-radius", str(args.aperture_radius),
             "--pole-radius", str(args.pole_radius),
             "--return-outer-radius", str(args.return_outer_radius),
             "--half-gap", str(args.half_gap),
             "--pole-thickness", str(args.pole_thickness),
             "--axial-pre-scale", str(args.axial_pre_scale)]
    if args.step_output is not None:
        child.extend(("--step-output", str(args.step_output.resolve())))
    if args.cubit_bin is not None:
        child.extend(("--cubit-bin", str(args.cubit_bin.resolve())))
    completed = subprocess.run(
        child, capture_output=True, text=True, errors="replace")
    child_report = None
    for line in completed.stdout.splitlines():
        if line.startswith("__RADIA_CUBIT_REPORT__"):
            child_report = json.loads(line.removeprefix(
                "__RADIA_CUBIT_REPORT__"))
    if not scaled_output.is_file() or scaled_output.stat().st_size == 0:
        tail = "\n".join((completed.stdout + completed.stderr).splitlines()[-40:])
        raise RuntimeError(
            "Cubit Sculpt/export did not produce the scaled HEX mesh "
            f"(exit={completed.returncode}):\n{tail}")
    from radia.topopt_cad import relabel_straight_mesh, rescale_netgen_vol_points
    rescale_netgen_vol_points(
        scaled_output, output, (1.0, 1.0, 1.0 / args.axial_pre_scale))
    import ngsolve as ng
    with ng.TaskManager():
        raw_mesh = ng.Mesh(str(output))
        # The Cubit-side inventory above proves Sculpt created the complete free
        # skin.  Some installed exporter builds omit free quads from Netgen VOL,
        # so rebuild that same skin exactly from volume-face ownership.  This is a
        # connectivity-preserving relabel, not a geometric approximation.
        mesh = relabel_straight_mesh(
            raw_mesh, lambda _center, _normal: "yoke_boundary",
            material_name="yoke")
        mesh.ngmesh.Save(str(output))
        mesh = ng.Mesh(str(output))
    cubit_skin_faces = 0
    if child_report is not None:
        cubit_skin_faces = max(
            (row["quads"] for row in child_report["sidesets"]), default=0)
    if cubit_skin_faces <= 0 or mesh.GetNE(ng.BND) != cubit_skin_faces:
        raise RuntimeError(
            "Cubit Sculpt skin / reconstructed Netgen boundary mismatch: "
            f"Cubit={cubit_skin_faces}, Netgen={mesh.GetNE(ng.BND)}")
    pole_intervals = []
    for element in mesh.Elements(ng.VOL):
        coordinates = np.asarray([
            tuple(mesh[vertex].point) for vertex in element.vertices],
            dtype=float)
        if np.max(np.linalg.norm(coordinates[:, :2], axis=1)) <= args.pole_radius:
            pole_intervals.append((float(np.min(coordinates[:, 2])),
                                   float(np.max(coordinates[:, 2]))))
    pole_intervals = np.asarray(pole_intervals, dtype=float)
    if (pole_intervals.size == 0 or np.any(
            (pole_intervals[:, 0] < 0.0)
            & (pole_intervals[:, 1] > 0.0))):
        raise RuntimeError(
            "Sculpt HEX cells bridge the beam gap inside the pole radius")
    negative_faces = pole_intervals[pole_intervals[:, 1] < 0.0, 1]
    positive_faces = pole_intervals[pole_intervals[:, 0] > 0.0, 0]
    if negative_faces.size == 0 or positive_faces.size == 0:
        raise RuntimeError("Sculpt mesh does not contain both pole packs")
    lower_gap_face = float(np.max(negative_faces))
    upper_gap_face = float(np.min(positive_faces))
    axial_size = args.mesh_size / args.axial_pre_scale
    if max(abs(lower_gap_face + args.half_gap),
           abs(upper_gap_face - args.half_gap)) > axial_size + 1.0e-12:
        raise RuntimeError(
            "Sculpt beam gap is under-resolved: requested faces at "
            f"+/-{args.half_gap:.6g} m, obtained {lower_gap_face:.6g} and "
            f"{upper_gap_face:.6g} m")
    from radia.topology_optimization import ngsolve_growth_topology
    with ng.TaskManager():
        topology = ngsolve_growth_topology(
            mesh, np.ones(mesh.ne, dtype=bool))
    if not topology.iron_connected:
        raise RuntimeError(
            "Sculpt C-yoke is disconnected; refine the radial return leg")
    report = {
        "schema": "radia.ffag-axisymmetric-c-yoke-hex/v1",
        "output": str(output),
        "step_output": str(
            (output.with_suffix(".step") if args.step_output is None
             else args.step_output.resolve())),
        "mesh_size_m": args.mesh_size,
        "axial_mesh_size_m": args.mesh_size / args.axial_pre_scale,
        "axial_pre_scale": args.axial_pre_scale,
        "aperture_radius_m": args.aperture_radius,
        "pole_radius_m": args.pole_radius,
        "return_outer_radius_m": args.return_outer_radius,
        "half_gap_m": args.half_gap,
        "pole_thickness_m": args.pole_thickness,
        "sculpt_mode": "stair-full",
        "hexes": int(mesh.ne),
        "nodes": int(mesh.nv),
        "boundary_faces": int(mesh.GetNE(ng.BND)),
        "resolved_lower_gap_face_m": lower_gap_face,
        "resolved_upper_gap_face_m": upper_gap_face,
        "iron_connected": bool(topology.iron_connected),
        "azimuthal_periodicity": {
            "representation": "explicit-full-annulus",
            "sector_cut": False,
            "periodic_hdiv_trace_required": False,
            "image_cyclic_used": False,
        },
        "cubit_sidesets": ([] if child_report is None
                            else child_report["sidesets"]),
        "bytes": int(output.stat().st_size),
        "cubit_child_exit_code": completed.returncode,
    }
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
