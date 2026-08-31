"""Full-ring BDM2 HDiv-VIM F/D closed-orbit validation.

This is a deliberately small, explicit six-cell ring of fixed/given
magnetization sources.  Every F and D pole is represented by upper/lower HEX
blocks, and the BDM2 source is evaluated by the persistent native C++ field
evaluator while the closed-orbit root uses the native tracker.  It validates
the F/D field, the HDiv field evaluator, and one-cell rotational closure as
one contract.

The fixture is a manufactured permanent-magnet field, not a Cubit C-yoke
engineering design.  The Curved-Q2 Cubit periodic-yoke validation remains the
separate geometry/material gate; this test keeps F/D closed-orbit mechanics
small enough to run before that full design is available.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import time
import traceback
from pathlib import Path

import numpy as np

_HEX_FACES = (
    (0, 3, 2, 1),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)


def _pole_vertices(*, angle, z0, z1, radius, radial_half_width, tangential_half_width):
    radial = np.asarray((math.cos(angle), math.sin(angle), 0.0))
    tangent = np.asarray((-math.sin(angle), math.cos(angle), 0.0))
    vertices = []
    for z_value in (z0, z1):
        for radial_offset, tangential_offset in (
            (-radial_half_width, -tangential_half_width),
            (radial_half_width, -tangential_half_width),
            (radial_half_width, tangential_half_width),
            (-radial_half_width, tangential_half_width),
        ):
            vertices.append(
                radial * (radius + radial_offset)
                + tangent * tangential_offset
                + np.asarray((0.0, 0.0, z_value))
            )
    return np.asarray(vertices, dtype=float)


def build_fd_ring_mesh(
    *,
    fold,
    radius,
    radial_half_width,
    tangential_half_width,
    half_gap,
    pole_height,
    f_phase,
    d_phase,
):
    """Build an explicit F/D ring from disconnected upper/lower HEX poles."""
    import ngsolve as ng
    from netgen.csg import Pnt
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint

    netmesh = Mesh(dim=3)
    descriptors = []
    for material, name in enumerate(("f_pole", "d_pole"), start=1):
        netmesh.SetMaterial(material, name)
        descriptor = netmesh.Add(
            FaceDescriptor(surfnr=material, domin=material, domout=0, bc=material)
        )
        netmesh.SetBCName(material - 1, "skin")
        descriptors.append(descriptor)
    cell_angle = 2.0 * math.pi / fold
    for sector in range(fold):
        start_angle = -0.5 * math.pi + sector * cell_angle
        for material, phase in ((1, f_phase), (2, d_phase)):
            angle = start_angle + phase * cell_angle
            for z0, z1 in (
                (half_gap, half_gap + pole_height),
                (-half_gap - pole_height, -half_gap),
            ):
                vertices = _pole_vertices(
                    angle=angle,
                    z0=z0,
                    z1=z1,
                    radius=radius,
                    radial_half_width=radial_half_width,
                    tangential_half_width=tangential_half_width,
                )
                points = [netmesh.Add(MeshPoint(Pnt(*vertex))) for vertex in vertices]
                netmesh.Add(Element3D(material, points))
                for face in _HEX_FACES:
                    netmesh.Add(
                        Element2D(
                            descriptors[material - 1], [points[index] for index in face]
                        )
                    )
    return ng.Mesh(netmesh)


class _SourceField:
    """Public B adapter for the C++ BDM2 source evaluator."""

    def __init__(self, source, mu0):
        self._source = source
        self._mu0 = float(mu0)

    def b_field(self, points):
        values = np.asarray(points, dtype=float)
        single = values.shape == (3,)
        packed = np.ascontiguousarray(values.reshape(-1, 3))
        field = self._mu0 * self._source.Field(packed, algorithm="direct")
        return field[0] if single else field


def _cell_phase(positions, cell_angle):
    angles = np.arctan2(positions[:, 1], positions[:, 0])
    return np.mod(angles + 0.5 * math.pi, cell_angle) / cell_angle


def _band_mean(values, phase, start, stop):
    selected = values[(phase >= start) & (phase <= stop)]
    if selected.size < 2:
        raise RuntimeError("F/D diagnostic band has too few orbit stations")
    return float(np.mean(selected))


def run_validation(args):
    import ngsolve as ng

    from radia import vim
    from radia.ffag_topopt import recover_periodic_planar_closed_orbit_native
    from radia.isochronous_topopt import MU0

    print("[1/6] Building explicit full-ring BDM2 F/D source mesh", flush=True)
    cell_angle = 2.0 * math.pi / args.fold
    mesh = build_fd_ring_mesh(
        fold=args.fold,
        radius=args.radius,
        radial_half_width=args.radial_half_width,
        tangential_half_width=args.tangential_half_width,
        half_gap=args.half_gap,
        pole_height=args.pole_height,
        f_phase=args.f_phase,
        d_phase=args.d_phase,
    )
    magnetization = mesh.MaterialCF(
        {
            "f_pole": ng.CF((0.0, 0.0, args.f_magnetization_am)),
            "d_pole": ng.CF((0.0, 0.0, args.d_magnetization_am)),
        }
    )
    started = time.perf_counter()
    ng.SetNumThreads(args.threads)
    with ng.TaskManager():
        source = vim.MagnetizationSource(
            mesh, magnetization, order=2, curve_gauss=args.curve_gauss
        )
    field = _SourceField(source, MU0)

    print("[2/6] Sampling nominal bend field", flush=True)

    nominal_phase = np.linspace(0.0, 1.0, args.nominal_probe_count)
    nominal_angles = -0.5 * math.pi + nominal_phase * cell_angle
    nominal_points = np.column_stack(
        (
            args.radius * np.cos(nominal_angles),
            args.radius * np.sin(nominal_angles),
            np.zeros(args.nominal_probe_count),
        )
    )
    nominal_bz = np.asarray(field.b_field(nominal_points), dtype=float)[:, 2]
    nominal_bz_mean = float(np.mean(nominal_bz))
    if nominal_bz_mean <= args.minimum_nominal_bz_t:
        raise RuntimeError(
            "manufactured F/D ring has no positive nominal bend field: "
            f"{nominal_bz_mean:.6e} T"
        )
    recovery_common = {
        "iron_evaluator": source._field_evaluator,
        "iron_scale": MU0 / (4.0 * math.pi),
        "magnetic_rigidity": nominal_bz_mean * args.radius,
        "cell_angle_rad": cell_angle,
        "n_segments": args.orbit_segments,
        "gradient_offset": args.gradient_offset,
        "iron_algorithm": args.iron_algorithm,
    }
    print("[3/6] Recovering a coarse native closed-orbit seed", flush=True)
    coarse_orbit = recover_periodic_planar_closed_orbit_native(
        field,
        initial_radius_m=args.radius,
        initial_incidence_angle_rad=args.initial_incidence_angle_rad,
        tracking_step_m=args.coarse_tracking_step_m,
        position_tolerance=args.coarse_position_tolerance,
        tangent_tolerance=args.coarse_tangent_tolerance,
        root_max_evaluations=args.root_max_evaluations,
        **recovery_common,
    )
    print("[4/6] Reclosing the seed with the accepted direct integration", flush=True)
    orbit = recover_periodic_planar_closed_orbit_native(
        field,
        initial_radius_m=coarse_orbit.entrance_radius_m,
        initial_incidence_angle_rad=coarse_orbit.entrance_incidence_angle_rad,
        tracking_step_m=args.tracking_step_m,
        position_tolerance=args.position_tolerance,
        tangent_tolerance=args.tangent_tolerance,
        root_max_evaluations=args.fine_root_max_evaluations,
        **recovery_common,
    )

    print("[5/6] Classifying F/D gradients and field-tree accuracy", flush=True)
    station_count = args.orbit_segments
    response = np.asarray(orbit.field_response, dtype=float)
    if response.shape != (2 * station_count,):
        raise RuntimeError("native closed-orbit field response has an unexpected shape")
    phase = _cell_phase(orbit.orbit.positions[:-1], cell_angle)
    bend = response[:station_count]
    gradient = response[station_count:]
    f_bend = _band_mean(bend, phase, args.f_band_start, args.f_band_stop)
    d_bend = _band_mean(bend, phase, args.d_band_start, args.d_band_stop)
    f_gradient = _band_mean(gradient, phase, args.f_band_start, args.f_band_stop)
    d_gradient = _band_mean(gradient, phase, args.d_band_start, args.d_band_stop)
    orbit_points = np.ascontiguousarray(orbit.orbit.positions[:-1], dtype=float)
    direct_field = MU0 * source.Field(orbit_points, algorithm="direct")
    tree_field = MU0 * source.Field(orbit_points, algorithm="tree")
    field_scale = max(
        float(np.max(np.linalg.norm(direct_field, axis=1))), np.finfo(float).tiny
    )
    tree_direct_field_relative_error = float(
        np.max(np.linalg.norm(tree_field - direct_field, axis=1)) / field_scale
    )
    gradient_signs_oppose = f_gradient * d_gradient < 0.0
    gradient_strength = min(abs(f_gradient), abs(d_gradient))
    gates = {
        "explicit_full_ring": bool(mesh.ne == 4 * args.fold),
        "bdm2_source_is_used": bool(source.order == 2),
        "source_projection_converged": bool(
            source.stats["projection_relative_residual"]
            <= args.maximum_projection_residual
        ),
        "nominal_bend_field_is_positive": bool(
            nominal_bz_mean > args.minimum_nominal_bz_t
        ),
        "native_closed_orbit_position": bool(
            orbit.periodic_position_residual_m <= args.position_tolerance
        ),
        "native_closed_orbit_tangent": bool(
            orbit.periodic_tangent_residual <= args.tangent_tolerance
        ),
        "tree_field_matches_direct_field": bool(
            tree_direct_field_relative_error
            <= args.maximum_tree_direct_field_relative_error
        ),
        "f_and_d_gradients_oppose": bool(gradient_signs_oppose),
        "f_and_d_gradients_are_resolved": bool(
            gradient_strength >= args.minimum_gradient_t_per_m
        ),
    }
    report = {
        "schema": "radia.ffag-bdm2-fd-closed-orbit/v1",
        "status": "pass" if all(gates.values()) else "fail",
        "scope": (
            "Explicit full-ring fixed/given permanent-magnet BDM2 HDiv-VIM "
            "F/D field with native one-cell closed-orbit recovery. This is a "
            "manufactured source/tracker gate, not the Cubit C-yoke design."
        ),
        "model": {
            "fold": int(args.fold),
            "hdiv_family": "BDM2",
            "geometry": "affine disconnected HEX permanent-magnet poles",
            "permanent_magnet_model": "fixed-given",
            "native_field_evaluator": "C++ direct element cloud",
            "native_orbit_tracker": True,
            "native_tracker_field_algorithm": args.iron_algorithm,
            "hysteresis": False,
        },
        "mesh": {
            "elements": int(mesh.ne),
            "vertices": int(mesh.nv),
            "source_dofs": int(source.space.ndof),
            "source_projection_relative_residual": float(
                source.stats["projection_relative_residual"]
            ),
        },
        "nominal_field": {
            "probe_points_m": nominal_points.tolist(),
            "bz_t": nominal_bz.tolist(),
            "mean_bz_t": nominal_bz_mean,
            "magnetic_rigidity_tm": float(orbit.magnetic_rigidity_tm),
        },
        "closed_orbit": {
            "entrance_radius_m": float(orbit.entrance_radius_m),
            "entrance_incidence_angle_rad": float(orbit.entrance_incidence_angle_rad),
            "position_residual_m": float(orbit.periodic_position_residual_m),
            "tangent_residual": float(orbit.periodic_tangent_residual),
            "vertical_position_residual_m": float(orbit.vertical_position_residual_m),
            "vertical_tangent_residual": float(orbit.vertical_tangent_residual),
            "native_root_track_evaluations": int(orbit.root_evaluations),
            "path_length_m": float(orbit.path_length_m),
        },
        "fd_field_diagnostics": {
            "f_band_phase": [args.f_band_start, args.f_band_stop],
            "d_band_phase": [args.d_band_start, args.d_band_stop],
            "f_mean_bend_t": f_bend,
            "d_mean_bend_t": d_bend,
            "f_mean_gradient_t_per_m": f_gradient,
            "d_mean_gradient_t_per_m": d_gradient,
            "minimum_abs_gradient_t_per_m": gradient_strength,
            "tree_direct_field_relative_error": tree_direct_field_relative_error,
        },
        "settings": vars(args) | {"output": str(args.output)},
        "runtime": {
            "host": platform.node(),
            "python": platform.python_version(),
            "ngsolve": ng.__version__,
            "wall_s": time.perf_counter() - started,
        },
        "gates": gates,
    }
    report["closed_orbit"]["coarse_seed"] = {
        "entrance_radius_m": float(coarse_orbit.entrance_radius_m),
        "entrance_incidence_angle_rad": float(
            coarse_orbit.entrance_incidence_angle_rad
        ),
        "position_residual_m": float(coarse_orbit.periodic_position_residual_m),
        "tangent_residual": float(coarse_orbit.periodic_tangent_residual),
        "native_root_track_evaluations": int(coarse_orbit.root_evaluations),
    }
    print("[6/6] Finalizing acceptance report", flush=True)
    return report


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fold", type=int, default=6)
    parser.add_argument("--radius", type=float, default=1.0)
    parser.add_argument("--radial-half-width", type=float, default=0.16)
    parser.add_argument("--tangential-half-width", type=float, default=0.10)
    parser.add_argument("--half-gap", type=float, default=0.06)
    parser.add_argument("--pole-height", type=float, default=0.10)
    parser.add_argument("--f-phase", type=float, default=0.30)
    parser.add_argument("--d-phase", type=float, default=0.73)
    parser.add_argument("--f-magnetization-am", type=float, default=1.4e6)
    parser.add_argument("--d-magnetization-am", type=float, default=0.7e6)
    parser.add_argument("--curve-gauss", type=int, default=8)
    parser.add_argument("--nominal-probe-count", type=int, default=9)
    parser.add_argument("--orbit-segments", type=int, default=48)
    parser.add_argument("--gradient-offset", type=float, default=0.01)
    parser.add_argument("--tracking-step-m", type=float, default=1.0e-3)
    parser.add_argument("--coarse-tracking-step-m", type=float, default=5.0e-3)
    parser.add_argument(
        "--iron-algorithm", choices=("direct", "tree"), default="direct"
    )
    parser.add_argument("--initial-incidence-angle-rad", type=float, default=0.0)
    parser.add_argument("--position-tolerance", type=float, default=3.0e-7)
    parser.add_argument("--tangent-tolerance", type=float, default=3.0e-7)
    parser.add_argument("--root-max-evaluations", type=int, default=24)
    parser.add_argument("--fine-root-max-evaluations", type=int, default=8)
    parser.add_argument("--minimum-nominal-bz-t", type=float, default=0.02)
    parser.add_argument("--minimum-gradient-t-per-m", type=float, default=0.10)
    parser.add_argument("--maximum-projection-residual", type=float, default=1.0e-12)
    parser.add_argument("--coarse-position-tolerance", type=float, default=2.0e-6)
    parser.add_argument("--coarse-tangent-tolerance", type=float, default=2.0e-6)
    parser.add_argument(
        "--maximum-tree-direct-field-relative-error", type=float, default=2.0e-5
    )
    parser.add_argument("--f-band-start", type=float, default=0.15)
    parser.add_argument("--f-band-stop", type=float, default=0.45)
    parser.add_argument("--d-band-start", type=float, default=0.58)
    parser.add_argument("--d-band-stop", type=float, default=0.88)
    parser.add_argument("--threads", type=int, default=8)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if (
        args.fold < 2
        or args.radius <= 0.0
        or args.radial_half_width <= 0.0
        or args.tangential_half_width <= 0.0
        or args.half_gap <= 0.0
        or args.pole_height <= 0.0
        or args.curve_gauss < 1
        or args.nominal_probe_count < 3
        or args.orbit_segments < 16
        or args.gradient_offset <= 0.0
        or args.tracking_step_m <= 0.0
        or args.position_tolerance <= 0.0
        or args.tangent_tolerance <= 0.0
        or args.coarse_position_tolerance <= 0.0
        or args.coarse_tangent_tolerance <= 0.0
        or args.coarse_tracking_step_m <= args.tracking_step_m
        or args.root_max_evaluations < 4
        or args.fine_root_max_evaluations < 4
        or args.threads < 1
        or args.minimum_nominal_bz_t <= 0.0
        or args.minimum_gradient_t_per_m <= 0.0
        or args.maximum_projection_residual <= 0.0
        or args.maximum_tree_direct_field_relative_error <= 0.0
        or not (0.0 <= args.f_band_start < args.f_band_stop <= 1.0)
        or not (0.0 <= args.d_band_start < args.d_band_stop <= 1.0)
    ):
        parser = argparse.ArgumentParser()
        parser.error("invalid BDM2 F/D closed-orbit validation settings")
    started = time.perf_counter()
    try:
        report = run_validation(args)
    except (RuntimeError, ValueError) as error:
        report = {
            "schema": "radia.ffag-bdm2-fd-closed-orbit/v1",
            "status": "error",
            "error": {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc(),
            },
            "settings": vars(args) | {"output": str(args.output)},
            "runtime": {
                "host": platform.node(),
                "python": platform.python_version(),
                "wall_s": time.perf_counter() - started,
            },
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
