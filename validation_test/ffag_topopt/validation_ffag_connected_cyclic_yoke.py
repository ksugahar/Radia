"""Connected-yoke FFAG cyclic-sector validation for BDM2 HEX.

The annular yoke crosses both azimuthal sector planes.  A tangential applied
field gives nonzero normal magnetization on those cuts, making an image-only
sector strongly wrong.  The production reduced route therefore combines:

* NGSolve PERIODIC point identifications on the two rotation-related faces;
* ``Compress(Periodic(HDiv))`` for the local normal-trace constraint;
* exclusion of the paired seam faces from the physical charge skin; and
* ``image_cyclic=N`` for the nonlocal full-ring interaction.

The reduced BDM2 result is compared with an explicit conforming full ring.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import time
from pathlib import Path

import numpy as np


_PHYSICAL_FACES = (
    (0, 3, 2, 1),
    (4, 5, 6, 7),
    (1, 2, 6, 5),
    (3, 0, 4, 7),
)
_MIN_FACE = (0, 1, 5, 4)
_MAX_FACE = (2, 3, 7, 6)


def build_connected_yoke_mesh(*, fold, inner_radius, outer_radius,
                              half_height, full_ring):
    import ngsolve as ng
    from netgen.csg import Pnt
    from netgen.meshing import (
        Element2D, Element3D, FaceDescriptor, IdentificationType,
        Mesh, MeshPoint,
    )

    netmesh = Mesh(dim=3)
    netmesh.SetMaterial(1, "yoke")
    descriptors = []
    for index, name in enumerate(
            ("skin", "periodic_min", "periodic_max"), start=1):
        descriptors.append(netmesh.Add(FaceDescriptor(
            surfnr=index, domin=1, domout=0, bc=index)))
        netmesh.SetBCName(index - 1, name)

    if full_ring:
        points = {}
        for sector in range(fold):
            angle = 2.0 * math.pi * sector / fold
            for radial_index, radius in enumerate(
                    (inner_radius, outer_radius)):
                for axial_index, z_value in enumerate(
                        (-half_height, half_height)):
                    points[sector, radial_index, axial_index] = netmesh.Add(
                        MeshPoint(Pnt(
                            radius * math.cos(angle),
                            radius * math.sin(angle), z_value)))
        for sector in range(fold):
            following = (sector + 1) % fold
            element = [
                points[sector, 0, 0], points[sector, 1, 0],
                points[following, 1, 0], points[following, 0, 0],
                points[sector, 0, 1], points[sector, 1, 1],
                points[following, 1, 1], points[following, 0, 1],
            ]
            netmesh.Add(Element3D(1, element))
            for face in _PHYSICAL_FACES:
                netmesh.Add(Element2D(
                    descriptors[0], [element[index] for index in face]))
        return ng.Mesh(netmesh)

    raw_points = []
    for angle in (0.0, 2.0 * math.pi / fold):
        for radius in (inner_radius, outer_radius):
            for z_value in (-half_height, half_height):
                raw_points.append(netmesh.Add(MeshPoint(Pnt(
                    radius * math.cos(angle), radius * math.sin(angle),
                    z_value))))
    element = [raw_points[index] for index in (0, 2, 6, 4, 1, 3, 7, 5)]
    netmesh.Add(Element3D(1, element))
    for face in _PHYSICAL_FACES:
        netmesh.Add(Element2D(
            descriptors[0], [element[index] for index in face]))
    netmesh.Add(Element2D(
        descriptors[1], [element[index] for index in _MIN_FACE]))
    netmesh.Add(Element2D(
        descriptors[2], [element[index] for index in _MAX_FACE]))
    for master, slave in zip(raw_points[:4], raw_points[4:]):
        netmesh.AddPointIdentification(
            master, slave, identnr=1, type=IdentificationType.PERIODIC)
    return ng.Mesh(netmesh)


def _mean_tangential_magnetization(result, mesh):
    import ngsolve as ng

    radius = ng.sqrt(ng.x * ng.x + ng.y * ng.y)
    tangent = ng.CF((-ng.y / radius, ng.x / radius, 0.0))
    volume = float(ng.Integrate(ng.CF(1.0), mesh))
    return float(ng.Integrate(result["gfM"] * tangent, mesh)) / volume


def run_validation(args):
    import ngsolve as ng
    import radia as rad
    import radia.vim as vim
    from radia.ffag_topopt import validate_ffag_cyclic_sector_contract

    contract = validate_ffag_cyclic_sector_contract(
        args.fold, body_crosses_periodic_planes=True,
        periodic_trace_identified=True)
    full = build_connected_yoke_mesh(
        fold=args.fold, inner_radius=args.inner_radius,
        outer_radius=args.outer_radius, half_height=args.half_height,
        full_ring=True)
    sector = build_connected_yoke_mesh(
        fold=args.fold, inner_radius=args.inner_radius,
        outer_radius=args.outer_radius, half_height=args.half_height,
        full_ring=False)
    radius = ng.sqrt(ng.x * ng.x + ng.y * ng.y)
    applied = ng.CF((
        -args.applied_field_am * ng.y / radius,
        args.applied_field_am * ng.x / radius, 0.0))
    started = time.perf_counter()
    ng.SetNumThreads(args.threads)
    with ng.TaskManager():
        rad.UtiDelAll()
        full_result = vim.Solve(
            full, mu_r=args.mu_r, H_ext=applied,
            order=args.hdiv_order, gram_eps=args.gram_eps, leaf=args.leaf,
            tol=args.solve_tol)
        rad.UtiDelAll()
        sector_result = vim.Solve(
            sector, mu_r=args.mu_r, H_ext=applied,
            order=args.hdiv_order, image_cyclic=args.fold,
            cyclic_periodic_boundaries=("periodic_min", "periodic_max"),
            gram_eps=args.gram_eps, leaf=args.leaf, tol=args.solve_tol)
        full_mean = _mean_tangential_magnetization(full_result, full)
        sector_mean = _mean_tangential_magnetization(sector_result, sector)
        probe_radius = np.linspace(
            0.8 * args.inner_radius, 1.2 * args.outer_radius,
            args.probe_count)
        probes = np.column_stack((
            probe_radius, np.zeros(args.probe_count),
            np.full(args.probe_count, 1.5 * args.half_height)))
        full_field = np.asarray(vim.FieldFromSolution(
            full_result, probes, algorithm="direct"), dtype=float)
        sector_field = np.asarray(vim.FieldFromSolution(
            sector_result, probes, algorithm="direct"), dtype=float)

    mean_relative_error = abs(sector_mean - full_mean) / abs(full_mean)
    field_relative_to_drive = float(np.max(np.linalg.norm(
        sector_field - full_field, axis=1)) / args.applied_field_am)
    gates = {
        "connected_periodic_contract": bool(
            contract.reduction_mode == "connected-periodic-fem-sector"),
        "bdm2_hex_is_used": bool(
            args.hdiv_order == 2 and set(
                len(element.vertices) for element in
                sector.Elements(ng.VOL)) == {8}),
        "periodic_trace_dofs_are_compressed": bool(
            int(sector_result["periodic_slave_dofs"]) > 0),
        "sector_dofs_scale_as_one_fold": bool(
            int(sector_result["ndof"]) * args.fold
            == int(full_result["ndof"])),
        "mean_magnetization_matches_full_ring": bool(
            mean_relative_error <= args.maximum_mean_relative_error),
        "external_field_matches_full_ring": bool(
            field_relative_to_drive <= args.maximum_field_relative_to_drive),
    }
    return {
        "schema": "radia.ffag-connected-cyclic-yoke/v1",
        "status": "pass" if all(gates.values()) else "fail",
        "scope": (
            "Connected annular HEX yoke with nonzero normal magnetization "
            "on the azimuthal cuts; periodic trace plus cyclic images."),
        "settings": vars(args) | {"output": str(args.output)},
        "contract": {
            "fold": contract.fold,
            "formulation": contract.formulation,
            "reduction_mode": contract.reduction_mode,
            "body_crosses_periodic_planes": True,
            "periodic_trace_identified": True,
            "periodic_charge_paired": False,
            "periodic_boundaries": ["periodic_min", "periodic_max"],
        },
        "mesh": {
            "full_ring_elements": int(full.ne),
            "sector_elements": int(sector.ne),
            "full_ring_dofs": int(full_result["ndof"]),
            "sector_dofs": int(sector_result["ndof"]),
            "periodic_slave_dofs": int(
                sector_result["periodic_slave_dofs"]),
        },
        "mean_tangential_magnetization_am": {
            "full_ring": full_mean,
            "periodic_sector": sector_mean,
            "relative_error": mean_relative_error,
        },
        "external_field": {
            "probe_points_m": probes.tolist(),
            "full_ring_demag_h_am": full_field.tolist(),
            "periodic_sector_demag_h_am": sector_field.tolist(),
            "maximum_difference_relative_to_applied_field": (
                field_relative_to_drive),
        },
        "timings": {
            "full_ring_total_wall_s": float(
                full_result["total_wall_s_internal"]),
            "periodic_sector_total_wall_s": float(
                sector_result["total_wall_s_internal"]),
            "validation_wall_s": time.perf_counter() - started,
        },
        "gates": gates,
        "runtime": {
            "host": platform.node(),
            "python": platform.python_version(),
        },
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fold", type=int, default=12)
    parser.add_argument("--hdiv-order", type=int, choices=(1, 2), default=2)
    parser.add_argument("--inner-radius", type=float, default=0.030)
    parser.add_argument("--outer-radius", type=float, default=0.050)
    parser.add_argument("--half-height", type=float, default=0.004)
    parser.add_argument("--mu-r", type=float, default=100.0)
    parser.add_argument("--applied-field-am", type=float, default=1.0e5)
    parser.add_argument("--gram-eps", type=float, default=1.0e-10)
    parser.add_argument("--leaf", type=int, default=32)
    parser.add_argument("--solve-tol", type=float, default=1.0e-10)
    parser.add_argument("--probe-count", type=int, default=5)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--maximum-mean-relative-error", type=float,
                        default=1.0e-8)
    parser.add_argument("--maximum-field-relative-to-drive", type=float,
                        default=2.0e-7)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    report = run_validation(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
