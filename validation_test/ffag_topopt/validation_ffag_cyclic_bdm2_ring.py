"""FFAG rotational-periodicity gate for the production BDM2 HEX route.

The fixture contains upper and lower iron pole packs wholly inside one FFAG
cell.  Their rotated copies are disjoint, so the one-cell model may use
``image_cyclic`` without an FE seam identification.  The reduced result is
compared with an explicitly meshed full ring made from exact rotations of the
same HEX cells.  Both the material magnetization and the gap-field samples are
checked; a lone cell must remain detectably different so the test cannot pass
when the cyclic images are a no-op.

This lane deliberately does not claim that a continuous annular return yoke
can be reduced the same way.  Such a yoke crosses the azimuthal cut planes and
requires rotation-related HDiv normal-trace identification in addition to the
nonlocal image sum.  ``validate_ffag_cyclic_sector_contract`` fails loudly
when that second contract is omitted.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import time
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


def _rotation(angle):
    cosine = math.cos(angle)
    sine = math.sin(angle)
    return np.asarray(((cosine, -sine, 0.0),
                       (sine, cosine, 0.0),
                       (0.0, 0.0, 1.0)))


def _pole_vertices(*, radius, radial_half_width, tangential_half_width,
                   half_gap, pole_thickness, upper):
    z0, z1 = ((half_gap, half_gap + pole_thickness) if upper else
              (-half_gap - pole_thickness, -half_gap))
    x0, x1 = radius - radial_half_width, radius + radial_half_width
    y0, y1 = -tangential_half_width, tangential_half_width
    return np.asarray((
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ), dtype=float)


def build_ffag_pole_ring_mesh(
        sector_indices, *, fold, radius, radial_half_width,
        tangential_half_width, half_gap, pole_thickness):
    """Build exact rotated copies of one two-pole HEX FFAG cell."""
    import ngsolve as ng
    from netgen.csg import Pnt
    from netgen.meshing import (
        Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint)

    netmesh = Mesh(dim=3)
    netmesh.SetMaterial(1, "pole")
    boundary = netmesh.Add(
        FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    netmesh.SetBCName(0, "pole_boundary")
    reference_poles = tuple(
        _pole_vertices(
            radius=radius, radial_half_width=radial_half_width,
            tangential_half_width=tangential_half_width,
            half_gap=half_gap, pole_thickness=pole_thickness,
            upper=upper)
        for upper in (False, True))
    for sector in sector_indices:
        rotation = _rotation(2.0 * math.pi * int(sector) / int(fold))
        for reference in reference_poles:
            coordinates = reference @ rotation.T
            points = [netmesh.Add(MeshPoint(Pnt(*point)))
                      for point in coordinates]
            netmesh.Add(Element3D(1, points))
            for face in _HEX_FACES:
                netmesh.Add(
                    Element2D(boundary, [points[index] for index in face]))
    return ng.Mesh(netmesh)


def _mean_mz(result, mesh):
    import ngsolve as ng

    volume = float(ng.Integrate(ng.CF(1.0), mesh))
    return float(ng.Integrate(result["gfM"][2], mesh)) / volume


def run_validation(args):
    import ngsolve as ng
    import radia as rad
    import radia.vim as vim
    from radia.ffag_topopt import validate_ffag_cyclic_sector_contract

    contract = validate_ffag_cyclic_sector_contract(args.fold)
    common = dict(
        fold=args.fold, radius=args.radius,
        radial_half_width=args.radial_half_width,
        tangential_half_width=args.tangential_half_width,
        half_gap=args.half_gap, pole_thickness=args.pole_thickness)
    started = time.perf_counter()
    ng.SetNumThreads(args.threads)
    with ng.TaskManager():
        ring = build_ffag_pole_ring_mesh(range(args.fold), **common)
        sector = build_ffag_pole_ring_mesh((0,), **common)
        external = ng.CF((0.0, 0.0, args.applied_field_am))

        rad.UtiDelAll()
        full = vim.Solve(
            ring, mu_r=args.mu_r, H_ext=external, order=args.hdiv_order,
            gram_eps=args.gram_eps, tol=args.solve_tol)
        rad.UtiDelAll()
        lone = vim.Solve(
            sector, mu_r=args.mu_r, H_ext=external,
            order=args.hdiv_order, gram_eps=args.gram_eps,
            tol=args.solve_tol)
        rad.UtiDelAll()
        cyclic = vim.Solve(
            sector, mu_r=args.mu_r, H_ext=external,
            order=args.hdiv_order, image_cyclic=args.fold,
            gram_eps=args.gram_eps, tol=args.solve_tol)

        full_mz = _mean_mz(full, ring)
        lone_mz = _mean_mz(lone, sector)
        cyclic_mz = _mean_mz(cyclic, sector)
        radial_offsets = np.linspace(
            -0.75 * args.radial_half_width,
            0.75 * args.radial_half_width, args.probe_count)
        probes = np.column_stack((
            args.radius + radial_offsets,
            np.zeros(args.probe_count), np.zeros(args.probe_count)))
        full_h = np.asarray(
            vim.FieldFromSolution(full, probes, algorithm="direct"), float)
        cyclic_h = np.asarray(
            vim.FieldFromSolution(cyclic, probes, algorithm="direct"), float)

    magnetization_gap = abs(lone_mz - full_mz)
    magnetization_residual = abs(cyclic_mz - full_mz)
    gap_fraction_closed = (
        1.0 - magnetization_residual / magnetization_gap
        if magnetization_gap > 0.0 else -math.inf)
    field_scale = max(float(np.max(np.linalg.norm(full_h, axis=1))), 1.0)
    field_relative_error = float(
        np.max(np.linalg.norm(cyclic_h - full_h, axis=1)) / field_scale)
    gates = {
        "lone_cell_is_discriminating": bool(
            magnetization_gap / max(abs(full_mz), 1.0) > 1.0e-3),
        "cyclic_closes_magnetization_gap": bool(
            gap_fraction_closed >= args.minimum_gap_fraction_closed),
        "cyclic_gap_field_matches_full_ring": bool(
            field_relative_error <= args.maximum_field_relative_error),
        "bdm2_hex_is_used": bool(
            args.hdiv_order == 2 and set(
                len(element.vertices) for element in sector.Elements(ng.VOL)) == {8}),
        "one_cell_reduces_dofs": bool(
            int(cyclic["ndof"]) * (args.fold - 1) < int(full["ndof"])),
    }
    report = {
        "schema": "radia.ffag-cyclic-bdm2-ring/v1",
        "status": "pass" if all(gates.values()) else "fail",
        "scope": (
            "Disjoint two-pole FFAG cells. Continuous return-yoke sector "
            "cuts require a separate rotation-related HDiv trace gate."),
        "contract": {
            "fold": contract.fold,
            "field_antiperiodic": contract.field_antiperiodic,
            "body_crosses_periodic_planes": (
                contract.body_crosses_periodic_planes),
            "periodic_trace_identified": contract.periodic_trace_identified,
            "reduction_mode": contract.reduction_mode,
        },
        "settings": vars(args) | {"output": str(args.output)},
        "mesh": {
            "full_ring_elements": int(ring.ne),
            "sector_elements": int(sector.ne),
            "full_ring_dofs": int(full["ndof"]),
            "sector_dofs": int(cyclic["ndof"]),
        },
        "magnetization_am": {
            "full_ring_mean_mz": full_mz,
            "lone_cell_mean_mz": lone_mz,
            "cyclic_cell_mean_mz": cyclic_mz,
            "lone_full_gap": magnetization_gap,
            "cyclic_full_residual": magnetization_residual,
            "gap_fraction_closed": gap_fraction_closed,
        },
        "gap_field": {
            "probe_points_m": probes.tolist(),
            "full_ring_demag_h_am": full_h.tolist(),
            "cyclic_demag_h_am": cyclic_h.tolist(),
            "maximum_relative_error": field_relative_error,
        },
        "gates": gates,
        "timings": {
            "full_ring_total_wall_s": float(full["total_wall_s_internal"]),
            "cyclic_total_wall_s": float(cyclic["total_wall_s_internal"]),
            "validation_wall_s": time.perf_counter() - started,
        },
        "runtime": {
            "host": platform.node(),
            "python": platform.python_version(),
        },
    }
    return report


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fold", type=int, default=6)
    parser.add_argument("--hdiv-order", type=int, default=2, choices=(1, 2))
    parser.add_argument("--radius", type=float, default=0.040)
    parser.add_argument("--radial-half-width", type=float, default=0.008)
    parser.add_argument("--tangential-half-width", type=float, default=0.006)
    parser.add_argument("--half-gap", type=float, default=0.004)
    parser.add_argument("--pole-thickness", type=float, default=0.006)
    parser.add_argument("--mu-r", type=float, default=500.0)
    parser.add_argument("--applied-field-am", type=float, default=1.0e5)
    parser.add_argument("--gram-eps", type=float, default=1.0e-10)
    parser.add_argument("--solve-tol", type=float, default=1.0e-10)
    parser.add_argument("--probe-count", type=int, default=5)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--minimum-gap-fraction-closed", type=float,
                        default=0.999)
    parser.add_argument("--maximum-field-relative-error", type=float,
                        default=2.0e-8)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if (args.fold < 2 or args.probe_count < 2 or args.threads < 1
            or args.mu_r <= 1.0 or args.gram_eps <= 0.0
            or args.solve_tol <= 0.0):
        raise ValueError("invalid FFAG cyclic validation settings")
    report = run_validation(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
