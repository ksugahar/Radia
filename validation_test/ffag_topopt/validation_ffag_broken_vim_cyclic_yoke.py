"""Curved-Q2 broken-HDiv VIM cyclic-face validation for an FFAG annular yoke.

The conforming FEM reduction identifies HDiv trace unknowns.  This validation
targets the distinct material-topology route: every volume element remains
broken, while the two identified sector skins contribute one pulled-back
surface-charge jump.  Rotational images then reproduce that quotient charge
over the full ring.

The direct oracle is an explicit full-ring broken-HDiv ChargeGram evaluated on
the same prescribed periodic magnetizations.  It avoids comparing an
unconstrained full-ring material solve, whose additional non-periodic broken
unknowns are not the algebraic symmetry restriction being validated here.
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
    (0, 3, 2, 1), (4, 5, 6, 7), (1, 2, 6, 5), (3, 0, 4, 7),
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
            for radial_index, radius in enumerate((inner_radius, outer_radius)):
                for axial_index, z_value in enumerate((-half_height, half_height)):
                    points[sector, radial_index, axial_index] = netmesh.Add(
                        MeshPoint(Pnt(radius * math.cos(angle),
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


def _curve_mesh(mesh, order):
    """Curve a programmatic 3D mesh under the NGSolve edge-descriptor contract."""
    from netgen.meshing import EdgeDescriptor

    ngmesh = mesh.ngmesh
    if not ngmesh.EdgeDescriptors():
        for index, face in enumerate(ngmesh.FaceDescriptors(), start=1):
            descriptor = EdgeDescriptor()
            descriptor.edgenr = index
            descriptor.surfnr = (face.surfnr, -1)
            descriptor.domin = face.domin
            descriptor.domout = face.domout
            descriptor.name = face.bcname or f"boundary_{index}"
            assert ngmesh.Add(descriptor) == index
    mesh.Curve(order)


def apply_cyclic_q2_geometry(mesh, args):
    """Install a rotation-equivariant, genuinely curved Q2 geometry map."""
    import ngsolve as ng

    _curve_mesh(mesh, args.geometry_order)
    space = ng.VectorH1(mesh, order=args.geometry_order)
    deformation = ng.GridFunction(space)
    radius_squared = ng.x * ng.x + ng.y * ng.y
    deformation.Interpolate(ng.CF((
        args.radial_deformation_scale * ng.x * radius_squared,
        args.radial_deformation_scale * ng.y * radius_squared,
        args.axial_deformation_scale * ng.z * radius_squared,
    )))
    mesh.SetDeformation(deformation)
    return deformation


def _operator(mesh, args, *, cyclic):
    import ngsolve as ng
    from radia.vim import DemagOperator

    fes = ng.HDiv(
        mesh, order=args.hdiv_order, discontinuous=True)
    options = dict(
        eps=args.gram_eps, leafsize=args.leaf,
        internal_interfaces=True)
    if cyclic:
        options.update(
            image_masks=(0,) * (args.fold - 1),
            image_signs=(1.0,) * (args.fold - 1),
            image_rot_angle=tuple(
                2.0 * math.pi * index / args.fold
                for index in range(1, args.fold)),
            cyclic_periodic_boundaries=("periodic_min", "periodic_max"))
    return fes, DemagOperator(fes, **options)


def _field_case(operator, coefficient, points):
    import ngsolve as ng
    from radia.isochronous_topopt import demag_field_from_solution

    field = ng.GridFunction(operator.space)
    field.Set(coefficient)
    values = np.asarray(demag_field_from_solution(
        operator, field, points, algorithm="direct"), dtype=float)
    charge = np.asarray(
        operator._B @ field.vec.FV().NumPy(), dtype=float).reshape(-1)
    return field, values, charge


def run_validation(args):
    import ngsolve as ng
    from radia.ffag_topopt import validate_ffag_cyclic_sector_contract
    from radia.vim._vim import _broken_hex_face_charge_basis

    contract = validate_ffag_cyclic_sector_contract(
        args.fold, body_crosses_periodic_planes=True,
        formulation="vim-broken", periodic_charge_paired=True)
    full_mesh = build_connected_yoke_mesh(
        fold=args.fold, inner_radius=args.inner_radius,
        outer_radius=args.outer_radius, half_height=args.half_height,
        full_ring=True)
    sector_mesh = build_connected_yoke_mesh(
        fold=args.fold, inner_radius=args.inner_radius,
        outer_radius=args.outer_radius, half_height=args.half_height,
        full_ring=False)
    full_deformation = apply_cyclic_q2_geometry(full_mesh, args)
    sector_deformation = apply_cyclic_q2_geometry(sector_mesh, args)
    ng.SetNumThreads(args.threads)
    started = time.perf_counter()
    with ng.TaskManager():
        full_space, full_operator = _operator(full_mesh, args, cyclic=False)
        sector_space, sector_operator = _operator(
            sector_mesh, args, cyclic=True)
        radius = ng.sqrt(ng.x * ng.x + ng.y * ng.y)
        cases = {
            "axial": ng.CF((0.0, 0.0, 1.0)),
            "tangential": ng.CF((-ng.y / radius, ng.x / radius, 0.0)),
        }
        probe_radius = np.linspace(
            0.8 * args.inner_radius, 1.2 * args.outer_radius,
            args.probe_count)
        probes = np.column_stack((
            probe_radius, np.zeros(args.probe_count),
            np.full(args.probe_count, 1.5 * args.half_height)))
        records = {}
        maximum_field_error = 0.0
        for name, coefficient in cases.items():
            _, full_field, _ = _field_case(
                full_operator, coefficient, probes)
            sector_gf, sector_field, _ = _field_case(
                sector_operator, coefficient, probes)
            scale = max(float(np.max(np.linalg.norm(
                full_field, axis=1))), 1.0e-300)
            relative = float(np.max(np.linalg.norm(
                sector_field - full_field, axis=1)) / scale)
            maximum_field_error = max(maximum_field_error, relative)
            records[name] = {
                "full_ring_demag_h_am_per_unit_m": full_field.tolist(),
                "periodic_sector_demag_h_am_per_unit_m": sector_field.tolist(),
                "maximum_relative_error": relative,
            }
            if name == "tangential":
                paired = _broken_hex_face_charge_basis(
                    sector_space, args.hdiv_order,
                    cyclic_periodic_boundaries=(
                        "periodic_min", "periodic_max"),
                    image_rot_angle=(2.0 * math.pi / args.fold,))
                charge = np.asarray(
                    paired["B"] @ sector_gf.vec.FV().NumPy(),
                    dtype=float).reshape(-1)
                block_size = (args.hdiv_order + 1) ** 2
                block = paired["facet_numbers"].index(
                    paired["periodic_master_facets"][0])
                seam = charge[
                    block * block_size:(block + 1) * block_size]
                seam_charge_max = float(np.max(np.abs(seam)))
                periodic_geometry_relative_residual = float(
                    paired["periodic_geometry_relative_residual"])
                face_nodes = np.asarray(
                    paired["face_nodes"], dtype=float)
                bilinear_centres = 0.25 * (
                    face_nodes[:, 0] + face_nodes[:, 2]
                    + face_nodes[:, 6] + face_nodes[:, 8])
                maximum_q2_face_curvature = float(np.max(np.linalg.norm(
                    face_nodes[:, 4] - bilinear_centres, axis=1)))

    expected_sector_charges = int(full_operator._B.shape[0]) // args.fold
    gates = {
        "broken_hdiv_spaces_are_used": bool(
            full_space.ndof > 0 and sector_space.ndof > 0),
        "cyclic_charge_count_matches_full_ring": bool(
            int(sector_operator._B.shape[0]) == expected_sector_charges),
        "periodic_surface_charge_jump_cancels": bool(
            seam_charge_max <= args.maximum_seam_charge),
        "curved_q2_geometry_is_used": bool(
            full_mesh.GetCurveOrder() == args.geometry_order == 2
            and sector_mesh.GetCurveOrder() == args.geometry_order
            and maximum_q2_face_curvature
            >= args.minimum_q2_face_curvature),
        "periodic_q2_geometry_matches_rotation": bool(
            periodic_geometry_relative_residual
            <= args.maximum_periodic_geometry_relative_error),
        "prescribed_periodic_fields_match_full_ring": bool(
            maximum_field_error <= args.maximum_field_relative_error),
    }
    return {
        "schema": "radia.ffag-broken-vim-cyclic-yoke/v2",
        "status": "pass" if all(gates.values()) else "fail",
        "scope": (
            "BDM2 broken-HDiv VIM periodic charge pairing on a connected "
            "curved-Q2 annular HEX yoke."),
        "settings": vars(args) | {"output": str(args.output)},
        "contract": {
            "fold": contract.fold,
            "formulation": contract.formulation,
            "reduction_mode": contract.reduction_mode,
            "periodic_trace_identified": False,
            "periodic_charge_paired": True,
            "periodic_boundaries": ["periodic_min", "periodic_max"],
        },
        "mesh": {
            "full_ring_elements": int(full_mesh.ne),
            "sector_elements": int(sector_mesh.ne),
            "full_ring_dofs": int(full_space.ndof),
            "sector_dofs": int(sector_space.ndof),
            "full_ring_charges": int(full_operator._B.shape[0]),
            "sector_charges": int(sector_operator._B.shape[0]),
        },
        "probe_points_m": probes.tolist(),
        "prescribed_magnetization": records,
        "periodic_seam_charge_max": seam_charge_max,
        "geometry": {
            "order": int(args.geometry_order),
            "maximum_q2_face_curvature_m": maximum_q2_face_curvature,
            "periodic_q2_node_relative_residual": (
                periodic_geometry_relative_residual),
            "full_deformation_dofs": int(full_deformation.space.ndof),
            "sector_deformation_dofs": int(sector_deformation.space.ndof),
        },
        "gates": gates,
        "timings": {"validation_wall_s": time.perf_counter() - started},
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
    parser.add_argument("--geometry-order", type=int, choices=(2,), default=2)
    parser.add_argument("--inner-radius", type=float, default=0.030)
    parser.add_argument("--outer-radius", type=float, default=0.050)
    parser.add_argument("--half-height", type=float, default=0.004)
    parser.add_argument("--gram-eps", type=float, default=1.0e-10)
    parser.add_argument("--leaf", type=int, default=32)
    parser.add_argument("--probe-count", type=int, default=5)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--radial-deformation-scale", type=float, default=2.5)
    parser.add_argument(
        "--axial-deformation-scale", type=float, default=1.5)
    parser.add_argument(
        "--minimum-q2-face-curvature", type=float, default=1.0e-6)
    parser.add_argument(
        "--maximum-periodic-geometry-relative-error",
        type=float, default=1.0e-12)
    parser.add_argument("--maximum-seam-charge", type=float, default=1.0e-14)
    parser.add_argument(
        "--maximum-field-relative-error", type=float, default=1.0e-11)
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
