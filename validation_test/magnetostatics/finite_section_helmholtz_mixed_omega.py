"""Heavy finite-section Helmholtz validation; not part of the CI test suite."""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

MU0 = 4.0e-7 * math.pi
R_INNER, R_OUTER = 10.0, 13.0
Z_INNER, Z_OUTER = 3.0, 6.0
CURRENT_DENSITY = 1.0e4
CORE_RADIUS, KELVIN_RADIUS = 4.0, 20.0
KELVIN_OFFSET = (60.0, 0.0, 0.0)
AXIS_Z = (0.0, 1.0, 2.0, 3.0)


def exact_bz(z_m: float, order: int = 96) -> float:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    radii = 0.5 * (R_INNER + R_OUTER) + 0.5 * (R_OUTER - R_INNER) * nodes
    axial = 0.5 * (Z_INNER + Z_OUTER) + 0.5 * (Z_OUTER - Z_INNER) * nodes
    radius, source_z = np.meshgrid(radii, axial, indexing="ij")
    weight = np.outer(weights, weights)
    kernel = radius**2 / (radius**2 + (z_m - source_z) ** 2) ** 1.5
    kernel += radius**2 / (radius**2 + (z_m + source_z) ** 2) ** 1.5
    return float(0.5 * MU0 * CURRENT_DENSITY * 0.25 * (R_OUTER - R_INNER) * (Z_OUTER - Z_INNER) * np.sum(weight * kernel))


def build_coil():
    import radia as rad
    from radia.coil_builder import CoilBuilder

    width, height = R_OUTER - R_INNER, Z_OUTER - Z_INNER
    radius = 0.5 * (R_INNER + R_OUTER)
    current = CURRENT_DENSITY * width * height
    objects = []
    for z_center in (-4.5, 4.5):
        builder = CoilBuilder(current).set_start((radius, 0.0, z_center)).set_cross_section(width, height)
        for _ in range(4):
            builder.add_arc(radius, 90.0)
        if not builder.is_closed:
            raise RuntimeError("CoilBuilder path did not close")
        objects.extend(builder.to_radia(arc_max_segment_length=0.5))
    return rad.ObjCnt(objects)


def build_mesh(maxh: float):
    import ngsolve as ng
    from netgen.meshing import IdentificationType
    from netgen.occ import Glue, OCCGeometry, Pnt, Sphere

    core = Sphere(Pnt(0, 0, 0), CORE_RADIUS)
    core.mat("core")
    core.maxh = 0.5 * maxh
    for face in core.faces:
        face.name = "core_air_interface"
    shell = Sphere(Pnt(0, 0, 0), KELVIN_RADIUS)
    shell.mat("air")
    for face in shell.faces:
        face.name = "kelvin_int"
    air = shell - core
    air.mat("air")
    exterior = Sphere(Pnt(*KELVIN_OFFSET), KELVIN_RADIUS)
    exterior.mat("kelvin")
    for face in exterior.faces:
        face.name = "kelvin_ext"
    for vertex in exterior.vertices:
        if vertex.p[2] > 0:
            vertex.name = "GND"
    geometry = Glue([air, core, exterior])
    named = {face.name: face for solid in geometry.solids for face in solid.faces if face.name}
    named["kelvin_int"].Identify(named["kelvin_ext"], "periodic", IdentificationType.PERIODIC)
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(geometry).GenerateMesh(maxh=maxh, grading=0.5))


def solve(maxh: float, order: int) -> dict:
    import ngsolve as ng
    import radia as rad
    from radia.kelvin_material import MU_0, make_kelvin_mu_cf
    from radia.kelvin_solver import project_source_interface_potential, solve_magnetostatic_mixed_total_reduced_omega_kelvin

    started = time.perf_counter()
    coil = build_coil()
    mesh = build_mesh(maxh)
    source = rad.RadiaField(coil, "h")
    with ng.TaskManager():
        core = project_source_interface_potential(mesh, source, "core_air_interface", order=order, relative_tolerance=0.01)
        outer = project_source_interface_potential(mesh, source, "kelvin_int", order=order, relative_tolerance=0.01)
        kelvin_mu = make_kelvin_mu_cf(mesh, KELVIN_RADIUS, KELVIN_OFFSET, kelvin_mats=("kelvin",), mu_r_by_material={})
        mu = mesh.MaterialCF({"core": MU_0, "air": MU_0, "kelvin": kelvin_mu}, default=MU_0)
        result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source, core["potential"], KELVIN_RADIUS, KELVIN_OFFSET,
            mu_cf=mu, reduced_materials=("air",), total_materials=("core", "kelvin"),
            interface_boundary="core_air_interface", order=order,
            kelvin_interface_boundary="kelvin_int", kelvin_source_potential=outer["potential"],
            total_source_h=None, total_source_materials=(),
        )
        field = np.asarray([result["B_cf"](mesh(0, 0, z)) for z in AXIS_Z], dtype=float)
    exact = np.asarray([exact_bz(z) for z in AXIS_Z])
    error = np.abs(field[:, 2] - exact) / exact
    return {
        "mesh_maxh_m": maxh,
        "order": order,
        "cell_count": mesh.ne,
        "vertex_count": mesh.nv,
        "ndof": result["fes"].ndof,
        "field_T": field.tolist(),
        "reference_bz_T": exact.tolist(),
        "maximum_relative_error": float(error.max()),
        "linear_residual": result.get("linear_residual"),
        "elapsed_seconds": time.perf_counter() - started,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxh", type=float, default=2.0)
    parser.add_argument("--order", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = solve(args.maxh, args.order)
    encoded = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
