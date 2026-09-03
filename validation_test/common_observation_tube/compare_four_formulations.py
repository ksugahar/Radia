"""Compare four magnetostatic formulations on one observation tube.

The benchmark is a finite 40 mm square linear iron yoke, 60 mm long, with a
20 mm square beam aperture.  A uniform vertical source field excites the yoke.
All engines sample the same 6 mm-radius transverse circle along the aperture:

* legacy ESRF Radia 2023 boundary-element reference;
* BDM1 HDiv-MMM on an iron-only mesh;
* reduced-A with an HCurl correction on a finite physical air domain;
* Omega-reduced-Omega with a periodic Kelvin exterior.

The JSON report contains pointwise B residuals, integrated-field residuals,
and normal/skew complex multipoles through order six.  The CSV files retain
the longitudinal main-field and integrated multipole values for plotting.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from radia.accelerator_field_validation import (
    CurvilinearObservationTube,
    MagnetFieldEngine,
    circular_transverse_offsets,
    compare_magnetic_flux_density,
    sample_field_engine,
    transverse_multipole_spectrum,
)

MU0 = 4.0e-7 * math.pi
MODEL = {
    "outer_cross_section_m": [0.040, 0.040],
    "aperture_m": [0.020, 0.020],
    "iron_length_m": 0.060,
    "mu_r": 100.0,
    "applied_b_t": 0.100,
}


def _iron_bars():
    from netgen.occ import Box, Pnt

    def box(lo, hi):
        solid = Box(Pnt(*lo), Pnt(*hi))
        solid.mat("iron")
        solid.maxh = 0.006
        for face in solid.faces:
            face.name = "iron_surf"
        return solid

    return [
        box((-0.020, 0.010, -0.030), (0.020, 0.020, 0.030)),
        box((-0.020, -0.020, -0.030), (0.020, -0.010, 0.030)),
        box((0.010, -0.010, -0.030), (0.020, 0.010, 0.030)),
        box((-0.020, -0.010, -0.030), (-0.010, 0.010, 0.030)),
    ]


def _iron_shape():
    from netgen.occ import Glue

    shape = Glue(_iron_bars())
    for solid in shape.solids:
        solid.name = "iron"
    return shape


def _generate_mesh(shape, maxh):
    import ngsolve as ng
    from netgen.occ import OCCGeometry

    with ng.TaskManager():
        ngmesh = OCCGeometry(shape).GenerateMesh(maxh=maxh, grading=0.4)
    return ng.Mesh(ngmesh)


def build_iron_mesh():
    return _generate_mesh(_iron_shape(), 0.006)


def build_reduced_a_mesh(radius=0.12):
    from netgen.occ import Glue, Pnt, Sphere

    iron = _iron_shape()
    outer = Sphere(Pnt(0, 0, 0), radius)
    outer.maxh = 0.018
    for face in outer.faces:
        face.name = "outer"
    air = outer - iron
    air.mat("air")
    shape = Glue([air, iron])
    return _generate_mesh(shape, 0.018)


def build_omega_kelvin_mesh(radius=0.12, offset=(0.30, 0.0, 0.0)):
    import ngsolve as ng
    from netgen.occ import Glue, IdentificationType, OCCGeometry, Pnt, Sphere, Vertex

    iron = _iron_shape()
    inner = Sphere(Pnt(0, 0, 0), radius)
    inner.maxh = 0.018
    for face in inner.faces:
        face.name = "kelvin_int"
    air_inner = inner - iron
    air_inner.mat("air_inner")

    outer = Sphere(Pnt(*offset), radius)
    outer.maxh = 0.020
    outer.mat("air_outer")
    for face in outer.faces:
        face.name = "kelvin_ext"
    gnd = Vertex(Pnt(*offset))
    gnd.name = "GND"

    shape = Glue([air_inner, iron, outer, gnd])
    for solid in shape.solids:
        if solid.name not in ("iron", "air_inner", "air_outer"):
            center = solid.center
            solid.name = "air_outer" if center.x > 0.18 else "air_inner"
    inner_faces = [face for face in shape.faces if face.name == "kelvin_int"]
    outer_faces = [face for face in shape.faces if face.name == "kelvin_ext"]
    if len(inner_faces) != 1 or len(outer_faces) != 1:
        raise RuntimeError(
            "Kelvin geometry must preserve one inner and one outer sphere face"
        )
    inner_faces[0].Identify(
        outer_faces[0], "kelvin_periodic", IdentificationType.PERIODIC
    )
    with ng.TaskManager():
        ngmesh = OCCGeometry(shape).GenerateMesh(maxh=0.020, grading=0.4)
    return ng.Mesh(ngmesh)


def observation_tube(stations=41, circle_points=24, radius=0.006):
    s = np.linspace(-0.055, 0.055, int(stations))
    return CurvilinearObservationTube(
        station_s=s,
        center=np.column_stack((np.zeros_like(s), np.zeros_like(s), s)),
        tangent=np.tile([0.0, 0.0, 1.0], (s.size, 1)),
        normal=np.tile([1.0, 0.0, 0.0], (s.size, 1)),
        binormal=np.tile([0.0, 1.0, 0.0], (s.size, 1)),
        transverse_offsets=circular_transverse_offsets(radius, circle_points),
    )


def _evaluate_cf(mesh, coefficient, points):
    values = np.empty_like(points, dtype=float)
    for index, point in enumerate(points):
        mip = mesh(*map(float, point))
        value = coefficient(mip)
        values[index] = [float(value[component]) for component in range(3)]
    return values


def _legacy_python38():
    explicit = os.environ.get("RADIA_LEGACY_PYTHON38")
    if explicit:
        executable = Path(explicit).expanduser().resolve()
    else:
        try:
            completed = subprocess.run(
                ["uv", "python", "find", "3.8"],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            raise RuntimeError(
                "Python 3.8 is required for the legacy Radia oracle; set "
                "RADIA_LEGACY_PYTHON38 to its python.exe"
            ) from exc
        executable = Path(completed.stdout.strip()).resolve()
    if not executable.is_file():
        raise RuntimeError(f"legacy Python 3.8 executable not found: {executable}")
    return executable


def radia_reference_engine(tube, output_dir, segmentation):
    points = tube.global_points().reshape(-1, 3)
    request_path = output_dir / "radia_reference_request.json"
    request_path.write_text(
        json.dumps(
            {
                "points_m": points.tolist(),
                "mu_r": MODEL["mu_r"],
                "applied_b_t": MODEL["applied_b_t"],
            }
        ),
        encoding="utf-8",
    )
    script = Path(__file__).with_name("radia_reference_py38.py")
    levels = sorted({max(1, int(segmentation) - 2), int(segmentation)})
    if len(levels) != 2:
        raise ValueError("radia_segmentation must provide two convergence levels")
    legacy_python = _legacy_python38()
    records = []
    fields = []
    for level in levels:
        result_path = output_dir / f"radia_reference_seg{level}.json"
        command = [
            str(legacy_python),
            str(script),
            "--input-json",
            str(request_path),
            "--output-json",
            str(result_path),
            "--segmentation",
            str(level),
        ]
        completed = subprocess.run(
            command, cwd=script.parent, check=True, capture_output=True, text=True
        )
        raw = json.loads(result_path.read_text(encoding="utf-8"))
        field = np.asarray(raw.pop("b_t"), dtype=float)
        if field.shape != points.shape:
            raise RuntimeError(
                f"legacy Radia returned {field.shape}, expected {points.shape}"
            )
        fields.append(field)
        records.append(raw)
    selected = fields[-1]
    convergence = []
    for index in range(1, len(fields)):
        difference = fields[index] - fields[index - 1]
        scale = max(float(np.sqrt(np.mean(fields[index] ** 2))), 1.0e-30)
        convergence.append(
            {
                "coarse_segmentation": levels[index - 1],
                "fine_segmentation": levels[index],
                "relative_rms_change": float(np.sqrt(np.mean(difference**2)) / scale),
                "maximum_vector_change_t": float(
                    np.max(np.linalg.norm(difference, axis=1))
                ),
            }
        )
    return MagnetFieldEngine("radia_reference", lambda query: selected.copy()), {
        "field_selection": "finest raw legacy-Radia discretization",
        "selected_segmentation": levels[-1],
        "levels": records,
        "convergence": convergence,
    }


def hdiv_mmm_engine():
    import ngsolve as ng

    from radia import vim

    mesh = build_iron_mesh()
    h0 = MODEL["applied_b_t"] / MU0
    with ng.TaskManager():
        result = vim.Solve(
            mesh,
            mu_r=MODEL["mu_r"],
            H_ext=ng.CF((0.0, h0, 0.0)),
            order=1,
            tol=1.0e-9,
            maxit=4000,
        )

    def evaluate(points):
        demagnetizing_h = vim.FieldFromSolution(result, points, algorithm="direct")
        total_h = demagnetizing_h + np.array([0.0, h0, 0.0])
        return MU0 * total_h

    return MagnetFieldEngine("hdiv_mmm", evaluate), {
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "ndof": int(result["ndof"]),
        "iterations": int(result["iters"]),
        "field_evaluator": result.get("field_evaluator_stats", {}),
    }


def reduced_a_engine(order=2):
    import ngsolve as ng

    from radia.vector_potential_solver import VectorPotentialSolver

    mesh = build_reduced_a_mesh()
    solver = VectorPotentialSolver(
        mesh, iron_domains="iron", mu_r=MODEL["mu_r"], order=order
    )
    solver.set_source_cf(ng.CF((0.0, MODEL["applied_b_t"], 0.0)))
    with ng.TaskManager():
        # ``.*`` would also match the named iron-air interfaces and would
        # incorrectly clamp tangential A_r there, eliminating the reaction
        # field.  Only the true truncation boundary is essential.
        solution = solver.solve_linear(dirichlet="outer", eps=1.0e-2, solver="direct")
    b_cf = solver.get_B()
    return MagnetFieldEngine(
        "reduced_a", lambda points: _evaluate_cf(mesh, b_cf, points)
    ), {
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "order": int(order),
        "ndof": int(solution.space.ndof),
        "outer_radius_m": 0.12,
    }


def omega_reduced_omega_engine(order=2):
    import ngsolve as ng
    from ngsolve import (
        BND,
        CF,
        H1,
        BilinearForm,
        GridFunction,
        LinearForm,
        Periodic,
        ds,
        dx,
        grad,
        specialcf,
        y,
    )

    from radia.kelvin_source import build_material_cf, kelvin_mu_factor_3d_cf

    radius = 0.12
    offset = (0.30, 0.0, 0.0)
    mesh = build_omega_kelvin_mesh(radius, offset)
    mesh.Curve(order)
    mu = build_material_cf(
        mesh,
        MU0,
        kelvin_mu_factor_3d_cf(center=offset, R=radius),
        outer_keyword="air_outer",
        overrides={"iron": MODEL["mu_r"] * MU0},
    )
    fes = Periodic(H1(mesh, order=order, dirichlet_bbnd="GND"))
    trial, test = fes.TnT()
    bilinear = BilinearForm(mu * grad(trial) * grad(test) * dx, symmetric=True)

    h0 = MODEL["applied_b_t"] / MU0
    omega_source = h0 * y
    source_b = CF((0.0, MODEL["applied_b_t"], 0.0))
    lift = GridFunction(fes, name="omega_source_lift")
    lift.Set(omega_source, BND, mesh.Boundaries("iron_surf"))
    linear = LinearForm(fes)
    linear += mu * grad(lift) * grad(test) * dx("air_inner")
    normal = -specialcf.normal(mesh.dim)
    linear += (normal * source_b) * test * ds("iron_surf")
    solution = GridFunction(fes, name="omega_total_reduced")
    with ng.TaskManager():
        bilinear.Assemble()
        linear.Assemble()
        solution.vec.data = (
            bilinear.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * linear.vec
        )

    air_space = H1(mesh, order=order, definedon="air_inner|air_outer")
    air_lift = GridFunction(air_space, name="omega_source_extension")
    air_lift.Set(omega_source, BND, mesh.Boundaries("iron_surf"))
    b_air = MU0 * (grad(solution) - grad(air_lift)) + source_b
    return MagnetFieldEngine(
        "omega_reduced_omega", lambda points: _evaluate_cf(mesh, b_air, points)
    ), {
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "order": int(order),
        "ndof": int(fes.ndof),
        "free_dofs": int(sum(bool(value) for value in fes.FreeDofs())),
        "kelvin_radius_m": radius,
        "kelvin_offset_m": list(offset),
    }


def _complex_rows(coefficients, s, reference_radius):
    integrated = np.trapezoid(coefficients, s, axis=0)
    main = integrated[0]
    rows = []
    for index, value in enumerate(integrated, start=1):
        normalized = (
            value * reference_radius ** (index - 1) / main
            if abs(main) > 1.0e-30
            else complex(float("nan"), float("nan"))
        )
        rows.append(
            {
                "order": index,
                "integrated_real_t_m_per_m_power": float(value.real),
                "integrated_imag_t_m_per_m_power": float(value.imag),
                "normal_units_at_reference_radius": float(1.0e4 * normalized.real),
                "skew_units_at_reference_radius": float(1.0e4 * normalized.imag),
            }
        )
    return rows


def _write_profile_csv(path, tube, samples, multipoles):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        header = ["s_m"]
        for sample in samples:
            header.extend(
                (f"{sample.engine_name}_C1_real_T", f"{sample.engine_name}_C1_imag_T")
            )
        writer.writerow(header)
        for station, s_value in enumerate(tube.station_s):
            row = [float(s_value)]
            for sample in samples:
                coefficient = multipoles[sample.engine_name][station, 0]
                row.extend((float(coefficient.real), float(coefficient.imag)))
            writer.writerow(row)


def _write_multipole_csv(path, report):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "engine",
                "order",
                "integrated_real_t_m_per_m_power",
                "integrated_imag_t_m_per_m_power",
                "normal_units_at_reference_radius",
                "skew_units_at_reference_radius",
            ),
        )
        writer.writeheader()
        for engine, rows in report.items():
            for row in rows:
                writer.writerow({"engine": engine, **row})


def run(
    output_dir,
    *,
    stations=41,
    circle_points=24,
    radia_segmentation=3,
    fe_order=2,
    gate_profile="retained",
):
    import ngsolve
    import radia

    stations = int(stations)
    circle_points = int(circle_points)
    radia_segmentation = int(radia_segmentation)
    fe_order = int(fe_order)
    if stations < 2:
        raise ValueError("stations must be at least two")
    if circle_points < 12:
        raise ValueError("circle_points must be at least twelve for order-six fitting")
    if radia_segmentation < 2:
        raise ValueError("radia_segmentation must be at least two")
    if fe_order < 1:
        raise ValueError("fe_order must be at least one")
    if gate_profile not in {"retained", "smoke"}:
        raise ValueError("gate_profile must be 'retained' or 'smoke'")
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    tube = observation_tube(stations, circle_points)
    constructors = [
        lambda: radia_reference_engine(tube, output_dir, radia_segmentation),
        hdiv_mmm_engine,
        lambda: reduced_a_engine(fe_order),
        lambda: omega_reduced_omega_engine(fe_order),
    ]
    samples = []
    diagnostics = {}
    timings = {}
    for constructor in constructors:
        t0 = time.perf_counter()
        engine, detail = constructor()
        samples.append(sample_field_engine(engine, tube))
        diagnostics[engine.name] = detail
        timings[engine.name] = time.perf_counter() - t0

    comparison = compare_magnetic_flux_density(samples)
    multipoles = {
        sample.engine_name: transverse_multipole_spectrum(sample, maximum_order=6)
        for sample in samples
    }
    reference_radius = float(np.linalg.norm(tube.transverse_offsets[0]))
    integrated_multipoles = {
        name: _complex_rows(values, tube.station_s, reference_radius)
        for name, values in multipoles.items()
    }
    reference_rows = integrated_multipoles["radia_reference"]
    multipole_comparison = {}
    for name, rows in integrated_multipoles.items():
        if name == "radia_reference":
            continue
        multipole_comparison[name] = [
            {
                "order": row["order"],
                "delta_normal_units": float(
                    row["normal_units_at_reference_radius"]
                    - reference["normal_units_at_reference_radius"]
                ),
                "delta_skew_units": float(
                    row["skew_units_at_reference_radius"]
                    - reference["skew_units_at_reference_radius"]
                ),
            }
            for row, reference in zip(rows, reference_rows)
        ]
    _write_profile_csv(
        output_dir / "longitudinal_main_field.csv", tube, samples, multipoles
    )
    _write_multipole_csv(
        output_dir / "integrated_multipoles.csv", integrated_multipoles
    )

    radia_pairs = {
        row["right"]: row
        for row in comparison["pairs"]
        if row["left"] == "radia_reference"
    }
    radia_convergence = diagnostics["radia_reference"]["convergence"][-1]
    checks = {
        "all_fields_finite": all(
            np.all(np.isfinite(sample.b_global)) for sample in samples
        ),
        "hdiv_solver_converged": 0 < diagnostics["hdiv_mmm"]["iterations"] < 4000,
    }
    accuracy_checks = {
        "radia_reference_relative_rms_change_below_3_percent": (
            radia_convergence["relative_rms_change"] < 0.03
        ),
        "radia_hdiv_relative_rms_below_5_percent": (
            radia_pairs["hdiv_mmm"]["relative_rms_error"] < 0.05
        ),
        "radia_reduced_a_relative_rms_below_5_percent": (
            radia_pairs["reduced_a"]["relative_rms_error"] < 0.05
        ),
        "radia_omega_relative_rms_below_5_percent": (
            radia_pairs["omega_reduced_omega"]["relative_rms_error"] < 0.05
        ),
        "all_radia_integrated_field_errors_below_7_percent": all(
            row["relative_integrated_error"] < 0.07 for row in radia_pairs.values()
        ),
    }
    if gate_profile == "retained":
        checks.update(accuracy_checks)
    report = {
        "schema": "radia.validation.common-observation-tube-four-formulations.v1",
        "status": "pass" if all(checks.values()) else "fail",
        "gate_profile": gate_profile,
        "machine": platform.node(),
        "python": platform.python_version(),
        "versions": {
            "python": platform.python_version(),
            "ngsolve": ngsolve.__version__,
            "radia": radia.__version__,
        },
        "model": MODEL,
        "observation_tube": {
            "station_range_m": [float(tube.station_s[0]), float(tube.station_s[-1])],
            "station_count": tube.station_count,
            "transverse_circle_radius_m": reference_radius,
            "transverse_point_count": tube.transverse_point_count,
            "local_components": ["x", "y", "s"],
        },
        "engine_diagnostics": diagnostics,
        "engine_runtime_s": timings,
        "field_comparison": comparison,
        "integrated_multipoles": integrated_multipoles,
        "multipole_comparison_to_radia": multipole_comparison,
        "checks": checks,
        "accuracy_diagnostics": accuracy_checks,
        "runtime_s": time.perf_counter() - started,
    }
    report_path = output_dir / "comparison_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "runtime_s": report["runtime_s"],
                "output": str(report_path),
                "checks": checks,
            },
            indent=2,
        )
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stations", type=int, default=41)
    parser.add_argument("--circle-points", type=int, default=24)
    parser.add_argument("--radia-segmentation", type=int, default=3)
    parser.add_argument("--fe-order", type=int, default=2)
    parser.add_argument(
        "--gate-profile", choices=("retained", "smoke"), default="retained")
    options = parser.parse_args()
    if options.stations < 2:
        parser.error("--stations must be at least two")
    if options.circle_points < 12:
        parser.error("--circle-points must be at least twelve")
    if options.radia_segmentation < 2:
        parser.error("--radia-segmentation must be at least two")
    if options.fe_order < 1:
        parser.error("--fe-order must be at least one")
    report = run(
        options.output_dir.resolve(),
        stations=options.stations,
        circle_points=options.circle_points,
        radia_segmentation=options.radia_segmentation,
        fe_order=options.fe_order,
        gate_profile=options.gate_profile,
    )
    raise SystemExit(0 if report["status"] == "pass" else 1)


if __name__ == "__main__":
    main()
