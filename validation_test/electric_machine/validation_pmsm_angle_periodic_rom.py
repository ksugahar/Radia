"""Curved 8-pole/24-slot PMSM angle-ROM held-out validation.

The independent planar A-formulation uses a fixed linear FEM operator, a
smooth rotating PM remanence, and 24 physical winding slots.  Thirty-three
angle samples build the periodic ROM.  The interlaced 33 angles are reserved
for Maxwell-stress, direct virtual-work, flux-linkage, and ROM comparisons.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import tempfile
import time

import numpy as np
import ngsolve as ng


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PANEL_VALIDATION = REPO / "validation_test" / "panels"
if str(PANEL_VALIDATION) not in sys.path:
    sys.path.insert(0, str(PANEL_VALIDATION))

from build_test_motor_mesh import build_pmsm_mesh  # noqa: E402


def _motor_rom_api():
    source = os.environ.get("RADIA_MOTOR_ROM_SOURCE")
    if source:
        spec = importlib.util.spec_from_file_location("radia_motor_rom_current", source)
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module
    from radia import motor_rom
    return motor_rom


motor_rom = _motor_rom_api()
AnglePeriodicMotorROM = motor_rom.AnglePeriodicMotorROM
MotorPortContract = motor_rom.MotorPortContract
PeriodicAngleTable = motor_rom.PeriodicAngleTable

MU0 = 4.0e-7 * math.pi
DEFAULT_RESULT = HERE / "validation_pmsm_angle_periodic_rom_summary.json"


def _mapping_quality(mesh, integration_order=10):
    minimum_det = float("inf")
    minimum_scaled = float("inf")
    samples = 0
    invalid = 0
    for element in mesh.Elements(ng.VOL):
        trafo = mesh.GetTrafo(element)
        for point in ng.IntegrationRule(element.type, integration_order):
            jacobian = np.asarray(trafo(point).jacobi, dtype=float)
            determinant = float(np.linalg.det(jacobian))
            norm_product = float(np.prod(np.linalg.norm(jacobian, axis=0)))
            scaled = determinant / norm_product if norm_product > 0.0 else float("-inf")
            minimum_det = min(minimum_det, determinant)
            minimum_scaled = min(minimum_scaled, scaled)
            invalid += int(not np.isfinite(determinant) or determinant <= 0.0)
            samples += 1
    return {
        "curve_order": int(mesh.GetCurveOrder()),
        "integration_order": int(integration_order),
        "samples": samples,
        "minimum_jacobian": minimum_det,
        "minimum_scaled_jacobian": minimum_scaled,
        "nonpositive_jacobian_samples": invalid,
    }


def _sampled_maxwell_torque(mesh, potential, *, radius, stack_length, points):
    phase = (np.arange(points) + 0.371) * 2.0 * np.pi / points
    x = radius * np.cos(phase)
    y = radius * np.sin(phase)
    gradient = np.asarray(ng.grad(potential)(mesh(x, y)), dtype=float)
    if gradient.shape == (points, 2):
        gradient = gradient.T
    bx = gradient[1]
    by = -gradient[0]
    br = bx * np.cos(phase) + by * np.sin(phase)
    bphi = -bx * np.sin(phase) + by * np.cos(phase)
    return float(
        stack_length / MU0 * radius**2 * (2.0 * np.pi / points)
        * np.sum(br * bphi)
    )


def build_summary(
    *,
    mesh_path: Path,
    maxh: float = 0.004,
    geometry_order: int = 4,
    fes_order: int = 3,
    angle_samples: int = 33,
    maxwell_points: int = 8191,
) -> dict:
    started = time.perf_counter()
    mesh_path.parent.mkdir(parents=True, exist_ok=True)
    build_pmsm_mesh(
        str(mesh_path),
        maxh=maxh,
        curved_geometry=True,
        curve_order=geometry_order,
        physical_slots=24,
    )
    mesh = ng.Mesh(str(mesh_path))
    quality = _mapping_quality(mesh)

    fes = ng.H1(mesh, order=fes_order, dirichlet="outer")
    trial, test = fes.TnT()
    materials = set(mesh.GetMaterials())
    reluctivity = mesh.MaterialCF(
        {
            material: 1.0 / (MU0 * 1000.0) if "iron" in material else 1.0 / MU0
            for material in materials
        }
    )
    bilinear = ng.BilinearForm(fes, symmetric=True)
    bilinear += reluctivity * ng.InnerProduct(ng.grad(trial), ng.grad(test)) * ng.dx
    with ng.TaskManager():
        bilinear.Assemble()
        inverse = bilinear.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")

    phase_regions = (
        ("stator_ind_Ap", "stator_ind_An"),
        ("stator_ind_Bp", "stator_ind_Bn"),
        ("stator_ind_Cp", "stator_ind_Cn"),
    )
    turns = 100.0
    slot_area = 1.0e-4
    stack_length = 0.05
    current_density = turns / slot_area

    def solve(linear_form):
        result = ng.GridFunction(fes)
        result.vec.data = inverse * linear_form.vec
        return result

    def phase_source(index):
        positive, negative = phase_regions[index]
        form = ng.LinearForm(fes)
        form += current_density * test * ng.dx(definedon=mesh.Materials(positive))
        form += -current_density * test * ng.dx(definedon=mesh.Materials(negative))
        form.Assemble()
        return form

    phase_sources = tuple(phase_source(index) for index in range(3))
    phase_potentials = tuple(solve(source) for source in phase_sources)

    def flux_linkage(potential):
        values = []
        for positive, negative in phase_regions:
            integral = ng.Integrate(
                potential, mesh, definedon=mesh.Materials(positive)
            ) - ng.Integrate(
                potential, mesh, definedon=mesh.Materials(negative)
            )
            values.append(stack_length * turns / slot_area * integral)
        return np.asarray(values, dtype=float)

    inductance = np.column_stack(
        [flux_linkage(potential) for potential in phase_potentials]
    )
    inductance = 0.5 * (inductance + inductance.T)

    def pm_solution(angle):
        phi = ng.atan2(ng.y, ng.x)
        amplitude = 1.2 * ng.cos(4.0 * (phi - float(angle)))
        remanence = mesh.MaterialCF(
            {
                "rotor_pm": ng.CoefficientFunction(
                    (amplitude * ng.cos(phi), amplitude * ng.sin(phi))
                )
            },
            default=ng.CoefficientFunction((0.0, 0.0)),
        )
        grad_test = ng.grad(test)
        form = ng.LinearForm(fes)
        form += (1.0 / MU0) * (
            -remanence[1] * grad_test[0] + remanence[0] * grad_test[1]
        ) * ng.dx
        form.Assemble()
        potential = solve(form)
        pm_flux = flux_linkage(potential)
        pm_coenergy = 0.5 * stack_length * float(ng.InnerProduct(form.vec, potential.vec))
        return form, potential, pm_flux, pm_coenergy

    mechanical_period = 2.0 * np.pi / 4.0
    train_angles = np.linspace(0.0, mechanical_period, angle_samples, endpoint=False)
    pm_flux_samples = []
    cogging_samples = []
    for angle in train_angles:
        _, _, pm_flux, pm_coenergy = pm_solution(angle)
        pm_flux_samples.append(pm_flux)
        cogging_samples.append(pm_coenergy)
    cogging_samples = np.asarray(cogging_samples)
    cogging_samples -= np.mean(cogging_samples)

    inductance_samples = np.repeat(
        inductance[None, :, :], angle_samples, axis=0
    )
    resistance_samples = np.repeat(
        (0.5 * np.eye(3))[None, :, :], angle_samples, axis=0
    )
    motor = AnglePeriodicMotorROM(
        MotorPortContract(("A", "B", "C")),
        PeriodicAngleTable(train_angles, inductance_samples, mechanical_period),
        PeriodicAngleTable(train_angles, resistance_samples, mechanical_period),
        PeriodicAngleTable(train_angles, np.asarray(pm_flux_samples), mechanical_period),
        cogging_coenergy_J=PeriodicAngleTable(
            train_angles, cogging_samples, mechanical_period
        ),
        inertia_kg_m2=1.0e-3,
    )

    currents = np.asarray((5.0, -2.5, -2.5))
    holdout_angles = train_angles + 0.5 * mechanical_period / angle_samples
    virtual_delta = 1.0e-5
    rows = []
    for angle in holdout_angles:
        _, pm_potential, pm_flux, pm_coenergy = pm_solution(angle)
        total = ng.GridFunction(fes)
        total.vec.data = pm_potential.vec
        for current, phase_potential in zip(currents, phase_potentials):
            total.vec.data += float(current) * phase_potential.vec
        torque_maxwell = _sampled_maxwell_torque(
            mesh,
            total,
            radius=0.05020,
            stack_length=stack_length,
            points=maxwell_points,
        )
        torque_rom = motor.torque_components(angle, currents)["total"]
        torque_rom_virtual = motor.virtual_work_torque(
            angle, currents, delta_angle_rad=virtual_delta
        )
        _, _, flux_minus, coenergy_minus = pm_solution(angle - virtual_delta)
        _, _, flux_plus, coenergy_plus = pm_solution(angle + virtual_delta)
        total_coenergy_minus = (
            0.5 * float(currents @ inductance @ currents)
            + float(currents @ flux_minus)
            + coenergy_minus
        )
        total_coenergy_plus = (
            0.5 * float(currents @ inductance @ currents)
            + float(currents @ flux_plus)
            + coenergy_plus
        )
        torque_fem_virtual = (
            total_coenergy_plus - total_coenergy_minus
        ) / (2.0 * virtual_delta)
        flux_relative_error = float(
            np.linalg.norm(motor.pm_flux(angle) - pm_flux)
            / max(np.linalg.norm(pm_flux), np.finfo(float).tiny)
        )
        rows.append(
            {
                "rotor_angle_rad": float(angle),
                "pm_flux_linkage_Wb": pm_flux.tolist(),
                "pm_flux_relative_error": flux_relative_error,
                "pm_coenergy_J": float(pm_coenergy),
                "torque_maxwell_Nm": torque_maxwell,
                "torque_fem_virtual_work_Nm": torque_fem_virtual,
                "torque_rom_derivative_Nm": torque_rom,
                "torque_rom_virtual_work_Nm": torque_rom_virtual,
            }
        )

    maxwell = np.asarray([row["torque_maxwell_Nm"] for row in rows])
    fem_virtual = np.asarray([row["torque_fem_virtual_work_Nm"] for row in rows])
    rom_torque = np.asarray([row["torque_rom_derivative_Nm"] for row in rows])
    rom_virtual = np.asarray([row["torque_rom_virtual_work_Nm"] for row in rows])
    torque_scale = max(float(np.linalg.norm(maxwell)), np.finfo(float).tiny)
    metrics = {
        "pm_flux_max_relative_error": max(row["pm_flux_relative_error"] for row in rows),
        "maxwell_vs_rom_torque_relative_rmse": float(
            np.linalg.norm(maxwell - rom_torque) / torque_scale
        ),
        "maxwell_vs_fem_virtual_torque_relative_rmse": float(
            np.linalg.norm(maxwell - fem_virtual) / torque_scale
        ),
        "fem_virtual_vs_rom_torque_relative_rmse": float(
            np.linalg.norm(fem_virtual - rom_torque)
            / max(float(np.linalg.norm(fem_virtual)), np.finfo(float).tiny)
        ),
        "rom_derivative_vs_virtual_max_abs_Nm": float(
            np.max(np.abs(rom_torque - rom_virtual))
        ),
        "minimum_inductance_eigenvalue_H": float(
            np.min(np.linalg.eigvalsh(inductance))
        ),
        "maximum_torque_abs_Nm": float(np.max(np.abs(maxwell))),
    }
    thresholds = {
        "pm_flux_max_relative_error": 1.0e-10,
        "maxwell_vs_rom_torque_relative_rmse": 3.0e-3,
        "maxwell_vs_fem_virtual_torque_relative_rmse": 3.0e-3,
        "fem_virtual_vs_rom_torque_relative_rmse": 2.0e-7,
        "rom_derivative_vs_virtual_max_abs_Nm": 2.0e-8,
    }
    checks = {
        name: metrics[name] <= threshold
        for name, threshold in thresholds.items()
    }
    checks["positive_curved_mapping"] = (
        quality["nonpositive_jacobian_samples"] == 0
        and quality["minimum_scaled_jacobian"] > 1.0e-6
    )
    checks["positive_inductance"] = metrics["minimum_inductance_eigenvalue_H"] > 0.0

    return {
        "schema": "radia.motor.pmsm_angle_periodic_rom_validation.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "validation_host": platform.node(),
        "python_version": platform.python_version(),
        "ngsolve_version": ng.__version__,
        "numpy_version": np.__version__,
        "geometry": {
            "machine": "8-pole 24-physical-slot PMSM benchmark",
            "pm_model": "smooth radial remanence fundamental",
            "mesh_elements": int(mesh.ne),
            "mesh_vertices": int(mesh.nv),
            "maxh_m": float(maxh),
            "geometry_order": int(geometry_order),
            "fes_order": int(fes_order),
            "maxwell_circle_radius_m": 0.05020,
            "maxwell_circle_points": int(maxwell_points),
        },
        "mapping_quality": quality,
        "rom": {
            "angle_samples": int(angle_samples),
            "held_out_angles": int(len(holdout_angles)),
            "mechanical_period_rad": float(mechanical_period),
            "generalized_current_order": list(motor.ports.generalized_names),
            "current_operating_point_A": currents.tolist(),
            "passive": motor.diagnostics()["passive"],
            "includes_cogging_coenergy": True,
        },
        "metrics": metrics,
        "thresholds": thresholds,
        "checks": checks,
        "passed": all(checks.values()),
        "elapsed_s": time.perf_counter() - started,
        "holdout_rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mesh",
        type=Path,
        default=Path(tempfile.gettempdir()) / "radia_pmsm_angle_rom.vol",
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--maxh", type=float, default=0.004)
    parser.add_argument("--geometry-order", type=int, default=4)
    parser.add_argument("--fes-order", type=int, default=3)
    parser.add_argument("--angle-samples", type=int, default=33)
    parser.add_argument("--maxwell-points", type=int, default=8191)
    args = parser.parse_args()
    result = build_summary(
        mesh_path=args.mesh,
        maxh=args.maxh,
        geometry_order=args.geometry_order,
        fes_order=args.fes_order,
        angle_samples=args.angle_samples,
        maxwell_points=args.maxwell_points,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({
        "passed": result["passed"],
        "metrics": result["metrics"],
        "checks": result["checks"],
        "elapsed_s": result["elapsed_s"],
    }, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
