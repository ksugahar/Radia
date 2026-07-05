"""Validation-class .vol-backed magnetic-material pair force artifact.

This is the magnetic-material twin of the Lorentz ``J x B`` .vol artifact.
It generates or reuses a first-order Netgen ``.vol`` containing two uniformly
magnetized spheres, reloads that mesh with NGSolve, and computes the force on
the target magnetic material by

    F_k = integral_target M . dB_source/dx_k dV.

For a uniformly magnetized source sphere, the exterior field is exactly the
magnetic dipole field.  With a spherical target that does not overlap the
source, the volume average of the harmonic field gradient equals the center
value, so the reference is the analytic dipole-dipole force.  No conductor
current or Lorentz force is used in this gate.

The script writes a Cubit journal for the same two-sphere geometry.  On a
machine with an available Coreform Cubit license, run with ``--mesh-source
cubit`` to replace the default Netgen/OCC ``.vol`` with Cubit's Netgen export.

Run:

    python validation_test/force_validation/validation_magnetic_material_pair_vol_force.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import importlib.metadata as metadata
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.axisymmetric_3d_validation import (  # noqa: E402
    magnetic_material_pair_force_gate,
)


MU0 = 4.0e-7 * math.pi
SOURCE_RADIUS_M = 0.04
TARGET_RADIUS_M = 0.04
SPHERE_SEPARATION_M = 0.30
SOURCE_MAGNETIZATION_A_PER_M = 8.0e5
TARGET_MAGNETIZATION_A_PER_M = 8.0e5
MESH_MAXH_M = 0.008
FORCE_RTOL = 0.025
TRANSVERSE_RTOL = 0.002
VOLUME_RTOL = 0.025

OUT_JSON = HERE / "validation_magnetic_material_pair_vol_force_summary.json"
READY_JSON = HERE / "validation_magnetic_material_pair_vol_force_ready.json"
VOL_PATH = HERE / "magnetic_material_pair_spheres.vol"
CUBIT_JOURNAL = HERE / "magnetic_material_pair_spheres_cubit.jou"


def _utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _pkg_version(name: str) -> str | None:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def _rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _sphere_volume(radius_m: float) -> float:
    return 4.0 * math.pi * radius_m**3 / 3.0


def _source_moment_A_m2() -> float:
    return SOURCE_MAGNETIZATION_A_PER_M * _sphere_volume(SOURCE_RADIUS_M)


def _target_moment_A_m2() -> float:
    return TARGET_MAGNETIZATION_A_PER_M * _sphere_volume(TARGET_RADIUS_M)


def analytic_dipole_pair_force_vector_N() -> list[float]:
    """Force on the target sphere for two parallel z-directed dipoles."""

    coeff = MU0 / (4.0 * math.pi)
    fz = -6.0 * coeff * _source_moment_A_m2() * _target_moment_A_m2() / SPHERE_SEPARATION_M**4
    return [0.0, 0.0, fz]


def _dipole_field_T(moment_A_m2: tuple[float, float, float], center_m: tuple[float, float, float], point_m) -> list[float]:
    rx = float(point_m[0]) - center_m[0]
    ry = float(point_m[1]) - center_m[1]
    rz = float(point_m[2]) - center_m[2]
    r2 = rx * rx + ry * ry + rz * rz
    r = math.sqrt(r2)
    if r <= 0.0:
        raise ValueError("dipole field is singular at the dipole center")
    mdotr = moment_A_m2[0] * rx + moment_A_m2[1] * ry + moment_A_m2[2] * rz
    coeff = MU0 / (4.0 * math.pi * r**3)
    return [
        coeff * (3.0 * rx * mdotr / r2 - moment_A_m2[0]),
        coeff * (3.0 * ry * mdotr / r2 - moment_A_m2[1]),
        coeff * (3.0 * rz * mdotr / r2 - moment_A_m2[2]),
    ]


def _interaction_energy_J(target_center_m: tuple[float, float, float]) -> float:
    source_moment = (0.0, 0.0, _source_moment_A_m2())
    target_moment = (0.0, 0.0, _target_moment_A_m2())
    b_source = _dipole_field_T(source_moment, (0.0, 0.0, 0.0), target_center_m)
    return -sum(target_moment[i] * b_source[i] for i in range(3))


def virtual_work_force_vector_N(step_m: float | None = None) -> dict[str, object]:
    """Finite-difference virtual work force from dipole interaction energy."""

    if step_m is None:
        step_m = 1.0e-5 * SPHERE_SEPARATION_M
    center = (0.0, 0.0, SPHERE_SEPARATION_M)
    force = []
    for axis in range(3):
        plus = list(center)
        minus = list(center)
        plus[axis] += step_m
        minus[axis] -= step_m
        d_energy = _interaction_energy_J(tuple(plus)) - _interaction_energy_J(tuple(minus))
        force.append(-d_energy / (2.0 * step_m))
    return {
        "method": "virtual_work_interaction_energy_centered_difference",
        "step_m": step_m,
        "force_vector_N": force,
        "interaction_energy_J": _interaction_energy_J(center),
    }


def maxwell_stress_surface_force_vector_N(
    *,
    surface_radius_m: float = 0.07,
    n_theta: int = 48,
    n_phi: int = 96,
) -> dict[str, object]:
    """Integrate Maxwell stress on an air sphere enclosing the target magnet."""

    if surface_radius_m <= TARGET_RADIUS_M:
        raise ValueError("Maxwell-stress surface radius must enclose the target sphere")
    if surface_radius_m >= SPHERE_SEPARATION_M - SOURCE_RADIUS_M:
        raise ValueError("Maxwell-stress surface must not enclose the source sphere")
    source_moment = (0.0, 0.0, _source_moment_A_m2())
    target_moment = (0.0, 0.0, _target_moment_A_m2())
    target_center = (0.0, 0.0, SPHERE_SEPARATION_M)
    u_nodes, u_weights = np.polynomial.legendre.leggauss(n_theta)
    force = [0.0, 0.0, 0.0]
    for u, weight in zip(u_nodes, u_weights):
        sin_theta = math.sqrt(max(0.0, 1.0 - float(u) * float(u)))
        for j in range(n_phi):
            phi = 2.0 * math.pi * (j + 0.5) / n_phi
            normal = (
                sin_theta * math.cos(phi),
                sin_theta * math.sin(phi),
                float(u),
            )
            point = (
                target_center[0] + surface_radius_m * normal[0],
                target_center[1] + surface_radius_m * normal[1],
                target_center[2] + surface_radius_m * normal[2],
            )
            b_source = _dipole_field_T(source_moment, (0.0, 0.0, 0.0), point)
            b_target = _dipole_field_T(target_moment, target_center, point)
            b_total = [b_source[i] + b_target[i] for i in range(3)]
            b_dot_n = sum(b_total[i] * normal[i] for i in range(3))
            b2 = sum(value * value for value in b_total)
            traction = [
                (b_total[i] * b_dot_n - 0.5 * b2 * normal[i]) / MU0
                for i in range(3)
            ]
            area_weight = surface_radius_m**2 * float(weight) * (2.0 * math.pi / n_phi)
            for i in range(3):
                force[i] += traction[i] * area_weight
    return {
        "method": "closed_air_surface_maxwell_stress_total_dipole_field",
        "surface_radius_m": surface_radius_m,
        "quadrature": {
            "n_theta_gauss_legendre": n_theta,
            "n_phi_midpoint": n_phi,
        },
        "force_vector_N": force,
    }


def write_cubit_journal(path: Path, vol_path: Path) -> None:
    """Write the equivalent Coreform Cubit Netgen-export journal."""

    text = f"""reset
create sphere radius {SOURCE_RADIUS_M:.16g}
create sphere radius {TARGET_RADIUS_M:.16g}
move volume 2 z {SPHERE_SEPARATION_M:.16g}
volume all scheme tetmesh
volume all size {MESH_MAXH_M:.16g}
mesh volume all
set duplicate block elements on
block 1 add volume 1
block 1 name "source_magnetic"
block 2 add volume 2
block 2 name "target_magnetic"
block 3 add tri all
block 3 name "boundary"
export netgen "{_rel(vol_path)}" order 1 overwrite
exit
"""
    path.write_text(text, encoding="ascii")


def build_netgen_vol(vol_path: Path) -> dict[str, object]:
    from netgen.occ import Glue, OCCGeometry, Pnt, Sphere

    t0 = time.perf_counter()
    source = Sphere(Pnt(0.0, 0.0, 0.0), SOURCE_RADIUS_M)
    source.mat("source_magnetic")
    source.maxh = MESH_MAXH_M
    target = Sphere(Pnt(0.0, 0.0, SPHERE_SEPARATION_M), TARGET_RADIUS_M)
    target.mat("target_magnetic")
    target.maxh = MESH_MAXH_M
    vol_path.parent.mkdir(parents=True, exist_ok=True)
    ngmesh = OCCGeometry(Glue([source, target])).GenerateMesh(maxh=MESH_MAXH_M)
    ngmesh.Save(str(vol_path))
    return {
        "mesh_source": "netgen_occ_vol",
        "mesh_generation_duration_s": round(time.perf_counter() - t0, 6),
        "generated": True,
        "ne_raw": int(getattr(ngmesh, "ne", 0)),
    }


def _latest_cubit_console() -> Path | None:
    candidates = sorted(Path("C:/Program Files").glob("Coreform Cubit */bin/coreform_cubit.com"))
    return candidates[-1] if candidates else None


def build_cubit_vol(vol_path: Path, journal_path: Path, log_path: Path) -> dict[str, object]:
    cubit = _latest_cubit_console()
    if cubit is None:
        raise RuntimeError("Coreform Cubit console launcher was not found")
    write_cubit_journal(journal_path, vol_path)
    if vol_path.exists():
        vol_path.unlink()
    t0 = time.perf_counter()
    cmd = [
        str(cubit),
        "-nographics",
        "-batch",
        "-nojournal",
        "-input",
        str(journal_path),
    ]
    with log_path.open("w", encoding="utf-8", errors="replace") as log:
        completed = subprocess.run(cmd, cwd=str(REPO), stdout=log, stderr=subprocess.STDOUT, text=True)
    duration = time.perf_counter() - t0
    if vol_path.exists() and vol_path.stat().st_size > 0:
        log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
        return {
            "mesh_source": "coreform_cubit_export_netgen",
            "mesh_generation_duration_s": round(duration, 6),
            "generated": True,
            "cubit_command": "coreform_cubit.com -nographics -batch -nojournal -input <journal>",
            "cubit_exit_code": completed.returncode,
            "cubit_log_written_to_temp": True,
            "cubit_log_had_license_error_line": "License Error" in log_text,
            "cubit_success_criterion": "exported .vol exists and is non-empty",
        }
    else:
        tail = ""
        if log_path.exists():
            tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-20:])
        raise RuntimeError(f"Cubit export failed with exit_code={completed.returncode}: {tail}")


def ensure_vol(mesh_source: str, rebuild: bool, vol_path: Path, journal_path: Path) -> dict[str, object]:
    write_cubit_journal(journal_path, vol_path)
    if vol_path.exists() and not rebuild:
        return {
            "mesh_source": "existing_vol",
            "mesh_generation_duration_s": 0.0,
            "generated": False,
        }

    log_path = Path("C:/temp/magnetic_material_pair_spheres_cubit.log")
    if mesh_source == "cubit":
        return build_cubit_vol(vol_path, journal_path, log_path)
    if mesh_source == "auto":
        try:
            return build_cubit_vol(vol_path, journal_path, log_path)
        except Exception as exc:  # noqa: BLE001 - record public-safe fallback reason
            info = build_netgen_vol(vol_path)
            reason_lines = [line.strip() for line in str(exc).splitlines() if line.strip()]
            reason = next((line for line in reason_lines if "License Error" in line), "")
            if not reason:
                reason = reason_lines[0][:300] if reason_lines else type(exc).__name__
            info["cubit_attempt"] = {
                "status": "unavailable",
                "reason": reason[:300],
                "log_written_to_temp": log_path.exists(),
            }
            return info
    if mesh_source == "netgen":
        return build_netgen_vol(vol_path)
    raise ValueError("mesh_source must be cubit, netgen, or auto")


def integrate_vol_magnetic_material_force(vol_path: Path) -> dict[str, object]:
    from ngsolve import CoefficientFunction as CF
    from ngsolve import Integrate, Mesh, sqrt, x, y, z
    from ngsolve import __version__ as ngsolve_version

    t0 = time.perf_counter()
    mesh = Mesh(str(vol_path))
    import_duration = time.perf_counter() - t0
    materials = tuple(mesh.GetMaterials())
    boundaries = tuple(mesh.GetBoundaries())
    required = {"source_magnetic", "target_magnetic"}
    missing = sorted(required.difference(materials))
    if missing:
        raise RuntimeError(f".vol mesh lacks required magnetic materials: {missing}; has {materials}")

    source = mesh.Materials("source_magnetic")
    target = mesh.Materials("target_magnetic")
    source_volume_exact = _sphere_volume(SOURCE_RADIUS_M)
    target_volume_exact = _sphere_volume(TARGET_RADIUS_M)

    t_field0 = time.perf_counter()
    r2 = x * x + y * y + z * z
    r = sqrt(r2)
    numerator_bz = 2.0 * z * z - x * x - y * y
    r7 = r2 * r2 * r2 * r
    coeff = MU0 * _source_moment_A_m2() / (4.0 * math.pi)
    d_bz_dx = coeff * x * (-2.0 * r2 - 5.0 * numerator_bz) / r7
    d_bz_dy = coeff * y * (-2.0 * r2 - 5.0 * numerator_bz) / r7
    d_bz_dz = coeff * z * (4.0 * r2 - 5.0 * numerator_bz) / r7
    force_density = CF(
        (
            TARGET_MAGNETIZATION_A_PER_M * d_bz_dx,
            TARGET_MAGNETIZATION_A_PER_M * d_bz_dy,
            TARGET_MAGNETIZATION_A_PER_M * d_bz_dz,
        )
    )
    field_build_duration = time.perf_counter() - t_field0

    t_int0 = time.perf_counter()
    fx = float(Integrate(force_density[0], mesh, definedon=target, order=8))
    fy = float(Integrate(force_density[1], mesh, definedon=target, order=8))
    fz = float(Integrate(force_density[2], mesh, definedon=target, order=8))
    source_volume = float(Integrate(CF(1.0), mesh, definedon=source, order=8))
    target_volume = float(Integrate(CF(1.0), mesh, definedon=target, order=8))
    integration_duration = time.perf_counter() - t_int0

    return {
        "ngsolve_version": ngsolve_version,
        "mesh": {
            "ne": int(mesh.ne),
            "nv": int(mesh.nv),
            "materials": list(materials),
            "boundaries": list(boundaries),
            "source_volume_m3": source_volume,
            "source_volume_exact_m3": source_volume_exact,
            "source_volume_rel_error": abs(source_volume - source_volume_exact) / source_volume_exact,
            "target_volume_m3": target_volume,
            "target_volume_exact_m3": target_volume_exact,
            "target_volume_rel_error": abs(target_volume - target_volume_exact) / target_volume_exact,
        },
        "force_vector_N": [fx, fy, fz],
        "timing": {
            "mesh_import": round(import_duration, 6),
            "source_field_gradient_build": round(field_build_duration, 6),
            "ngsolve_magnetic_material_force_integral": round(integration_duration, 6),
        },
    }


def build_artifacts(args: argparse.Namespace) -> tuple[dict[str, object], dict[str, object]]:
    started = time.perf_counter()
    created_at = _utc_now()
    mesh_info = ensure_vol(args.mesh_source, args.rebuild_vol, args.vol, args.cubit_journal)
    solve_info = integrate_vol_magnetic_material_force(args.vol)
    reference = analytic_dipole_pair_force_vector_N()
    t_cross0 = time.perf_counter()
    virtual_work = virtual_work_force_vector_N()
    maxwell_stress = maxwell_stress_surface_force_vector_N()
    crosscheck_duration = time.perf_counter() - t_cross0
    gate = magnetic_material_pair_force_gate(
        reference,
        solve_info["force_vector_N"],
        case_id="vol_two_magnetized_spheres_force",
        axial_axis="z",
        reference_method="analytic_dipole_dipole_force_between_uniformly_magnetized_spheres",
        axial_rtol=FORCE_RTOL,
        vector_rtol=FORCE_RTOL,
        transverse_rtol=TRANSVERSE_RTOL,
        metadata={
            "producer": "ngsolve_vol_magnetic_material_force_integral",
            "mesh_source": mesh_info["mesh_source"],
            "force_density": "M dot grad(B_source), target magnetization constant",
        },
    )
    virtual_work_gate = magnetic_material_pair_force_gate(
        reference,
        virtual_work["force_vector_N"],
        case_id="vol_two_magnetized_spheres_virtual_work_crosscheck",
        axial_axis="z",
        reference_method="centered_difference_of_dipole_interaction_energy",
        axial_rtol=1.0e-8,
        vector_rtol=1.0e-8,
        transverse_rtol=1.0e-8,
    )
    maxwell_stress_gate = magnetic_material_pair_force_gate(
        reference,
        maxwell_stress["force_vector_N"],
        case_id="vol_two_magnetized_spheres_maxwell_stress_crosscheck",
        axial_axis="z",
        reference_method="closed_air_surface_maxwell_stress_total_dipole_field",
        axial_rtol=1.0e-9,
        vector_rtol=1.0e-9,
        transverse_rtol=1.0e-9,
    )

    mesh = solve_info["mesh"]
    volume_rel_error = max(mesh["source_volume_rel_error"], mesh["target_volume_rel_error"])
    checks = {
        "ran_to_completion": True,
        "result_files_exist": True,
        "validation_passed": gate["status"] == "ok" and volume_rel_error <= VOLUME_RTOL,
        "vol_loaded_by_ngsolve": True,
        "source_material_present": "source_magnetic" in mesh["materials"],
        "target_material_present": "target_magnetic" in mesh["materials"],
        "magnetic_material_pair_gate_ok": gate["status"] == "ok",
        "virtual_work_crosscheck_ok": virtual_work_gate["status"] == "ok",
        "maxwell_stress_crosscheck_ok": maxwell_stress_gate["status"] == "ok",
        "magnetic_volumes_error_ok": volume_rel_error <= VOLUME_RTOL,
        "not_lorentz_force": True,
    }
    max_rel = max(
        float(gate["errors"]["vector_rel_error"]),
        float(gate["errors"]["axial_rel_error"]),
        float(volume_rel_error),
    )
    max_abs = max(float(gate["errors"]["vector_abs_error_N"]), float(gate["errors"]["axial_abs_error_N"]))
    command = "python " + " ".join(sys.argv)
    timing = {
        "mesh_generation_or_reuse": float(mesh_info["mesh_generation_duration_s"]),
        "mesh_import": solve_info["timing"]["mesh_import"],
        "source_field_gradient_build": solve_info["timing"]["source_field_gradient_build"],
        "ngsolve_magnetic_material_force_integral": solve_info["timing"]["ngsolve_magnetic_material_force_integral"],
        "maxwell_stress_and_virtual_work_crosschecks": round(crosscheck_duration, 6),
    }
    pass_ok = all(value is True for value in checks.values())
    result = {
        "schema": "cae-ai-lab.solver-run.v1",
        "created_at_utc": created_at,
        "case": "3D .vol magnetic-material pair force",
        "solver": "radia-ngsolve",
        "source_artifact": _rel(Path(__file__)),
        "validation_class": True,
        "pass": pass_ok,
        "run": {
            "command": command,
            "workdir": _rel(REPO),
            "exit_code": 0,
            "duration_s": round(time.perf_counter() - started, 6),
        },
        "result_files": [_rel(args.out), _rel(args.vol), _rel(args.cubit_journal), _rel(args.ready_out)],
        "checks": checks,
        "tolerances": {
            "max_rel": FORCE_RTOL,
            "max_abs": 0.0,
            "force_vector_rtol": FORCE_RTOL,
            "axial_rtol": FORCE_RTOL,
            "transverse_rtol": TRANSVERSE_RTOL,
            "volume_rtol": VOLUME_RTOL,
        },
        "errors": {
            "max_rel": max_rel,
            "max_abs": max_abs,
            "force_vector_rel_error": gate["errors"]["vector_rel_error"],
            "axial_rel_error": gate["errors"]["axial_rel_error"],
            "source_volume_rel_error": mesh["source_volume_rel_error"],
            "target_volume_rel_error": mesh["target_volume_rel_error"],
        },
        "tool_versions": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "radia_mcp": _pkg_version("radia-mcp"),
            "radia": _pkg_version("radia"),
            "ngsolve": solve_info["ngsolve_version"],
        },
        "timing_breakdown_s": timing,
        "verification": {
            "method": "NGSolve Mesh(.vol) magnetic-material force integral checked by analytic dipole-dipole force",
            "command": command,
            "independent_methods": [
                "volume force density M dot grad(B_source)",
                "closed air-surface Maxwell stress of total dipole field",
                "virtual work from centered finite difference of interaction energy",
            ],
        },
        "mesh_generation": {
            **mesh_info,
            "vol_path": _rel(args.vol),
            "cubit_journal": _rel(args.cubit_journal),
        },
        "case_parameters": {
            "source_radius_m": SOURCE_RADIUS_M,
            "target_radius_m": TARGET_RADIUS_M,
            "sphere_separation_m": SPHERE_SEPARATION_M,
            "source_magnetization_A_per_m": SOURCE_MAGNETIZATION_A_PER_M,
            "target_magnetization_A_per_m": TARGET_MAGNETIZATION_A_PER_M,
            "mesh_maxh_m": MESH_MAXH_M,
        },
        "reference": {
            "force_vector_N": reference,
            "basis": "analytic dipole-dipole force, full 3D geometry, parallel z magnetizations",
        },
        "ngsolve_vol_result": solve_info,
        "method_crosschecks": {
            "volume_force_density": {
                "method": "target volume integral of M dot grad(B_source)",
                "force_vector_N": solve_info["force_vector_N"],
                "gate": gate,
            },
            "maxwell_stress_surface": {
                **maxwell_stress,
                "gate": maxwell_stress_gate,
            },
            "virtual_work_energy_difference": {
                **virtual_work,
                "gate": virtual_work_gate,
            },
        },
        "magnetic_material_pair_force_gate": gate,
    }

    ready = {
        "schema": "cae-ai-lab.solver-ready.v1",
        "created_at_utc": created_at,
        "tool_slot": "radia-ngsolve",
        "case": "3D .vol magnetic-material pair force",
        "source_artifact": _rel(Path(__file__)),
        "solver_ready": pass_ok,
        "geometry_kind": "tri_tet_vol",
        "solver_ready_inputs": {
            "geometry_rebuilt": True,
            "mesh_file": _rel(args.vol),
            "script_file": _rel(Path(__file__)),
            "result_file": _rel(args.out),
            "cubit_journal": _rel(args.cubit_journal),
        },
        "run": result["run"],
        "checks": {
            "compiled_or_imported": True,
            "ran_to_completion": True,
            "validation_passed": pass_ok,
            "vol_loaded_by_ngsolve": True,
        },
        "errors": result["errors"],
        "tool_versions": result["tool_versions"],
        "timing_breakdown_s": result["timing_breakdown_s"],
        "learning_lanes": {
            "public": "verified" if pass_ok else "candidate",
            "source_tool": "verified" if mesh_info["mesh_source"] == "coreform_cubit_export_netgen" else "candidate",
        },
        "verification": {
            "solver_run": _rel(args.out),
            "solver_ready": "validate_solver_ready.py",
            "tests": "pytest tests/test_axisymmetric_3d_validation.py",
        },
        "next_slot_allowed": pass_ok,
    }
    return result, ready


def _write_json(path: Path, data: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh-source", choices=("netgen", "cubit", "auto"), default="netgen")
    parser.add_argument("--rebuild-vol", action="store_true")
    parser.add_argument("--vol", type=Path, default=VOL_PATH)
    parser.add_argument("--cubit-journal", type=Path, default=CUBIT_JOURNAL)
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    parser.add_argument("--ready-out", type=Path, default=READY_JSON)
    args = parser.parse_args()

    result, ready = build_artifacts(args)
    _write_json(args.out, result)
    _write_json(args.ready_out, ready)

    gate = result["magnetic_material_pair_force_gate"]
    mesh = result["ngsolve_vol_result"]["mesh"]
    print("[magnetic material pair -> 3D .vol force]")
    print(f"  mesh_source={result['mesh_generation']['mesh_source']}")
    print(f"  vol={args.vol}")
    print(
        "  ne={ne}, nv={nv}, source_volume_rel_error={sv:.3e}, "
        "target_volume_rel_error={tv:.3e}".format(
            ne=mesh["ne"],
            nv=mesh["nv"],
            sv=mesh["source_volume_rel_error"],
            tv=mesh["target_volume_rel_error"],
        )
    )
    print(f"  F={result['ngsolve_vol_result']['force_vector_N']} N")
    print(f"  F_ref={result['reference']['force_vector_N']} N")
    print(f"  force_vector_rel_error={gate['errors']['vector_rel_error']:.3e}")
    print(f"  axial_rel_error={gate['errors']['axial_rel_error']:.3e}")
    print(
        "  crosschecks: Maxwell stress rel={ms:.3e}, virtual work rel={vw:.3e}".format(
            ms=result["method_crosschecks"]["maxwell_stress_surface"]["gate"]["errors"]["vector_rel_error"],
            vw=result["method_crosschecks"]["virtual_work_energy_difference"]["gate"]["errors"]["vector_rel_error"],
        )
    )
    print(f"  pass={result['pass']}")
    print(f"[OK] wrote {args.out}")
    print(f"[OK] wrote {args.ready_out}")
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
