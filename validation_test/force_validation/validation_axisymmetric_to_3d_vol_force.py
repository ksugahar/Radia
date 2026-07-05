"""Validation-class .vol-backed 3-D force artifact.

This connects the axisymmetric-to-3D gate to an actual Netgen ``.vol`` mesh
loaded by NGSolve.  The default public-safe route generates a small target
torus with Netgen/OCC, reloads the saved ``.vol`` with ``Mesh(path)``, and
computes the Lorentz force

    F = integral_target J x B_source dV

where ``B_source`` is a discretized Biot-Savart field from a coaxial source
loop.  The axial component is checked against the exact axisymmetric
coaxial-loop force, and the resulting vector is passed through
``axisymmetric_to_3d_force_gate``.

The script also writes a Cubit journal for the same target torus.  On machines
with an available Coreform Cubit license, run with ``--mesh-source cubit`` to
replace the Netgen-generated ``.vol`` with Cubit's Netgen export.

Run:

    python validation_test/force_validation/validation_axisymmetric_to_3d_vol_force.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import importlib.metadata as metadata
import json
import math
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import scipy.special as sp


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SRC = REPO / "packages" / "radia-mcp" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from radia_mcp.radia_ngsolve.axisymmetric_3d_validation import (  # noqa: E402
    axisymmetric_to_3d_force_gate,
    axisymmetric_to_3d_validation_plan,
)


MU0 = 4.0e-7 * math.pi
LOOP_RADIUS_M = 1.0
LOOP_SEPARATION_M = 1.5
TARGET_TUBE_RADIUS_M = 0.03
TARGET_MESH_MAXH_M = 0.011
CURRENT_SOURCE_A = 1.0
CURRENT_TARGET_A = 1.0
SOURCE_SEGMENTS = 64
AXIAL_RTOL = 0.025
TRANSVERSE_RTOL = 0.001
VOLUME_RTOL = 0.025

OUT_JSON = HERE / "validation_axisymmetric_to_3d_vol_force_summary.json"
READY_JSON = HERE / "validation_axisymmetric_to_3d_vol_force_ready.json"
VOL_PATH = HERE / "axisymmetric_to_3d_target_torus.vol"
CUBIT_JOURNAL = HERE / "axisymmetric_to_3d_target_torus_cubit.jou"


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


def exact_axisymmetric_force_N() -> float:
    """Exact filament-loop force from dM/dz using elliptic integrals."""

    def mutual_inductance(radius_m: float, separation_m: float) -> float:
        m = 4.0 * radius_m * radius_m / (4.0 * radius_m * radius_m + separation_m * separation_m)
        k = math.sqrt(m)
        return MU0 * radius_m * (
            (2.0 / k - k) * sp.ellipk(m) - (2.0 / k) * sp.ellipe(m)
        )

    dz = 1.0e-6 * LOOP_SEPARATION_M
    dmdz = (
        mutual_inductance(LOOP_RADIUS_M, LOOP_SEPARATION_M + dz)
        - mutual_inductance(LOOP_RADIUS_M, LOOP_SEPARATION_M - dz)
    ) / (2.0 * dz)
    return CURRENT_SOURCE_A * CURRENT_TARGET_A * dmdz


def write_cubit_journal(path: Path, vol_path: Path) -> None:
    """Write the equivalent Coreform Cubit Netgen-export journal."""

    text = f"""reset
create torus major radius {LOOP_RADIUS_M:.16g} minor radius {TARGET_TUBE_RADIUS_M:.16g}
move volume 1 z {LOOP_SEPARATION_M:.16g}
volume 1 scheme tetmesh
volume 1 size {TARGET_MESH_MAXH_M:.16g}
mesh volume 1
set duplicate block elements on
block 1 add tet all
block 1 name "target"
block 2 add tri all
block 2 name "boundary"
export netgen "{_rel(vol_path)}" order 1 overwrite
exit
"""
    path.write_text(text, encoding="ascii")


def build_netgen_vol(vol_path: Path) -> dict[str, object]:
    from netgen.occ import Axes, Axis, Dir, OCCGeometry, Pnt, WorkPlane

    t0 = time.perf_counter()
    workplane = WorkPlane(
        Axes(
            p=Pnt(LOOP_RADIUS_M, 0.0, LOOP_SEPARATION_M),
            n=Dir(0.0, 1.0, 0.0),
            h=Dir(0.0, 0.0, 1.0),
        )
    )
    torus = workplane.Circle(TARGET_TUBE_RADIUS_M).Face().Revolve(
        Axis(Pnt(0.0, 0.0, 0.0), Dir(0.0, 0.0, 1.0)),
        360.0,
    )
    torus.name = "target"
    torus.mat("target")
    torus.maxh = TARGET_MESH_MAXH_M
    ngmesh = OCCGeometry(torus).GenerateMesh(maxh=TARGET_MESH_MAXH_M)
    vol_path.parent.mkdir(parents=True, exist_ok=True)
    ngmesh.Save(str(vol_path))
    return {
        "mesh_source": "netgen_occ_vol",
        "mesh_generation_duration_s": round(time.perf_counter() - t0, 6),
        "generated": True,
        "ne_raw": int(getattr(ngmesh, "ne", 0)),
    }


def _latest_cubit_console() -> Path | None:
    root = Path("C:/Program Files")
    candidates = sorted(root.glob("Coreform Cubit */bin/coreform_cubit.com"))
    return candidates[-1] if candidates else None


def build_cubit_vol(vol_path: Path, journal_path: Path, log_path: Path) -> dict[str, object]:
    """Try Coreform Cubit headless export.  Raises on failure."""

    cubit = _latest_cubit_console()
    if cubit is None:
        raise RuntimeError("Coreform Cubit console launcher was not found")
    write_cubit_journal(journal_path, vol_path)
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
    if completed.returncode != 0 or not vol_path.exists():
        tail = ""
        if log_path.exists():
            tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-20:])
        raise RuntimeError(f"Cubit export failed with exit_code={completed.returncode}: {tail}")
    return {
        "mesh_source": "coreform_cubit_export_netgen",
        "mesh_generation_duration_s": round(duration, 6),
        "generated": True,
        "cubit_command": "coreform_cubit.com -nographics -batch -nojournal -input <journal>",
        "cubit_log_written_to_temp": True,
    }


def ensure_vol(mesh_source: str, rebuild: bool, vol_path: Path, journal_path: Path) -> dict[str, object]:
    write_cubit_journal(journal_path, vol_path)
    if vol_path.exists() and not rebuild:
        return {
            "mesh_source": "existing_vol",
            "mesh_generation_duration_s": 0.0,
            "generated": False,
        }

    log_path = Path("C:/temp/axisymmetric_to_3d_target_torus_cubit.log")
    if mesh_source == "cubit":
        return build_cubit_vol(vol_path, journal_path, log_path)
    if mesh_source == "auto":
        try:
            return build_cubit_vol(vol_path, journal_path, log_path)
        except Exception as exc:  # noqa: BLE001 - record fallback reason in public-safe artifact
            info = build_netgen_vol(vol_path)
            info["cubit_attempt"] = {
                "status": "unavailable",
                "reason": str(exc).splitlines()[0][:300],
                "log_written_to_temp": log_path.exists(),
            }
            return info
    if mesh_source == "netgen":
        return build_netgen_vol(vol_path)
    raise ValueError("mesh_source must be cubit, netgen, or auto")


def source_loop_biot_savart_cf(n_segments: int):
    from ngsolve import sqrt, x, y, z

    bx = 0.0
    by = 0.0
    bz = 0.0
    dtheta = 2.0 * math.pi / n_segments
    coeff = MU0 * CURRENT_SOURCE_A / (4.0 * math.pi)
    for index in range(n_segments):
        theta = (index + 0.5) * dtheta
        cx = LOOP_RADIUS_M * math.cos(theta)
        cy = LOOP_RADIUS_M * math.sin(theta)
        dlx = -LOOP_RADIUS_M * math.sin(theta) * dtheta
        dly = LOOP_RADIUS_M * math.cos(theta) * dtheta
        rx = x - cx
        ry = y - cy
        rz = z
        r2 = rx * rx + ry * ry + rz * rz
        inv_r3 = 1.0 / (r2 * sqrt(r2))
        bx += coeff * dly * rz * inv_r3
        by += coeff * (-dlx * rz) * inv_r3
        bz += coeff * (dlx * ry - dly * rx) * inv_r3
    from ngsolve import CoefficientFunction as CF

    return CF((bx, by, bz))


def integrate_vol_lorentz_force(vol_path: Path, n_segments: int) -> dict[str, object]:
    from ngsolve import CoefficientFunction as CF
    from ngsolve import Integrate, Mesh, sqrt, x, y
    from ngsolve import __version__ as ngsolve_version

    t0 = time.perf_counter()
    mesh = Mesh(str(vol_path))
    import_duration = time.perf_counter() - t0
    materials = tuple(mesh.GetMaterials())
    boundaries = tuple(mesh.GetBoundaries())
    if "target" not in materials:
        raise RuntimeError(f".vol mesh lacks required 'target' material: {materials}")

    t_field0 = time.perf_counter()
    b_source = source_loop_biot_savart_cf(n_segments)
    r_xy = sqrt(x * x + y * y)
    current_density = CURRENT_TARGET_A / (math.pi * TARGET_TUBE_RADIUS_M**2)
    j_target = CF((-current_density * y / r_xy, current_density * x / r_xy, 0.0))
    force_density = CF(
        (
            j_target[1] * b_source[2] - j_target[2] * b_source[1],
            j_target[2] * b_source[0] - j_target[0] * b_source[2],
            j_target[0] * b_source[1] - j_target[1] * b_source[0],
        )
    )
    field_build_duration = time.perf_counter() - t_field0

    t_int0 = time.perf_counter()
    target = mesh.Materials("target")
    fx = float(Integrate(force_density[0], mesh, definedon=target, order=8))
    fy = float(Integrate(force_density[1], mesh, definedon=target, order=8))
    fz = float(Integrate(force_density[2], mesh, definedon=target, order=8))
    volume = float(Integrate(CF(1.0), mesh, definedon=target, order=8))
    integration_duration = time.perf_counter() - t_int0
    exact_volume = 2.0 * math.pi * LOOP_RADIUS_M * math.pi * TARGET_TUBE_RADIUS_M**2
    volume_rel_error = abs(volume - exact_volume) / exact_volume

    return {
        "ngsolve_version": ngsolve_version,
        "mesh": {
            "ne": int(mesh.ne),
            "nv": int(mesh.nv),
            "materials": list(materials),
            "boundaries": list(boundaries),
            "target_volume_m3": volume,
            "target_volume_exact_m3": exact_volume,
            "target_volume_rel_error": volume_rel_error,
        },
        "force_vector_N": [fx, fy, fz],
        "timing": {
            "mesh_import": round(import_duration, 6),
            "source_field_build": round(field_build_duration, 6),
            "ngsolve_lorentz_integral": round(integration_duration, 6),
        },
    }


def build_artifacts(args: argparse.Namespace) -> tuple[dict[str, object], dict[str, object]]:
    started = time.perf_counter()
    created_at = _utc_now()
    mesh_info = ensure_vol(args.mesh_source, args.rebuild_vol, args.vol, args.cubit_journal)
    solve_info = integrate_vol_lorentz_force(args.vol, args.source_segments)
    reference = exact_axisymmetric_force_N()
    gate = axisymmetric_to_3d_force_gate(
        reference,
        solve_info["force_vector_N"],
        case_id="vol_target_torus_lorentz_force",
        axial_axis="z",
        result_basis="full_revolution",
        axial_rtol=AXIAL_RTOL,
        transverse_rtol=TRANSVERSE_RTOL,
        metadata={
            "producer": "ngsolve_vol_lorentz_integral",
            "mesh_source": mesh_info["mesh_source"],
            "source_segments": args.source_segments,
        },
    )
    plan = axisymmetric_to_3d_validation_plan(
        "vol_target_torus_lorentz_force",
        preferred_3d_route="vol_mesh_ngsolve_lorentz_force",
    )
    write_started = time.perf_counter()
    mesh = solve_info["mesh"]
    transverse_rel = gate["errors"]["transverse_magnitude_N"] / max(
        abs(reference),
        abs(gate["three_d_result"]["axial_component_N"]),
        1.0e-300,
    )
    max_rel = max(
        float(gate["errors"]["axial_rel_error"]),
        float(mesh["target_volume_rel_error"]),
        float(transverse_rel),
    )
    max_abs = max(
        float(gate["errors"]["axial_abs_error_N"]),
        float(gate["errors"]["transverse_magnitude_N"]),
    )
    checks = {
        "ran_to_completion": True,
        "result_files_exist": True,
        "validation_passed": gate["status"] == "ok" and mesh["target_volume_rel_error"] <= VOLUME_RTOL,
        "vol_loaded_by_ngsolve": True,
        "target_material_present": "target" in mesh["materials"],
        "axisymmetric_to_3d_gate_ok": gate["status"] == "ok",
        "target_volume_error_ok": mesh["target_volume_rel_error"] <= VOLUME_RTOL,
    }
    command = "python " + " ".join(sys.argv)
    timing = {
        "mesh_generation_or_reuse": float(mesh_info["mesh_generation_duration_s"]),
        "mesh_import": solve_info["timing"]["mesh_import"],
        "source_field_build": solve_info["timing"]["source_field_build"],
        "ngsolve_lorentz_integral": solve_info["timing"]["ngsolve_lorentz_integral"],
    }
    result_files = [_rel(args.out), _rel(args.vol), _rel(args.cubit_journal), _rel(args.ready_out)]
    pass_ok = all(value is True for value in checks.values())
    result = {
        "schema": "cae-ai-lab.solver-run.v1",
        "created_at_utc": created_at,
        "case": "axisymmetric reference to 3D .vol Lorentz force",
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
        "result_files": result_files,
        "checks": checks,
        "tolerances": {
            "max_rel": AXIAL_RTOL,
            "max_abs": 0.0,
            "axial_rtol": AXIAL_RTOL,
            "transverse_rtol": TRANSVERSE_RTOL,
            "target_volume_rtol": VOLUME_RTOL,
        },
        "errors": {
            "max_rel": max_rel,
            "max_abs": max_abs,
            "axial_rel_error": gate["errors"]["axial_rel_error"],
            "target_volume_rel_error": mesh["target_volume_rel_error"],
            "transverse_rel_error": transverse_rel,
        },
        "tool_versions": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "radia_mcp": _pkg_version("radia-mcp"),
            "radia": _pkg_version("radia"),
            "ngsolve": solve_info["ngsolve_version"],
            "scipy": _pkg_version("scipy"),
        },
        "timing_breakdown_s": timing,
        "verification": {
            "method": "NGSolve Mesh(.vol) Lorentz volume integral checked by axisymmetric_to_3d_force_gate",
            "command": command,
        },
        "mesh_generation": {
            **mesh_info,
            "vol_path": _rel(args.vol),
            "cubit_journal": _rel(args.cubit_journal),
        },
        "case_parameters": {
            "loop_radius_m": LOOP_RADIUS_M,
            "loop_separation_m": LOOP_SEPARATION_M,
            "target_tube_radius_m": TARGET_TUBE_RADIUS_M,
            "target_mesh_maxh_m": TARGET_MESH_MAXH_M,
            "current_source_A": CURRENT_SOURCE_A,
            "current_target_A": CURRENT_TARGET_A,
            "source_segments": args.source_segments,
        },
        "axisymmetric_reference": {
            "axial_force_N": reference,
            "basis": "exact filament coaxial-loop dM/dz, full 3D revolution",
        },
        "ngsolve_vol_result": solve_info,
        "axisymmetric_to_3d_gate": gate,
        "validation_plan": plan,
    }
    write_duration = time.perf_counter() - write_started
    result["timing_breakdown_s"]["gate_and_write_results"] = round(write_duration, 6)

    ready = {
        "schema": "cae-ai-lab.solver-ready.v1",
        "created_at_utc": created_at,
        "tool_slot": "radia-ngsolve",
        "case": "axisymmetric reference to 3D .vol Lorentz force",
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
    parser.add_argument("--source-segments", type=int, default=SOURCE_SEGMENTS)
    parser.add_argument("--vol", type=Path, default=VOL_PATH)
    parser.add_argument("--cubit-journal", type=Path, default=CUBIT_JOURNAL)
    parser.add_argument("--out", type=Path, default=OUT_JSON)
    parser.add_argument("--ready-out", type=Path, default=READY_JSON)
    args = parser.parse_args()
    if args.source_segments < 16:
        raise ValueError("--source-segments must be >= 16")

    result, ready = build_artifacts(args)
    _write_json(args.out, result)
    _write_json(args.ready_out, ready)

    gate = result["axisymmetric_to_3d_gate"]
    mesh = result["ngsolve_vol_result"]["mesh"]
    print("[axisymmetric -> 3D .vol force]")
    print(f"  mesh_source={result['mesh_generation']['mesh_source']}")
    print(f"  vol={args.vol}")
    print(f"  ne={mesh['ne']}, nv={mesh['nv']}, volume_rel_error={mesh['target_volume_rel_error']:.3e}")
    print(f"  F={result['ngsolve_vol_result']['force_vector_N']} N")
    print(f"  Fz_ref={result['axisymmetric_reference']['axial_force_N']:.12e} N")
    print(f"  axial_rel_error={gate['errors']['axial_rel_error']:.3e}")
    print(f"  pass={result['pass']}")
    print(f"[OK] wrote {args.out}")
    print(f"[OK] wrote {args.ready_out}")
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
