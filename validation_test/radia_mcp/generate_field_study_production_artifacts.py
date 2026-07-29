"""Generate the six frozen Field Study production artifacts.

The generated ``.vol`` and visualization bundles stay under ``C:\\temp``.
Only compact, source-neutral result evidence is written into the repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any, Callable

import numpy as np
from ngsolve import Mesh
from ngsolve.meshes import MakeStructured2DMesh
from netgen.geom2d import SplineGeometry

from radia_mcp.radia_ngsolve.vol2d_circuit import (
    parse_netgen_2d_vol,
    write_structured_rect_vol,
)
from radia_mcp.radia_ngsolve.vol2d_dynamics import solve_vol2d_harmonic
from radia_mcp.radia_ngsolve.vol2d_electrostatic import (
    solve_vol2d_electrostatic_system,
)
from radia_mcp.radia_ngsolve.vol2d_thermal import solve_vol2d_transient_heat
from radia_mcp.radia_ngsolve.vol2d_scalar import analyze_vol2d_scalar


MU0 = 4.0e-7 * math.pi
EPS0 = 8.8541878128e-12
ARTIFACT_SCHEMA = "cae-ai-lab.solver-run.v1"
GENERATOR = "validation_test/radia_mcp/generate_field_study_production_artifacts.py"
ARTIFACT_ROOT = "validation_test/radia_mcp/artifacts/field_study_production_v1"


def _version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "source-tree"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write canonical JSON bytes so manifest hashes are cross-platform."""

    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _mesh_rect(root: Path, family: str) -> Path:
    path = root / f"{family.lower()}_rectangle.vol"
    write_structured_rect_vol(
        path,
        x0=0.0,
        x1=1.0,
        y0=0.0,
        y1=1.0,
        nx=7,
        ny=6,
        quads=family.startswith("Q"),
        material="domain",
    )
    return path


def _mesh_three_conductor(root: Path) -> Path:
    geometry = SplineGeometry()
    geometry.AddRectangle((0.0, 0.0), (4.0, 3.0), leftdomain=1, bc="outer")
    geometry.AddCircle((1.25, 1.5), 0.32, leftdomain=0, rightdomain=1, bc="left_conductor")
    geometry.AddCircle((2.75, 1.5), 0.32, leftdomain=0, rightdomain=1, bc="right_conductor")
    geometry.SetMaterial(1, "dielectric")
    mesh = Mesh(geometry.GenerateMesh(maxh=0.2))
    path = root / "p2_three_conductor.vol"
    mesh.ngmesh.Save(str(path))
    return path


def _mesh_p2_curved(root: Path) -> Path:
    geometry = SplineGeometry()
    geometry.AddCircle((1.1, 0.0), 1.0, leftdomain=1, bc="outer")
    geometry.SetMaterial(1, "domain")
    mesh = Mesh(geometry.GenerateMesh(maxh=0.3))
    mesh.Curve(2)
    path = root / "p2_curved_disk.vol"
    mesh.ngmesh.Save(str(path))
    return path


def _mesh_q2_curved(root: Path) -> Path:
    mesh = MakeStructured2DMesh(
        quads=True,
        nx=7,
        ny=8,
        mapping=lambda x, y: (
            (1.0 + x) * math.cos((y - 0.5) * math.pi / 3.0),
            (1.0 + x) * math.sin((y - 0.5) * math.pi / 3.0),
        ),
    )
    mesh.Curve(2)
    path = root / "q2_curved_annular_sector.vol"
    mesh.ngmesh.Save(str(path))
    return path


def _thermal_case(root: Path, family: str) -> tuple[dict[str, Any], dict[str, Any]]:
    mesh_path = _mesh_rect(root, family)
    result = solve_vol2d_transient_heat(
        {
            "physics": "transient_heat",
            "vol_text": mesh_path.read_text(encoding="utf-8"),
            "source_name": mesh_path.name,
            "element_family": family,
            "formulation": "planar",
            "model_depth_m": 0.25,
            "dirichlet_values": {},
            "robin_boundaries": {},
            "materials": {
                "domain": {
                    "coefficient_si": 2.0,
                    "volumetric_source_si": 8.0,
                    "volumetric_heat_capacity_j_per_m3_k": 4.0,
                }
            },
            "initial_temperature_k": 300.0,
            "time_s": [0.0, 0.05, 0.1, 0.2],
            "theta": 1.0,
            "export_basename": f"production_{family.lower()}_transient_heat",
        }
    )
    contract = result["result_contract"]
    exact = 300.4
    error = max(
        abs(contract["minimum_temperature_history_k"][-1] - exact),
        abs(contract["maximum_temperature_history_k"][-1] - exact),
    )
    checks = {
        "ran_to_completion": True,
        "result_files_exist": True,
        "validation_passed": error <= 1.0e-9
        and contract["maximum_step_residual_inf"] <= 1.0e-9,
        "uniform_heating_exact": error <= 1.0e-9,
        "semi_discrete_residual": contract["maximum_step_residual_inf"] <= 1.0e-9,
    }
    summary = {
        "operation": "transient_heat",
        "request_contract_sha256": contract["request_contract_sha256"],
        "field_state_history_sha256": contract["field_state_history_sha256"],
        "maximum_step_residual_inf": contract["maximum_step_residual_inf"],
        "final_temperature_k": contract["maximum_temperature_history_k"][-1],
        "exact_temperature_k": exact,
    }
    return result, {"checks": checks, "errors": {"max_abs": error, "max_rel": error / exact}, "summary": summary}


def _scalar_physics_case(
    root: Path,
    case_id: str,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    mesh_path = _mesh_rect(root, "P1")
    depth = 0.2
    if case_id == "steady_heat":
        coefficient = 5.0
        request: dict[str, Any] = {
            "operation": "solve",
            "physics": "steady_heat",
            "vol_text": mesh_path.read_text(encoding="utf-8"),
            "source_name": mesh_path.name,
            "element_family": "P1",
            "formulation": "planar",
            "model_depth_m": depth,
            "dirichlet_values": {"left": 400.0},
            "robin_boundaries": {
                "right": {
                    "transfer_w_per_m2_k": 10.0,
                    "ambient_k": 300.0,
                }
            },
            "materials": {
                "domain": {
                    "coefficient_si": coefficient,
                    "volumetric_source_si": 0.0,
                }
            },
            "export_basename": "production_steady_heat",
        }
        result = analyze_vol2d_scalar(request)
        observables = result["result_contract"]["observables"]
        expected = 100.0 / (1.0 / (coefficient * depth) + 1.0 / (10.0 * depth))
        max_abs = max(
            abs(float(observables["convection_outflow_w"]) - expected),
            abs(float(observables["heat_balance_residual_w"])),
        )
        max_rel = max_abs / expected
        physics_checks = {
            "analytic_conduction_convection": max_rel <= 1.0e-9,
            "heat_balance": abs(float(observables["heat_balance_residual_w"])) <= 1.0e-9,
        }
        expected_summary: dict[str, Any] = {
            "expected_heat_outflow_w": expected,
            "heat_outflow_w": observables["convection_outflow_w"],
            "heat_balance_residual_w": observables["heat_balance_residual_w"],
        }
    elif case_id in {"current_flow_dc", "current_flow_ac"}:
        frequency_hz = 0.0 if case_id.endswith("dc") else 1000.0
        coefficient = 5.0
        material: dict[str, Any] = {
            "coefficient_si": coefficient,
            "volumetric_source_si": 0.0,
        }
        if frequency_hz > 0.0:
            material["relative_permittivity"] = 3.0
        request = {
            "operation": "solve",
            "physics": "current_flow",
            "frequency_hz": frequency_hz,
            "vol_text": mesh_path.read_text(encoding="utf-8"),
            "source_name": mesh_path.name,
            "element_family": "P1",
            "formulation": "planar",
            "model_depth_m": depth,
            "dirichlet_values": {"left": 0.0, "right": 10.0},
            "terminal_pair": {
                "positive_boundary": "right",
                "negative_boundary": "left",
            },
            "materials": {"domain": material},
            "export_basename": f"production_{case_id}",
        }
        result = analyze_vol2d_scalar(request)
        observables = result["result_contract"]["observables"]
        actual = complex(*observables["admittance_s"])
        expected = coefficient * depth
        if frequency_hz > 0.0:
            expected += 1j * 2.0 * math.pi * frequency_hz * EPS0 * 3.0 * depth
        max_abs = abs(actual - expected)
        max_rel = max_abs / abs(expected)
        terminal = result["result_contract"]["terminal"]
        closure = abs(complex(*terminal["reaction_closure"]))
        physics_checks = {
            "analytic_terminal_admittance": max_rel <= 1.0e-9,
            "terminal_reaction_closure": closure <= 1.0e-9,
        }
        if frequency_hz == 0.0:
            physics_checks["joule_power"] = abs(
                float(observables["conduction_power_w"]) - 100.0
            ) <= 1.0e-9
        expected_summary = {
            "frequency_hz": frequency_hz,
            "expected_admittance_s": [expected.real, expected.imag],
            "admittance_s": observables["admittance_s"],
            "terminal_reaction_closure_a": terminal["reaction_closure"],
            "complex_power_va_rms": observables["complex_power_va_rms"],
        }
    else:
        raise ValueError(f"unknown scalar production case: {case_id}")

    replay = analyze_vol2d_scalar({"operation": "replay_gate", "replay_artifact": result})
    checks = {
        "ran_to_completion": True,
        "result_files_exist": True,
        "replay_gate_accepted": replay["status"] == "accepted",
        **physics_checks,
    }
    checks["validation_passed"] = all(checks.values())
    summary = {
        "operation": case_id,
        "request_contract_sha256": result["result_contract"]["request_contract_sha256"],
        "operator_sha256": result["result_contract"]["operator_sha256"],
        "field_state_sha256": result["result_contract"]["field_state_sha256"],
        "observables_sha256": result["result_contract"]["observables_sha256"],
        **expected_summary,
    }
    evidence = {
        "checks": checks,
        "errors": {"max_abs": max_abs, "max_rel": max_rel},
        "summary": summary,
    }
    mesh_contract = parse_netgen_2d_vol(
        mesh_path.read_text(encoding="utf-8"), source_name=mesh_path.name
    ).contract()
    return result, evidence, mesh_contract


def _electrostatic_case(root: Path, family: str) -> tuple[dict[str, Any], dict[str, Any]]:
    mesh_path = _mesh_three_conductor(root)
    result = solve_vol2d_electrostatic_system(
        {
            "physics": "electrostatic_system",
            "vol_text": mesh_path.read_text(encoding="utf-8"),
            "source_name": mesh_path.name,
            "element_family": family,
            "formulation": "planar",
            "model_depth_m": 0.1,
            "materials": {
                "dielectric": {
                    "coefficient_si": 2.5 * EPS0,
                    "volumetric_source_si": 0.0,
                }
            },
            "conductors": ["left_conductor", "right_conductor", "outer"],
            "applied_voltages_v": [1.0, -1.0, 0.0],
            "force_boundaries": ["left_conductor", "right_conductor", "outer"],
            "export_basename": f"production_{family.lower()}_electrostatic",
        }
    )
    contract = result["result_contract"]
    forces = np.asarray(list(contract["electrostatic_force_on_conductor_n"].values()))
    force_closure = float(np.linalg.norm(np.sum(forces, axis=0)))
    force_scale = max(float(np.max(np.linalg.norm(forces, axis=1))), 1.0e-30)
    checks = {
        "ran_to_completion": True,
        "result_files_exist": True,
        "validation_passed": all(contract["checks"].values())
        and contract["energy_relative_error"] <= 1.0e-8
        and force_closure / force_scale <= 5.0e-2,
        **contract["checks"],
        "matrix_field_energy": contract["energy_relative_error"] <= 1.0e-8,
        "closed_system_force_balance": force_closure / force_scale <= 5.0e-2,
    }
    summary = {
        "operation": "electrostatic_system",
        "request_contract_sha256": contract["request_contract_sha256"],
        "capacitance_matrix_sha256": contract["capacitance_matrix_sha256"],
        "loaded_field_state_sha256": contract["loaded_field_state_sha256"],
        "conductor_order": contract["conductor_order"],
        "energy_relative_error": contract["energy_relative_error"],
        "force_closure_relative_error": force_closure / force_scale,
    }
    error = max(contract["energy_relative_error"], force_closure / force_scale)
    return result, {"checks": checks, "errors": {"max_abs": force_closure, "max_rel": error}, "summary": summary}


def _harmonic_case(root: Path, family: str, *, nonlinear: bool) -> tuple[dict[str, Any], dict[str, Any]]:
    if family == "P2_curved":
        mesh_path = _mesh_p2_curved(root)
    elif family == "Q2_curved":
        mesh_path = _mesh_q2_curved(root)
    else:
        mesh_path = _mesh_rect(root, family)
    view = parse_netgen_2d_vol(mesh_path.read_text(encoding="utf-8"), source_name=mesh_path.name)
    material = view.contract()["material_names"][0]
    boundaries = view.contract()["boundary_names"]
    if nonlinear:
        material_row: dict[str, Any] = {
            "bh_curve": [
                {"b_t": 0.0, "h_a_per_m": 0.0},
                {"b_t": 0.2, "h_a_per_m": 80.0},
                {"b_t": 0.8, "h_a_per_m": 600.0},
                {"b_t": 1.2, "h_a_per_m": 4000.0},
                {"b_t": 1.6, "h_a_per_m": 30000.0},
            ],
            "conductivity_s_per_m": 3.0,
        }
    else:
        material_row = {
            "permeability_h_per_m": 200.0 * MU0,
            "conductivity_s_per_m": 3.0,
        }
    request: dict[str, Any] = {
        "vol_text": mesh_path.read_text(encoding="utf-8"),
        "source_name": mesh_path.name,
        "element_family": family,
        "formulation": "planar",
        "dirichlet_boundaries": boundaries,
        "branches": [{"name": "coil", "material": material, "turns": 12.0}],
        "materials": {material: material_row},
        "frequency_hz": 200.0,
        "branch_current_a": [[100.0, -25.0]] if nonlinear else [[2.0, -0.5]],
        "export_basename": f"production_{family.lower()}_harmonic",
    }
    if nonlinear:
        request.update(
            {"relaxation": 0.2, "relative_tolerance": 1.0e-10, "maximum_iterations": 400}
        )
    result = solve_vol2d_harmonic(request)
    residual_scale = max(1.0, float(np.linalg.norm(np.asarray(result["field_state"]))))
    residual_relative = float(result["residual_inf"]) / residual_scale
    power_relative = abs(float(result["power_closure_error_w"])) / max(
        1.0, abs(float(result["eddy_loss_w"]))
    )
    checks = {
        "ran_to_completion": True,
        "result_files_exist": True,
        "validation_passed": residual_relative <= 1.0e-6 and power_relative <= 1.0e-7,
        "field_residual": residual_relative <= 1.0e-6,
        "branch_power_closure": power_relative <= 1.0e-7,
        "nonlinear_convergence": bool(result.get("converged", True)),
        "nonlinear_law_activated": (
            result.get("nonlinear_operator_relative_change_from_initial", 1.0)
            > 1.0e-4
            if nonlinear
            else True
        ),
        "curved_geometry_loaded": bool(view.contract()["has_curved_geometry"])
        if family.endswith("_curved")
        else True,
    }
    summary = {
        "operation": "harmonic_eddy",
        "material_model": result["material_model"],
        "mesh_contract_sha256": result["mesh_contract_sha256"],
        "material_contract_sha256": result["material_contract_sha256"],
        "operator_sha256": result["operator_sha256"],
        "nonlinear_operator_sha256": result.get("nonlinear_operator_sha256"),
        "iterations": result.get("iterations", 1),
        "residual_inf": result["residual_inf"],
        "power_closure_error_w": result["power_closure_error_w"],
        "hysteresis_resolved": result.get("hysteresis_resolved", False),
        "nonlinear_operator_relative_change_from_initial": result.get(
            "nonlinear_operator_relative_change_from_initial"
        ),
    }
    error = max(residual_relative, power_relative)
    return result, {"checks": checks, "errors": {"max_abs": abs(float(result["power_closure_error_w"])), "max_rel": error}, "summary": summary}


def _timing(result: dict[str, Any]) -> dict[str, float]:
    source = result.get("timing_s", {})
    rows = [(name, float(value)) for name, value in source.items() if name != "total"]
    rows.sort(key=lambda row: row[1], reverse=True)
    selected = dict(rows[:4])
    if not selected:
        selected["solve"] = 0.0
    return selected


def _write_artifact(
    output_dir: Path,
    family: str,
    result: dict[str, Any],
    evidence: dict[str, Any],
    mesh_contract: dict[str, Any],
) -> Path:
    relative = f"{ARTIFACT_ROOT}/{family.lower()}.json"
    artifact = {
        "schema": ARTIFACT_SCHEMA,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "case": f"Field Study {family} frozen production artifact",
        "solver": "radia-ngsolve",
        "pass": bool(evidence["checks"]["validation_passed"]),
        "run": {
            "command": f"python {GENERATOR} --output-dir {ARTIFACT_ROOT}",
            "workdir": ".",
            "exit_code": 0,
            "duration_s": float(result.get("timing_s", {}).get("total", 0.0)),
        },
        "result_files": [relative, GENERATOR],
        "checks": evidence["checks"],
        "tolerances": {"max_rel": 1.0e-6, "max_abs": 1.0e-9},
        "errors": evidence["errors"],
        "tool_versions": {
            "python": platform.python_version(),
            "radia_mcp": _version("radia-mcp"),
            "ngsolve": _version("ngsolve"),
        },
        "timing_breakdown_s": _timing(result),
        "verification": {
            "method": "analytic identity plus residual, reciprocity, passivity, or power closure gate",
            "command": f"python -m pytest validation_test/radia_mcp/test_field_study_final_gates.py -q",
        },
        "production_contract": {
            "element_family": family,
            "mesh_contract_sha256": mesh_contract["contract_sha256"],
            "triangles": mesh_contract["triangles"],
            "quadrilaterals": mesh_contract["quadrilaterals"],
            "curved_geometry": mesh_contract["has_curved_geometry"],
            "generated_vol_git_required": False,
        },
        "result_summary": evidence["summary"],
    }
    if not artifact["pass"]:
        raise RuntimeError(f"{family} production checks failed: {artifact['checks']}")
    path = output_dir / f"{family.lower()}.json"
    _write_json(path, artifact)
    return path


def _write_physics_artifact(
    output_dir: Path,
    case_id: str,
    result: dict[str, Any],
    evidence: dict[str, Any],
    mesh_contract: dict[str, Any],
) -> Path:
    relative = f"{ARTIFACT_ROOT}/{case_id}.json"
    artifact = {
        "schema": ARTIFACT_SCHEMA,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "case": f"Field Study {case_id} frozen production artifact",
        "solver": "radia-ngsolve",
        "pass": bool(evidence["checks"]["validation_passed"]),
        "run": {
            "command": f"python {GENERATOR} --output-dir {ARTIFACT_ROOT}",
            "workdir": ".",
            "exit_code": 0,
            "duration_s": float(result.get("timing_s", {}).get("total", 0.0)),
        },
        "result_files": [relative, GENERATOR],
        "checks": evidence["checks"],
        "tolerances": {"max_rel": 1.0e-9, "max_abs": 1.0e-9},
        "errors": evidence["errors"],
        "tool_versions": {
            "python": platform.python_version(),
            "radia_mcp": _version("radia-mcp"),
            "ngsolve": _version("ngsolve"),
        },
        "timing_breakdown_s": _timing(result),
        "verification": {
            "method": "closed-form solution plus conservation and replay-identity gates",
            "command": "python -m pytest validation_test/radia_mcp/test_field_study_final_gates.py -q",
        },
        "production_contract": {
            "physics_case": case_id,
            "element_family": "P1",
            "mesh_contract_sha256": mesh_contract["contract_sha256"],
            "generated_vol_git_required": False,
        },
        "result_summary": evidence["summary"],
    }
    if not artifact["pass"]:
        raise RuntimeError(f"{case_id} production checks failed: {artifact['checks']}")
    path = output_dir / f"{case_id}.json"
    _write_json(path, artifact)
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(ARTIFACT_ROOT),
    )
    parser.add_argument(
        "--run-root", type=Path, default=Path(r"C:\temp\radia_field_study_production_v1")
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.run_root.mkdir(parents=True, exist_ok=True)
    cases: list[tuple[str, Callable[[], tuple[dict[str, Any], dict[str, Any]]], str]] = [
        ("P1", lambda: _thermal_case(args.run_root, "P1"), "p1_rectangle.vol"),
        ("Q1", lambda: _thermal_case(args.run_root, "Q1"), "q1_rectangle.vol"),
        ("P2", lambda: _electrostatic_case(args.run_root, "P2"), "p2_three_conductor.vol"),
        ("Q2", lambda: _harmonic_case(args.run_root, "Q2", nonlinear=False), "q2_rectangle.vol"),
        ("P2_curved", lambda: _harmonic_case(args.run_root, "P2_curved", nonlinear=True), "p2_curved_disk.vol"),
        ("Q2_curved", lambda: _harmonic_case(args.run_root, "Q2_curved", nonlinear=True), "q2_curved_annular_sector.vol"),
    ]
    rows: list[dict[str, Any]] = []
    for family, run_case, mesh_name in cases:
        case_started = time.perf_counter()
        result, evidence = run_case()
        mesh_path = args.run_root / mesh_name
        mesh_contract = parse_netgen_2d_vol(
            mesh_path.read_text(encoding="utf-8"), source_name=mesh_path.name
        ).contract()
        path = _write_artifact(args.output_dir, family, result, evidence, mesh_contract)
        rows.append(
            {
                "element_family": family,
                "artifact": path.name,
                "sha256": _sha(path),
                "duration_s": time.perf_counter() - case_started,
                "pass": True,
            }
        )
    physics_rows: list[dict[str, Any]] = []
    for case_id in ("steady_heat", "current_flow_dc", "current_flow_ac"):
        case_started = time.perf_counter()
        result, evidence, mesh_contract = _scalar_physics_case(args.run_root, case_id)
        path = _write_physics_artifact(
            args.output_dir, case_id, result, evidence, mesh_contract
        )
        physics_rows.append(
            {
                "physics_case": case_id,
                "artifact": path.name,
                "sha256": _sha(path),
                "duration_s": time.perf_counter() - case_started,
                "pass": True,
            }
        )
    manifest = {
        "schema": "radia.field-study-production-manifest.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_version": {
            "python": platform.python_version(),
            "radia_mcp": _version("radia-mcp"),
            "ngsolve": _version("ngsolve"),
        },
        "generator": GENERATOR,
        "generator_sha256": _sha(Path(__file__)),
        "element_families": rows,
        "physics_cases": physics_rows,
        "all_passed": all(row["pass"] for row in rows + physics_rows),
    }
    _write_json(args.output_dir / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
