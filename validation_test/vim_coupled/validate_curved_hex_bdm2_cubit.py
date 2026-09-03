"""Validate true Cubit order-2 curved HEX BDM2 on hibino first or idle-CI mdx.

Generate the temporary mesh on a Cubit host first::

    python validation_test/cubit/build_curved_hex_bdm2_cylinder.py

Then copy the ``.vol`` and sibling ``.vol.json`` to the compute host and run
this script there.  This is a heavy validation lane, not a CI test.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from datetime import UTC, datetime
from pathlib import Path

import ngsolve as ng
import numpy as np
from cubit_mesh_export.check import check_consistency

import radia
from radia.vim import FieldFromSolution, MagnetizationSource, Solve
from radia.vim._vim import _hex_mapping_affinity_report

DEFAULT_MESH = Path(r"C:\temp\radia_hex_bdm2\curved_cylinder_hex_q2.vol")
DEFAULT_OUTPUT = Path(__file__).with_name("curved_hex_bdm2_cubit_summary.json")
FIELD_POINTS = np.asarray(
    ((0.05, 0.0, 0.0), (0.0, 0.0, 0.06), (0.04, 0.03, 0.02)),
    dtype=float,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _solve_summary(result: dict) -> dict:
    return {
        "converged": bool(result["last_solve_converged"]),
        "iterations": int(result["iters"]),
        "final_relative_residual": float(result["last_solve_final_relative_residual"]),
        "linear_solver": str(result["linear_solver"]),
        "M_avg": np.asarray(result["M_avg"], dtype=float).tolist(),
        "timing": dict(result.get("timing", {})),
    }


def _source_reference(mesh: ng.Mesh) -> dict:
    source = MagnetizationSource(mesh, (1.0e5, 2.0e5, 3.0e5), order=2, curve_order=2)
    actual = np.asarray(source.Field(FIELD_POINTS, algorithm="direct"), dtype=float)
    position = ng.CF((ng.x, ng.y, ng.z))
    sigma = ng.InnerProduct(source.magnetization, ng.specialcf.normal(3))
    reference = []
    for point in FIELD_POINTS:
        delta = ng.CF(tuple(point)) - position
        radius = ng.sqrt(ng.InnerProduct(delta, delta))
        reference.append(
            ng.Integrate(
                sigma * delta / (4.0 * np.pi * radius**3),
                mesh,
                ng.BND,
                order=14,
            )
        )
    reference_array = np.asarray(reference, dtype=float)
    relative = np.linalg.norm(actual - reference_array, axis=1) / np.maximum(
        np.linalg.norm(reference_array, axis=1), np.finfo(float).tiny
    )
    return {
        "actual": actual.tolist(),
        "ngsolve_boundary_reference": reference_array.tolist(),
        "relative_error": relative.tolist(),
        "maximum_relative_error": float(relative.max()),
        "source_stats": source.stats,
    }


def run(mesh_path: Path) -> dict:
    started = time.perf_counter()
    mesh_path = mesh_path.resolve()
    if not mesh_path.is_file():
        raise FileNotFoundError(mesh_path)
    mesh_check = check_consistency(
        mesh_path, min_curve_order=2, required_materials=("iron",)
    )
    if not mesh_check["passed"]:
        raise RuntimeError(f"check-vol failed: {mesh_check}")
    mesh = ng.Mesh(str(mesh_path))
    affinity = _hex_mapping_affinity_report(mesh)
    if {len(element.vertices) for element in mesh.Elements(ng.VOL)} != {8}:
        raise RuntimeError("curved BDM2 production validation requires pure HEX")

    solve_options = {
        "H_ext": ng.CF((0.0, 0.0, 1000.0)),
        "order": 2,
        "curve_order": 2,
        "gram_eps": 1.0e-14,
        "leaf": 4096,
        "tol": 1.0e-9,
        "maxit": 2000,
    }
    mu0 = 4.0e-7 * np.pi
    linear_bh = np.asarray(
        ((0.0, 0.0), (1.0e3, mu0 * 100.0e3), (1.0e5, mu0 * 100.0e5)),
        dtype=float,
    )
    with ng.TaskManager():
        linear = Solve(mesh, mu_r=100.0, **solve_options)
        linear_field = np.asarray(
            FieldFromSolution(linear, FIELD_POINTS, algorithm="direct"), dtype=float
        )
        nonlinear = Solve(mesh, bh_table=linear_bh, nl_maxit=30, **solve_options)
        source = _source_reference(mesh)

    derivative_rejected = False
    derivative_message = ""
    gram = linear["_charge_gram"]
    stored = gram.hex_stored_nodes()
    cell_nodes = np.asarray(stored["cell_nodes"]).reshape(-1, 27, 3)
    face_nodes = np.asarray(stored["face_nodes"]).reshape(-1, 9, 3)
    try:
        gram.hex_charge_gram_directional_derivative(
            np.zeros_like(cell_nodes), np.zeros_like(face_nodes)
        )
    except RuntimeError as exc:
        derivative_message = str(exc)
        derivative_rejected = "mapped HEX BDM2 shape derivatives" in derivative_message

    linear_m = np.asarray(linear["M_avg"], dtype=float)
    nonlinear_m = np.asarray(nonlinear["M_avg"], dtype=float)
    material_difference = float(
        np.linalg.norm(linear_m - nonlinear_m)
        / max(np.linalg.norm(linear_m), np.finfo(float).tiny)
    )
    sidecar = Path(str(mesh_path) + ".json")
    checks = {
        "cubit_vol_check_passed": bool(mesh_check["passed"]),
        "true_curve_order_two": int(mesh.GetCurveOrder()) == 2,
        "all_hex_cells_are_nonaffine": (
            affinity["cell_count"] == 4 and affinity["nonaffine_cell_count"] == 4
        ),
        "linear_material_solve_converged": bool(linear["last_solve_converged"]),
        "nonlinear_energy_newton_converged": (
            bool(nonlinear["last_solve_converged"])
            and nonlinear["linear_solver"] == "energy-newton-cpp"
        ),
        "linear_and_equivalent_nonlinear_material_match": material_difference <= 1.0e-6,
        "radfld_is_finite": bool(np.isfinite(linear_field).all()),
        "prescribed_field_matches_ngsolve_boundary_integral": (
            source["maximum_relative_error"] <= 1.0e-8
        ),
        "mapped_shape_derivative_fails_loud": derivative_rejected,
    }
    return {
        "schema": "radia.validation.curved-hex-bdm2-cubit.v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "machine": platform.node(),
        "versions": {
            "radia": getattr(radia, "__version__", "unknown"),
            "ngsolve": getattr(ng, "__version__", "unknown"),
            "python": platform.python_version(),
            "cubit_mesh_export": mesh_check["versions"]["cubit_mesh_export"],
        },
        "mesh": {
            "path": str(mesh_path),
            "sha256": _sha256(mesh_path),
            "sidecar_sha256": _sha256(sidecar) if sidecar.is_file() else None,
            "elements": int(mesh.ne),
            "curve_order": int(mesh.GetCurveOrder()),
            "affinity": affinity,
            "minimum_scaled_jacobian": float(
                mesh_check["quality"]["minimum_scaled_jacobian"]
            ),
            "volume_error_percent": float(mesh_check["materials"][0]["error_pct"]),
        },
        "linear": _solve_summary(linear),
        "nonlinear": _solve_summary(nonlinear),
        "equivalent_material_relative_difference": material_difference,
        "field_points": FIELD_POINTS.tolist(),
        "linear_field": linear_field.tolist(),
        "prescribed_source": source,
        "shape_derivative": {
            "rejected": derivative_rejected,
            "message": derivative_message,
        },
        "checks": checks,
        "pass": all(checks.values()),
        "timing_s": {"total": time.perf_counter() - started},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, default=DEFAULT_MESH)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = run(args.mesh)
    args.output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "pass": summary["pass"],
                "checks": summary["checks"],
                "timing_s": summary["timing_s"],
            },
            indent=2,
        ),
        flush=True,
    )
    if not summary["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
