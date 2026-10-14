"""Isolate the fixed-source projection from the HDiv ChargeGram setup cost.

This diagnostic is deliberately a linear zero-load solve on the response iron
mesh.  It preserves the exact HDiv order, curve order, Gram tolerance, and
TaskManager setting of the hybrid-undulator three-engine validation while
removing the prescribed PM field.  The resulting timings identify whether a
slow three-engine run is in ChargeGram construction or source projection.
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent


def _configure_source_overlay() -> Path | None:
    """Use an opt-in Python overlay while retaining the installed native wheel."""
    raw = os.environ.get("RADIA_SOURCE_OVERLAY")
    if raw:
        overlay = Path(raw).resolve()
    else:
        overlay = HERE.parents[1] / "src"
    package = overlay / "radia"
    if not (package / "__init__.py").is_file():
        return None
    installed = importlib.import_module("radia")
    if str(package) not in installed.__path__:
        installed.__path__.insert(0, str(package))
    return package


SOURCE_OVERLAY = _configure_source_overlay()

import ngsolve as ng
import numpy as np
from radia import vim
from radia.esrf_examples import build_esrf_fixed_magnetization_source


def _source_accuracy_points() -> np.ndarray:
    centreline = np.asarray(
        [[0.0, float(y), 0.0] for y in np.linspace(-0.065, 0.065, 27)],
        dtype=float,
    )
    off_axis = np.asarray(
        [
            [x, y, z]
            for x in (-0.015, 0.015)
            for y in np.linspace(-0.060, 0.060, 13)
            for z in (-0.003, 0.003)
        ],
        dtype=float,
    )
    return np.vstack((centreline, off_axis))


def _relative_rms(reference: np.ndarray, candidate: np.ndarray) -> float:
    scale = float(np.sqrt(np.mean(np.sum(reference * reference, axis=1))))
    if scale <= 0.0:
        raise RuntimeError("fixed-source direct-field reference is zero")
    return float(np.sqrt(np.mean(np.sum((candidate - reference) ** 2, axis=1))) / scale)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iron-mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--phase",
        choices=("zero-load-solve", "source-projection", "source-field-compare"),
        default="zero-load-solve",
    )
    parser.add_argument("--source-mesh", type=Path)
    parser.add_argument("--reference-source-mesh", type=Path)
    parser.add_argument("--source-tree-theta", type=float, default=0.14)
    parser.add_argument(
        "--source-field-algorithm",
        choices=("auto", "direct", "tree"),
        default="auto",
    )
    parser.add_argument("--order", choices=(1, 2), type=int, default=2)
    parser.add_argument("--curve-order", choices=(0, 2), type=int, default=2)
    parser.add_argument("--gram-eps", type=float, default=1.0e-12)
    parser.add_argument(
        "--gram-backend",
        choices=("hmat", "exact-dense"),
        default="hmat",
        help="explicit ChargeGram backend; exact-dense requires a memory cap",
    )
    parser.add_argument(
        "--exact-dense-memory-mb",
        type=int,
        default=None,
        help="mandatory memory cap when --gram-backend=exact-dense",
    )
    parser.add_argument("--threads", type=int, default=0)
    options = parser.parse_args(argv)
    if not 0.0 < options.gram_eps < 1.0:
        raise ValueError("--gram-eps must lie in (0, 1)")
    if options.gram_backend == "exact-dense":
        if options.exact_dense_memory_mb is None or options.exact_dense_memory_mb <= 0:
            raise ValueError(
                "--exact-dense-memory-mb must be positive with "
                "--gram-backend=exact-dense"
            )
    elif options.exact_dense_memory_mb is not None:
        raise ValueError(
            "--exact-dense-memory-mb requires --gram-backend=exact-dense"
        )
    if options.threads > 0:
        ng.SetNumThreads(options.threads)
    mesh = ng.Mesh(str(options.iron_mesh.resolve()))
    if options.phase == "source-field-compare":
        if options.source_mesh is None or options.reference_source_mesh is None:
            raise ValueError(
                "--source-mesh and --reference-source-mesh are required for "
                "--phase source-field-compare"
            )
        candidate_mesh = ng.Mesh(str(options.source_mesh.resolve()))
        reference_mesh = ng.Mesh(str(options.reference_source_mesh.resolve()))
        candidate = build_esrf_fixed_magnetization_source(
            candidate_mesh,
            3,
            order=options.order,
            curve_order=(None if options.curve_order == 0 else options.curve_order),
        )
        reference = build_esrf_fixed_magnetization_source(
            reference_mesh,
            3,
            order=options.order,
            curve_order=(None if options.curve_order == 0 else options.curve_order),
        )
        points = _source_accuracy_points()
        candidate_field = np.asarray(candidate.Field(points, algorithm="direct"), dtype=float)
        reference_field = np.asarray(reference.Field(points, algorithm="direct"), dtype=float)
        payload = {
            "schema": "radia.validation.esrf-hybrid-undulator-fixed-source-field.v1",
            "case": "ESRF Example #3 fixed-magnetization source mesh comparison",
            "source_present": True,
            "candidate_source_mesh": str(options.source_mesh.resolve()),
            "reference_source_mesh": str(options.reference_source_mesh.resolve()),
            "order": options.order,
            "curve_order": options.curve_order,
            "sample_count": int(points.shape[0]),
            "relative_rms": _relative_rms(reference_field, candidate_field),
            "maximum_absolute_difference_A_per_m": float(
                np.max(np.linalg.norm(candidate_field - reference_field, axis=1))
            ),
            "candidate_source": candidate.stats,
            "reference_source": reference.stats,
        }
        output = options.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"event": "profile_complete", "output": str(output),
                          "relative_rms": payload["relative_rms"]}, sort_keys=True))
        return 0

    if options.phase == "source-projection":
        if options.source_mesh is None:
            raise ValueError("--source-mesh is required for --phase source-projection")
        source_mesh = ng.Mesh(str(options.source_mesh.resolve()))
        source = build_esrf_fixed_magnetization_source(
            source_mesh,
            3,
            order=options.order,
            curve_order=(None if options.curve_order == 0 else options.curve_order),
            field_cf_algorithm=(
                None if options.source_field_algorithm == "auto"
                else options.source_field_algorithm
            ),
            field_tree_options={"theta": options.source_tree_theta},
        )
        fes = ng.HDiv(mesh, order=options.order)
        load = ng.LinearForm(fes)
        load += source.field_cf * fes.TestFunction() * ng.dx
        started = time.perf_counter()
        with ng.TaskManager():
            load.Assemble()
        payload = {
            "schema": "radia.validation.esrf-hybrid-undulator-source-projection.v1",
            "case": "ESRF Example #3 fixed-magnetization source-to-iron weak load",
            "source_present": True,
            "source_overlay": None if SOURCE_OVERLAY is None else str(SOURCE_OVERLAY),
            "order": options.order,
            "curve_order": options.curve_order,
            "source_tree_theta": options.source_tree_theta,
            "source_field_algorithm_requested": options.source_field_algorithm,
            "threads": options.threads,
            "source_projection_wall_s": time.perf_counter() - started,
            "iron_ndof": int(fes.ndof),
            "source_projection": source.stats,
            "load_norm": float(ng.Norm(load.vec)),
        }
        output = options.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"event": "profile_complete", "output": str(output),
                          "source_projection_wall_s": payload["source_projection_wall_s"]},
                         sort_keys=True))
        return 0

    started = time.perf_counter()
    with ng.TaskManager():
        result = vim.Solve(
            mesh,
            mu_r=1000.0,
            H_ext=ng.CoefficientFunction((0.0, 0.0, 0.0)),
            order=options.order,
            curve_order=(None if options.curve_order == 0 else options.curve_order),
            gram_eps=options.gram_eps,
            tol=1.0e-8,
            maxit=4000,
            preconditioner="mass-riesz",
            gram_backend=options.gram_backend,
            exact_dense_memory_mb=options.exact_dense_memory_mb,
        )
    payload = {
        "schema": "radia.validation.esrf-hybrid-undulator-hdiv-setup.v1",
        "case": "ESRF Example #3 response-iron zero-load setup",
        "source_present": False,
        "order": options.order,
        "curve_order": options.curve_order,
        "gram_eps": options.gram_eps,
        "gram_backend": result.get("gram_backend"),
        "exact_dense_normalized_gram": bool(
            result.get("exact_dense_normalized_gram", False)
        ),
        "exact_dense_memory_mb": options.exact_dense_memory_mb,
        "threads": options.threads,
        "elapsed_wall_s": time.perf_counter() - started,
        "timings": {
            name: result.get(name)
            for name in (
                "fes_wall_s",
                "charge_gram_wall_s",
                "charge_basis_wall_s",
                "charge_gram_cpp_wall_s",
                "hex_state_check_wall_s",
                "projection_wall_s",
                "demag_probe_wall_s",
                "setup_wall_s",
                "solve_wall_s",
                "post_wall_s",
                "total_wall_s_internal",
            )
        },
        "ndof": int(result["ndof"]),
        "n_charge": int(result["n_charge"]),
        "iterations": int(result["iters"]),
        "hmat_stats": result.get("hmat_stats"),
    }
    output = options.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"event": "profile_complete", "output": str(output),
                      "elapsed_wall_s": payload["elapsed_wall_s"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
