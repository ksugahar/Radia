"""Three-way nonlinear validation for ESRF Example #3 hybrid undulator.

The permanent magnets are a prescribed C++ ``MagnetizationSource``.  They are
never material-response unknowns: the HDiv iron response, HCurl reduced-A, and
TOSCA-style mixed total/reduced Omega paths evaluate the same source field.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
MU0 = 4.0e-7 * math.pi
MIXED_DOMAIN_LABEL = "H1 TOSCA mixed total/reduced Omega"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _progress(event: str, **values: object) -> None:
    print(json.dumps({"event": event, **values}, sort_keys=True), flush=True)


def _configure_source_overlay() -> Path | None:
    """Prepend source-only Python modules without copying a native extension."""
    raw = os.environ.get("RADIA_SOURCE_OVERLAY")
    if raw:
        overlay = Path(raw).resolve()
    else:
        # A repository invocation should naturally exercise its checkout.  A
        # compute-node invocation has no source tree here and must opt in via
        # RADIA_SOURCE_OVERLAY rather than silently mixing revisions.
        local_source = HERE.parents[1] / "src"
        if not (local_source / "radia" / "__init__.py").is_file():
            return None
        overlay = local_source.resolve()
    package = overlay / "radia"
    if not (package / "__init__.py").is_file():
        raise FileNotFoundError(
            "RADIA_SOURCE_OVERLAY must contain radia/__init__.py; "
            f"got {overlay}"
        )
    # Load the installed wheel first, including its native extension.  Then make
    # only Python modules resolve from the source overlay.  Copying a .pyd into
    # a validation directory is forbidden because it hides ABI drift.
    installed = importlib.import_module("radia")
    if str(package) not in installed.__path__:
        installed.__path__.insert(0, str(package))
    return package


SOURCE_OVERLAY = _configure_source_overlay()

import ngsolve as ng
import numpy as np
import radia as rad

from radia import vim
from radia.electromagnet_validation import (
    require_static_electromagnet_three_engine_contract,
)
from radia.esrf_examples import (
    build_esrf_fixed_magnetization_source,
    get_esrf_bh_table,
)
from radia.kelvin_identify_ngsolve import detect_kelvin_offset, has_kelvin_identification
from radia.static_electromagnet import (
    StaticElectromagnetMixedDomain,
    solve_static_electromagnet_mixed_total_reduced_omega,
)
from radia.vector_potential_solver import VectorPotentialSolver


MIXED_DOMAIN = StaticElectromagnetMixedDomain(
    reduced_materials=("air",),
    total_materials=("iron", "kelvin"),
    nonlinear_materials=("iron",),
)


def observation_points() -> np.ndarray:
    """Air-gap centreline samples from the original U46 field profile."""
    return np.asarray(
        [[0.0, float(y), 0.0] for y in np.linspace(-0.065, 0.065, 27)],
        dtype=float,
    )


def source_accuracy_points() -> np.ndarray:
    """Centrelines and off-axis points inside the physical undulator gap.

    The native tree is used while NGSolve assembles volume forms, rather than
    only at the centreline probes reported by the original Radia example.  Its
    direct-reference gate must therefore include a small three-dimensional
    air-gap stencil, without sampling a magnet surface where either kernel is
    intentionally singular.
    """
    centreline = observation_points()
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


def _evaluate_cf(field, mesh: ng.Mesh, points: np.ndarray) -> np.ndarray:
    values = []
    for point in points:
        values.append(np.asarray(field(mesh(*map(float, point))), dtype=float))
    return np.asarray(values, dtype=float).reshape(-1, 3)


def _relative_rms(reference: np.ndarray, candidate: np.ndarray) -> float:
    denominator = float(np.sqrt(np.mean(np.sum(reference * reference, axis=1))))
    if denominator <= 0.0:
        raise RuntimeError("hybrid-undulator comparison has zero reference field")
    return float(
        np.sqrt(np.mean(np.sum((candidate - reference) ** 2, axis=1))) / denominator
    )


def _source_tree_accuracy(source, points: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    """Evaluate exact and tree fields, preserving their cost as evidence."""
    direct_started = time.perf_counter()
    direct = np.asarray(source.Field(points, algorithm="direct"), dtype=float)
    direct_wall_s = time.perf_counter() - direct_started
    tree_started = time.perf_counter()
    tree = np.asarray(source.Field(points, algorithm="tree"), dtype=float)
    tree_wall_s = time.perf_counter() - tree_started
    return direct, {
        "tree_vs_direct_relative_rms": _relative_rms(direct, tree),
        "direct_field_wall_s": direct_wall_s,
        "tree_field_wall_s": tree_wall_s,
        "tree_speedup": direct_wall_s / max(tree_wall_s, np.finfo(float).tiny),
        "sample_count": int(points.shape[0]),
    }


def _source_mesh_direct_accuracy(candidate, reference, points: np.ndarray) -> dict[str, float]:
    """Certify a coarsened fixed-M source against its retained direct source.

    A source mesh is independent from the iron-response mesh.  Coarsening it is
    therefore allowed only when the immutable C++ direct fields agree on the
    three-dimensional gap stencil.  This is deliberately not a treecode check:
    it validates the physical source discretization before any acceleration is
    considered.
    """
    reference_started = time.perf_counter()
    reference_field = np.asarray(
        reference.Field(points, algorithm="direct"), dtype=float
    )
    reference_wall_s = time.perf_counter() - reference_started
    candidate_started = time.perf_counter()
    candidate_field = np.asarray(
        candidate.Field(points, algorithm="direct"), dtype=float
    )
    candidate_wall_s = time.perf_counter() - candidate_started
    return {
        "algorithm": "direct",
        "relative_rms": _relative_rms(reference_field, candidate_field),
        "maximum_absolute_difference_A_per_m": float(
            np.max(np.linalg.norm(candidate_field - reference_field, axis=1))
        ),
        "reference_field_wall_s": reference_wall_s,
        "candidate_field_wall_s": candidate_wall_s,
        "sample_count": int(points.shape[0]),
    }


def _pairwise_metrics(fields: dict[str, np.ndarray]) -> dict[str, object]:
    values: dict[str, object] = {}
    names = tuple(fields)
    for index, left in enumerate(names):
        for right in names[index + 1 :]:
            left_field, right_field = fields[left], fields[right]
            values[f"{left}__vs__{right}"] = {
                "relative_rms": _relative_rms(left_field, right_field),
                "maximum_absolute_difference_T": float(
                    np.max(np.linalg.norm(left_field - right_field, axis=1))
                ),
            }
    return values


def _checkpoint_contract(**values: object) -> dict[str, object]:
    return dict(values)


def _read_checkpoint(path: Path, contract: dict[str, object]):
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("contract") != contract:
        raise RuntimeError(f"checkpoint contract changed: remove {path}")
    return np.asarray(payload["field_T"], dtype=float), dict(payload["diagnostics"])


def _write_checkpoint(path: Path, contract: dict[str, object], field: np.ndarray,
                      diagnostics: dict[str, object]) -> None:
    path.write_text(
        json.dumps(
            {
                "schema": "radia.validation.esrf-hybrid-undulator-checkpoint.v1",
                "contract": contract,
                "field_T": np.asarray(field, dtype=float).tolist(),
                "diagnostics": diagnostics,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )


def _source_stats(source) -> dict[str, object]:
    return {
        str(name): value.item() if isinstance(value, np.generic) else value
        for name, value in dict(source.stats).items()
    }


def _solve_hdiv(iron_mesh: ng.Mesh, source, bh_table, points: np.ndarray, *,
                order: int, gram_eps: float, nonlinear_tolerance: float,
                nonlinear_maximum_iterations: int,
                nonlinear_solver: str, gram_backend: str,
                exact_dense_memory_mb: int | None) -> tuple[np.ndarray, dict[str, object]]:
    started = time.perf_counter()
    with ng.TaskManager():
        result = vim.Solve(
            iron_mesh,
            bh_table=bh_table,
            magnetization_sources=source,
            order=int(order),
            curve_order=2,
            gram_eps=float(gram_eps),
            tol=1.0e-8,
            maxit=12000,
            nl_tol=float(nonlinear_tolerance),
            nl_maxit=int(nonlinear_maximum_iterations),
            nonlinear_solver=nonlinear_solver,
            preconditioner="mass-riesz",
            gram_backend=gram_backend,
            exact_dense_memory_mb=exact_dense_memory_mb,
        )
        demag_h = np.asarray(
            vim.FieldFromSolution(result, points, algorithm="direct"), dtype=float
        )
    source_h = np.asarray(source.Field(points, algorithm="direct"), dtype=float)
    nonlinear_stats = dict(result.get("nonlinear_solve_stats", {}))
    nonlinear_stats["converged"] = bool(
        nonlinear_stats.get("nonlinear_converged_final_stage", False)
    )
    return MU0 * (source_h + demag_h), {
        "formulation": "HDiv-MMM",
        "discretization": f"BDM{int(order)} curved-Q2 HEX response",
        "open_boundary": "Coulomb charge Gram; iron-only response mesh",
        "fixed_magnetization_source": "independent C++ MagnetizationSource",
        "mesh_elements": int(iron_mesh.ne),
        "mesh_vertices": int(iron_mesh.nv),
        "ndof": int(result["ndof"]),
        "linear_iterations": int(result["iters"]),
        "nonlinear_stats": nonlinear_stats,
        "preconditioner": result.get("preconditioner"),
        "nonlinear_solver": nonlinear_solver,
        "linear_solver": result.get("linear_solver"),
        "gram_eps": float(gram_eps),
        "gram_backend": result.get("gram_backend"),
        "exact_dense_normalized_gram": bool(
            result.get("exact_dense_normalized_gram", False)),
        "exact_dense_memory_mb": exact_dense_memory_mb,
        "runtime_s": float(time.perf_counter() - started),
    }


def _solve_reduced_a(mesh: ng.Mesh, source, bh_table, points: np.ndarray, *,
                     order: int, kelvin_radius: float, kelvin_center,
                     linear_solver: str, nonlinear_tolerance: float,
                     nonlinear_maximum_iterations: int,
                     nonlinear_relaxation: float) -> tuple[np.ndarray, dict[str, object]]:
    started = time.perf_counter()
    solver = VectorPotentialSolver(
        mesh,
        iron_domains="iron",
        mu_r=1000.0,
        order=int(order),
        kelvin_region="kelvin",
        kelvin_radius=float(kelvin_radius),
        kelvin_center=kelvin_center,
    )
    solver.set_source_cf(MU0 * source.field_cf)
    with ng.TaskManager():
        solution = solver.solve_nonlinear(
            bh_table,
            tol=float(nonlinear_tolerance),
            maxiter=int(nonlinear_maximum_iterations),
            relax=float(nonlinear_relaxation),
            dirichlet="GND",
            verbose=False,
            solver=linear_solver,
        )
    return _evaluate_cf(solver.get_B(), mesh, points), {
        "formulation": "HCurl reduced-A",
        "open_boundary": "periodic spherical Kelvin transform",
        "source_contract": "same C++ fixed-magnetization B source; iron contrast RHS",
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "ndof": int(solution.space.ndof),
        "nonlinear_stats": dict(getattr(solver, "_last_nonlinear_stats", {})),
        "linear_solver": linear_solver,
        "runtime_s": float(time.perf_counter() - started),
    }


def _solve_mixed_omega(mesh: ng.Mesh, source, bh_table, points: np.ndarray, *,
                       order: int, kelvin_radius: float, kelvin_center,
                       nonlinear_tolerance: float,
                       nonlinear_maximum_iterations: int,
                       nonlinear_relaxation: float,
                       source_potential_tolerance: float | None) -> tuple[np.ndarray, dict[str, object]]:
    started = time.perf_counter()
    with ng.TaskManager():
        result = solve_static_electromagnet_mixed_total_reduced_omega(
            mesh,
            source.field_cf,
            MIXED_DOMAIN,
            float(kelvin_radius),
            kelvin_center,
            order=int(order),
            bh_table=bh_table,
            source_potential_contract="global_physical",
            source_trace_tolerance=source_potential_tolerance,
            nonlinear_tolerance=float(nonlinear_tolerance),
            nonlinear_max_iterations=int(nonlinear_maximum_iterations),
            nonlinear_relaxation=float(nonlinear_relaxation),
        )
    return _evaluate_cf(result["B_cf"], mesh, points), {
        "formulation": MIXED_DOMAIN_LABEL,
        "open_boundary": "periodic spherical Kelvin transform",
        "source_contract": "one global physical scalar potential for current-free fixed M",
        "mesh_elements": int(mesh.ne),
        "mesh_vertices": int(mesh.nv),
        "ndof": int(result["fes"].ndof),
        "nonlinear_stats": dict(result.get("nonlinear_stats", {})),
        "source_potential": dict(
            result["static_electromagnet_contract"]["source_trace"]
        ),
        "runtime_s": float(time.perf_counter() - started),
    }


def _verify_overlay() -> dict[str, str]:
    modules = {
        name: importlib.import_module(name)
        for name in (
            "radia.electromagnet_validation",
            "radia.esrf_examples",
            "radia.kelvin_solver",
            "radia.static_electromagnet",
        )
    }
    paths = {name: str(Path(module.__file__).resolve()) for name, module in modules.items()}
    if SOURCE_OVERLAY is not None:
        for name, path in paths.items():
            if not Path(path).is_relative_to(SOURCE_OVERLAY):
                raise RuntimeError(f"{name} did not load from RADIA_SOURCE_OVERLAY: {path}")
    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assets-dir", type=Path, required=True)
    parser.add_argument("--fem-mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--source-mesh",
        type=Path,
        default=None,
        help="independent fixed-M source .vol; defaults to the asset mesh",
    )
    parser.add_argument(
        "--source-reference-mesh",
        type=Path,
        default=None,
        help="retained direct-source .vol required whenever --source-mesh is set",
    )
    parser.add_argument(
        "--source-mesh-relative-rms-tolerance",
        type=float,
        default=2.0e-5,
        help="direct-field source-mesh coarsening acceptance limit",
    )
    parser.add_argument("--hdiv-order", choices=(1, 2), type=int, default=2)
    parser.add_argument("--fem-order", type=int, default=1)
    parser.add_argument("--hdiv-gram-eps", type=float, default=1.0e-12)
    parser.add_argument(
        "--hdiv-gram-backend",
        choices=("hmat", "exact-dense"),
        default="hmat",
        help="explicit HDiv charge-Gram backend; exact-dense requires a memory cap",
    )
    parser.add_argument(
        "--hdiv-exact-dense-memory-mb",
        type=int,
        default=None,
        help="mandatory memory cap when --hdiv-gram-backend=exact-dense",
    )
    parser.add_argument("--nonlinear-tolerance", type=float, default=2.0e-5)
    parser.add_argument("--nonlinear-maximum-iterations", type=int, default=80)
    parser.add_argument("--nonlinear-relaxation", type=float, default=0.2)
    parser.add_argument(
        "--hdiv-nonlinear-solver",
        choices=("picard-mass-riesz", "energy-newton", "picard-energy"),
        default="picard-mass-riesz",
        help="explicit nonlinear HDiv solver; Picard mass-Riesz is the forward validation baseline",
    )
    parser.add_argument("--reduced-a-solver", choices=("direct", "bddc", "ams", "auto"), default="direct")
    parser.add_argument("--source-potential-tolerance", type=float, default=None)
    parser.add_argument(
        "--source-tree-theta",
        type=float,
        default=0.25,
        help="validated native tree opening angle for the fixed PM source",
    )
    parser.add_argument("--source-tree-relative-tolerance", type=float, default=2.0e-5)
    parser.add_argument(
        "--source-field-algorithm",
        choices=("direct", "tree"),
        default="direct",
        help="native field path while assembling each formulation",
    )
    parser.add_argument("--relative-rms-tolerance", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--preflight",
        action="store_true",
        help="validate meshes, source projection, overlay provenance, and gap samples only",
    )
    options = parser.parse_args(argv)
    if options.fem_order < 1:
        raise ValueError("--fem-order must be positive")
    if not 0.0 < options.hdiv_gram_eps < 1.0:
        raise ValueError("--hdiv-gram-eps must lie in (0, 1)")
    if options.hdiv_gram_backend == "exact-dense":
        if (options.hdiv_exact_dense_memory_mb is None
                or options.hdiv_exact_dense_memory_mb <= 0):
            raise ValueError(
                "--hdiv-exact-dense-memory-mb must be positive with "
                "--hdiv-gram-backend=exact-dense")
    elif options.hdiv_exact_dense_memory_mb is not None:
        raise ValueError(
            "--hdiv-exact-dense-memory-mb requires --hdiv-gram-backend=exact-dense")
    if not 0.0 < options.nonlinear_relaxation <= 1.0:
        raise ValueError("--nonlinear-relaxation must lie in (0, 1]")
    if options.source_tree_relative_tolerance <= 0.0:
        raise ValueError("--source-tree-relative-tolerance must be positive")
    if options.source_mesh_relative_rms_tolerance <= 0.0:
        raise ValueError("--source-mesh-relative-rms-tolerance must be positive")
    if options.source_tree_theta <= 0.0:
        raise ValueError("--source-tree-theta must be positive")
    if options.threads > 0:
        ng.SetNumThreads(int(options.threads))

    assets_dir = options.assets_dir.resolve()
    iron_vol = assets_dir / "model.vol"
    default_source_vol = assets_dir / "fixed_magnetization_sources" / "magnet_source.vol"
    source_vol = (
        default_source_vol
        if options.source_mesh is None
        else options.source_mesh.resolve()
    )
    source_mesh_overridden = options.source_mesh is not None
    source_reference_vol = (
        None
        if options.source_reference_mesh is None
        else options.source_reference_mesh.resolve()
    )
    if source_mesh_overridden and source_reference_vol is None:
        raise ValueError(
            "--source-reference-mesh is required with --source-mesh; "
            "a coarsened fixed-M source needs a direct-field certificate"
        )
    fem_mesh_path = options.fem_mesh.resolve()
    for path in (iron_vol, source_vol, fem_mesh_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    fem_report = fem_mesh_path.with_suffix(".mesh.json")
    if not fem_report.is_file():
        raise FileNotFoundError(fem_report)
    fem_contract = json.loads(fem_report.read_text(encoding="utf-8"))
    if fem_contract.get("identification_count", 0) <= 0:
        raise RuntimeError("FEM mesh report has no Kelvin identification")

    implementation_paths = _verify_overlay()
    iron_mesh = ng.Mesh(str(iron_vol))
    source_mesh = ng.Mesh(str(source_vol))
    fem_mesh = ng.Mesh(str(fem_mesh_path))
    if not has_kelvin_identification(fem_mesh):
        raise RuntimeError("FEM mesh has no Kelvin periodic identification")
    kelvin_offset = tuple(float(value) for value in detect_kelvin_offset(fem_mesh))
    source = build_esrf_fixed_magnetization_source(
        source_mesh,
        3,
        order=options.hdiv_order,
        curve_order=2,
        field_cf_algorithm=options.source_field_algorithm,
        field_tree_options={"theta": options.source_tree_theta},
    )
    bh_table = get_esrf_bh_table(3)
    points = observation_points()
    source_points = source_accuracy_points()
    output = options.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    source_h_samples = np.asarray(source.Field(source_points, algorithm="direct"), dtype=float)
    if source.field_cf_algorithm == "tree":
        _, source_tree_metrics = _source_tree_accuracy(source, source_points)
        source_tree_relative_rms = source_tree_metrics["tree_vs_direct_relative_rms"]
        if source_tree_relative_rms > options.source_tree_relative_tolerance:
            raise RuntimeError(
                "native tree source field differs from direct source field by "
                f"{source_tree_relative_rms:.3e}, exceeding "
                f"{options.source_tree_relative_tolerance:.3e}"
            )
    else:
        source_tree_metrics = {
            "algorithm": "direct",
            "tree_vs_direct_relative_rms": None,
            "sample_count": int(source_points.shape[0]),
        }

    source_mesh_accuracy = None
    if source_reference_vol is not None:
        if not source_reference_vol.is_file():
            raise FileNotFoundError(source_reference_vol)
        reference_source = build_esrf_fixed_magnetization_source(
            ng.Mesh(str(source_reference_vol)),
            3,
            order=options.hdiv_order,
            curve_order=2,
            field_cf_algorithm="direct",
        )
        source_mesh_accuracy = _source_mesh_direct_accuracy(
            source, reference_source, source_points
        )
        if (
            source_mesh_accuracy["relative_rms"]
            > options.source_mesh_relative_rms_tolerance
        ):
            raise RuntimeError(
                "fixed-M source mesh differs from retained direct source by "
                f"{source_mesh_accuracy['relative_rms']:.3e}, exceeding "
                f"{options.source_mesh_relative_rms_tolerance:.3e}"
            )
    if options.preflight:
        payload = {
            "schema": "radia.validation.esrf-hybrid-undulator-preflight.v1",
            "passed": True,
            "source_overlay": None if SOURCE_OVERLAY is None else str(SOURCE_OVERLAY),
            "implementation_paths": implementation_paths,
            "source_projection": _source_stats(source),
            "source_field_gap_rms_A_per_m": float(
                np.sqrt(np.mean(np.sum(source_h_samples * source_h_samples, axis=1)))
            ),
            "source_tree_accuracy": source_tree_metrics,
            "source_tree_relative_tolerance": options.source_tree_relative_tolerance,
            "source_mesh": {
                "path": str(source_vol),
                "sha256": _sha256(source_vol),
                "overridden": source_mesh_overridden,
                "reference_path": (
                    None if source_reference_vol is None else str(source_reference_vol)
                ),
                "reference_sha256": (
                    None if source_reference_vol is None else _sha256(source_reference_vol)
                ),
                "direct_field_accuracy": source_mesh_accuracy,
                "relative_rms_tolerance": options.source_mesh_relative_rms_tolerance,
            },
            "fem_mesh_contract": fem_contract,
            "kelvin_offset_m": list(kelvin_offset),
            "observation_points_m": points.tolist(),
            "source_accuracy_points_m": source_points.tolist(),
        }
        output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        _progress("preflight_complete", output=str(output), passed=True)
        return 0

    common = _checkpoint_contract(
        assets={"iron": _sha256(iron_vol), "source": _sha256(source_vol)},
        fem_mesh_sha256=_sha256(fem_mesh_path),
        bh_table_sha256=_json_digest(bh_table),
        observation_points_m=points.tolist(),
        hdiv_order=options.hdiv_order,
        fem_order=options.fem_order,
        hdiv_gram_eps=options.hdiv_gram_eps,
        hdiv_gram_backend=options.hdiv_gram_backend,
        hdiv_exact_dense_memory_mb=options.hdiv_exact_dense_memory_mb,
        source_tree_theta=options.source_tree_theta,
        source_field_algorithm=options.source_field_algorithm,
        hdiv_nonlinear_solver=options.hdiv_nonlinear_solver,
        source_reference_sha256=(
            None if source_reference_vol is None else _sha256(source_reference_vol)
        ),
        nonlinear_tolerance=options.nonlinear_tolerance,
        nonlinear_maximum_iterations=options.nonlinear_maximum_iterations,
        nonlinear_relaxation=options.nonlinear_relaxation,
    )
    fields: dict[str, np.ndarray] = {}
    diagnostics: dict[str, dict[str, object]] = {}
    engine_specs = (
        ("hdiv_mmm", lambda: _solve_hdiv(
            iron_mesh, source, bh_table, points, order=options.hdiv_order,
            gram_eps=options.hdiv_gram_eps,
            nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=options.nonlinear_maximum_iterations,
            nonlinear_solver=options.hdiv_nonlinear_solver,
            gram_backend=options.hdiv_gram_backend,
            exact_dense_memory_mb=options.hdiv_exact_dense_memory_mb)),
        ("reduced_a", lambda: _solve_reduced_a(
            fem_mesh, source, bh_table, points, order=options.fem_order,
            kelvin_radius=0.18, kelvin_center=kelvin_offset,
            linear_solver=options.reduced_a_solver,
            nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=options.nonlinear_maximum_iterations,
            nonlinear_relaxation=options.nonlinear_relaxation)),
        ("mixed_total_reduced_omega", lambda: _solve_mixed_omega(
            fem_mesh, source, bh_table, points, order=options.fem_order,
            kelvin_radius=0.18, kelvin_center=kelvin_offset,
            nonlinear_tolerance=options.nonlinear_tolerance,
            nonlinear_maximum_iterations=options.nonlinear_maximum_iterations,
            nonlinear_relaxation=options.nonlinear_relaxation,
            source_potential_tolerance=options.source_potential_tolerance)),
    )
    for engine, solve in engine_specs:
        checkpoint = output.with_suffix(f".{engine}.checkpoint.json")
        contract = _checkpoint_contract(**common, engine=engine)
        resumed = _read_checkpoint(checkpoint, contract) if options.resume else None
        if resumed is not None:
            fields[engine], diagnostics[engine] = resumed
            diagnostics[engine]["resumed_from_checkpoint"] = True
            _progress("engine_resumed", engine=engine, checkpoint=str(checkpoint))
            continue
        _progress("engine_start", engine=engine)
        fields[engine], diagnostics[engine] = solve()
        _write_checkpoint(checkpoint, contract, fields[engine], diagnostics[engine])
        _progress("engine_complete", engine=engine, checkpoint=str(checkpoint),
                  runtime_s=diagnostics[engine]["runtime_s"])

    three_engine_contract = require_static_electromagnet_three_engine_contract(diagnostics)
    pairwise = _pairwise_metrics(fields)
    maximum_pairwise = max(row["relative_rms"] for row in pairwise.values())
    nonlinear_converged = all(
        bool(row.get("nonlinear_stats", {}).get("converged", False))
        for row in diagnostics.values()
    )
    passed = nonlinear_converged and maximum_pairwise <= options.relative_rms_tolerance
    result = {
        "schema": "radia.validation.esrf-hybrid-undulator-three-engine.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "passed": bool(passed),
        "machine": platform.node(),
        "python": sys.version,
        "radia_version": importlib.metadata.version("radia"),
        "radia_module": str(Path(rad.__file__).resolve()),
        "source_overlay": None if SOURCE_OVERLAY is None else str(SOURCE_OVERLAY),
        "implementation_paths": implementation_paths,
        "formulation_contract": three_engine_contract,
        "shared_input_contract": {
            "case": "ESRF Example #3 hybrid permanent-magnet undulator",
            "fixed_magnetization_source": "same native C++ MagnetizationSource",
            "source_projection": _source_stats(source),
            "source_tree_accuracy": source_tree_metrics,
            "source_tree_relative_tolerance": options.source_tree_relative_tolerance,
            "source_mesh_sha256": _sha256(source_vol),
            "source_mesh_overridden": source_mesh_overridden,
            "source_reference_mesh_sha256": (
                None if source_reference_vol is None else _sha256(source_reference_vol)
            ),
            "source_mesh_direct_field_accuracy": source_mesh_accuracy,
            "source_mesh_relative_rms_tolerance": options.source_mesh_relative_rms_tolerance,
            "source_field_algorithm": options.source_field_algorithm,
            "hdiv_nonlinear_solver": options.hdiv_nonlinear_solver,
            "hdiv_gram_backend": options.hdiv_gram_backend,
            "hdiv_exact_dense_memory_mb": options.hdiv_exact_dense_memory_mb,
            "iron_response_mesh_sha256": _sha256(iron_vol),
            "fem_mesh_sha256": _sha256(fem_mesh_path),
            "fem_mesh_contract": fem_contract,
            "bh_table_sha256": _json_digest(bh_table),
            "finite_outer_air_box_forbidden": True,
            "global_reduced_omega_acceptance_forbidden": True,
        },
        "observation_points_m": points.tolist(),
        "engines": diagnostics,
        "fields_T": {name: value.tolist() for name, value in fields.items()},
        "pairwise": pairwise,
        "maximum_pairwise_relative_rms": maximum_pairwise,
        "relative_rms_tolerance": options.relative_rms_tolerance,
        "nonlinear_converged": nonlinear_converged,
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    _progress("validation_complete", output=str(output), passed=bool(passed),
              maximum_pairwise_relative_rms=maximum_pairwise)
    if not passed:
        raise RuntimeError(f"three-engine comparison did not pass; see {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
