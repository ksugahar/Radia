"""Tune HDiv-VIM HACApK parameters with an accuracy-constrained Optuna study.

This is a compute-host validation benchmark, not a pytest test.  It builds the
NGSolve HDiv space and Radia charge map once, then repeatedly rebuilds the same
C++ ChargeGram H-matrix.  Each trial minimizes a measured build/apply workload
while constraining the physical ``B.T @ G @ B`` action against the accepted
current parameters.  A tighter H-matrix remains an independent diagnostic.

Example (run on hibino first, or on mdx only behind an idle CI queue)::

    python validation_test/feec/tune_hdiv_hmatrix_optuna.py \
        --topology tet --order 2 --cells 5 --cyclic-sectors 4 \
        --trials 24 --expected-applies 200 --output C:/temp/hdiv_hmat_tune.json

Install the pinned upstream optimizer with ``pip install radia[optuna-upstream]``.
The standalone ``radia-optuna`` distribution supplies the MATLAB implementation;
this NGSolve-object benchmark uses upstream Python Optuna 4.9.0 directly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
import warnings
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import ngsolve as ng
import numpy as np
import radia._radia_pybind as radia_pybind
from ngsolve.meshes import MakeStructured3DMesh

import radia
from radia.vim import ChargeGram

SCHEMA = "radia.validation.hdiv-hmatrix-optuna.v1"
SUPPORTED_OPTUNA_VERSION = "4.9.0"


def _csv_ints(value: str) -> list[int]:
    values = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not values or any(item <= 0 for item in values):
        raise argparse.ArgumentTypeError("expected comma-separated positive integers")
    return values


def _csv_floats(value: str) -> list[float]:
    values = [float(item.strip()) for item in value.split(",") if item.strip()]
    if not values or any(not math.isfinite(item) or item <= 0.0 for item in values):
        raise argparse.ArgumentTypeError("expected comma-separated positive finite values")
    return values


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _compute_host(*, allow_local_diagnostic: bool) -> tuple[str, bool]:
    hostname = platform.node()
    short = hostname.lower().split(".", 1)[0]
    accepted = short in {"mdx", "hibino"}
    if not accepted and not allow_local_diagnostic:
        raise SystemExit(
            "HDiv H-matrix tuning is a validation workload and must run on "
            "hibino first, or on mdx only when its CI queue is idle. Pass "
            "--allow-local-diagnostic only for a bounded workflow "
            f"smoke test; got hostname {hostname!r}.")
    return hostname, not accepted


def _mesh(topology: str, cells: int, curve_order: int, cyclic_sectors: int):
    kwargs: dict[str, Any] = {"nx": cells, "ny": cells, "nz": cells}
    if topology == "tet":
        kwargs["hexes"] = False
    elif topology == "hex":
        kwargs["hexes"] = True
    elif topology == "wedge":
        kwargs["prism"] = True
    else:
        raise ValueError(f"unsupported topology {topology!r}")

    if cyclic_sectors > 1:
        # Keep the retained body away from the z axis so every cyclic image is
        # geometrically distinct.  The benchmark tunes the real image kernel;
        # it does not rely on a synthetic point-cloud surrogate.
        kwargs["mapping"] = lambda x, y, z: (
            1.2 + 0.6 * x, -0.25 + 0.5 * y, -0.25 + 0.5 * z)
    mesh = MakeStructured3DMesh(**kwargs)
    if curve_order > 1:
        mesh.Curve(curve_order)
    return mesh


def _cyclic_images(sectors: int, alternating: bool) -> dict[str, list[float] | list[int]]:
    if sectors < 1:
        raise ValueError("cyclic_sectors must be positive")
    if alternating and sectors % 2:
        raise ValueError("alternating cyclic images require an even sector count")
    if sectors == 1:
        return {"image_masks": [], "image_signs": [], "image_rot_angle": []}
    return {
        "image_masks": [0] * (sectors - 1),
        "image_signs": [(-1.0 if alternating and index % 2 else 1.0)
                        for index in range(1, sectors)],
        "image_rot_angle": [2.0 * math.pi * index / sectors
                            for index in range(1, sectors)],
    }


def _apply_probes(gram, probes: np.ndarray) -> np.ndarray:
    return np.vstack([
        np.asarray(gram.apply_configured_demag(probe, True), dtype=float)
        for probe in probes
    ])


def _relative_action_error(candidate: np.ndarray, reference: np.ndarray) -> float:
    denominator = float(np.linalg.norm(reference))
    if denominator == 0.0:
        raise RuntimeError("tight H-matrix baseline returned a zero probe action")
    return float(np.linalg.norm(candidate - reference) / denominator)


def _timed_apply(gram, probes: np.ndarray, repeats: int) -> tuple[np.ndarray, dict[str, float]]:
    first_start = time.perf_counter()
    first = _apply_probes(gram, probes)
    first_s = time.perf_counter() - first_start
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        repeated = _apply_probes(gram, probes)
        samples.append(time.perf_counter() - start)
        if not np.array_equal(first, repeated):
            raise RuntimeError("repeated symmetric demag action changed numerically")
    return first, {
        "probe_batch_first_s": first_s,
        "probe_batch_median_s": float(np.median(samples)),
        "apply_per_probe_median_s": float(np.median(samples)) / len(probes),
    }


def _trial_record(trial) -> dict[str, Any]:
    return {
        "number": int(trial.number),
        "state": trial.state.name,
        "params": _jsonable(trial.params),
        "value": None if trial.value is None else float(trial.value),
        "user_attrs": _jsonable(trial.user_attrs),
        "datetime_start": (
            None if trial.datetime_start is None else trial.datetime_start.isoformat()),
        "datetime_complete": (
            None if trial.datetime_complete is None else trial.datetime_complete.isoformat()),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", choices=("tet", "hex", "wedge"), default="tet")
    parser.add_argument("--order", type=int, choices=(1, 2), default=2)
    parser.add_argument("--curve-order", type=int, choices=(1, 2), default=1)
    parser.add_argument("--cells", type=int, default=3)
    parser.add_argument("--cyclic-sectors", type=int, default=1)
    parser.add_argument("--alternating", action="store_true")
    parser.add_argument("--trials", type=int, default=16)
    parser.add_argument("--sampler", choices=("tpe", "grid"), default="tpe")
    parser.add_argument("--seed", type=int, default=20260826)
    parser.add_argument("--probes", type=int, default=4)
    parser.add_argument("--apply-repeats", type=int, default=3)
    parser.add_argument("--expected-applies", type=int, default=100)
    parser.add_argument("--max-action-relative-delta", type=float, default=1.0e-8)
    parser.add_argument("--minimum-speedup", type=float, default=1.03)
    parser.add_argument("--leaf-sizes", type=_csv_ints, default=_csv_ints("16,24,32,48,64,96,128"))
    parser.add_argument("--etas", type=_csv_floats, default=_csv_floats("1,1.5,2,3,4"))
    parser.add_argument("--eps-values", type=_csv_floats,
                        default=_csv_floats("1e-10,3e-11,1e-11,3e-12,1e-12"))
    parser.add_argument("--reference-eps", type=float, default=1.0e-14)
    parser.add_argument("--reference-leaf-size", type=int, default=128)
    parser.add_argument("--reference-eta", type=float, default=2.0)
    parser.add_argument("--current-eps", type=float, default=1.0e-12)
    parser.add_argument("--current-leaf-size", type=int, default=32)
    parser.add_argument("--current-eta", type=float, default=2.0)
    parser.add_argument("--study-name", default="hdiv-hmatrix")
    parser.add_argument("--storage", help="Optional Optuna storage URL, for example sqlite:///study.db")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-local-diagnostic", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    if min(args.cells, args.trials, args.probes, args.apply_repeats) <= 0:
        raise SystemExit("cells, trials, probes, and apply-repeats must be positive")
    if (args.expected_applies < 0 or args.max_action_relative_delta <= 0.0
            or args.minimum_speedup < 1.0):
        raise SystemExit(
            "expected-applies must be non-negative, the delta limit positive, "
            "and minimum-speedup at least one")

    hostname, local_diagnostic = _compute_host(
        allow_local_diagnostic=args.allow_local_diagnostic)
    try:
        import optuna
    except ImportError as exc:
        raise SystemExit(
            "Optuna is required; install pip install radia[optuna-upstream]") from exc
    if optuna.__version__ != SUPPORTED_OPTUNA_VERSION:
        raise SystemExit(
            f"expected pinned optuna=={SUPPORTED_OPTUNA_VERSION}, got {optuna.__version__}")

    images = _cyclic_images(args.cyclic_sectors, args.alternating)
    kernel_path = Path(radia_pybind.__file__).resolve()
    kernel_sha256 = _sha256(kernel_path)
    script_sha256 = _sha256(Path(__file__).resolve())
    mesh = _mesh(args.topology, args.cells, args.curve_order, args.cyclic_sectors)
    vertex_counts = sorted({len(element.vertices) for element in mesh.Elements(ng.VOL)})
    expected_vertices = {"tet": [4], "hex": [8], "wedge": [6]}[args.topology]
    if vertex_counts != expected_vertices:
        raise RuntimeError(
            f"{args.topology}: expected volume vertex counts {expected_vertices}, got {vertex_counts}")

    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=args.order)
        basis_start = time.perf_counter()
        _, gram, _ = ChargeGram(
            fes,
            eps=args.reference_eps,
            leafsize=args.reference_leaf_size,
            eta=args.reference_eta,
            curve_order=(args.curve_order if args.curve_order > 1 else None),
            _build_hmatrix=False,
            **images,
        )
        basis_s = time.perf_counter() - basis_start

    rng = np.random.default_rng(args.seed)
    probes = np.ascontiguousarray(
        rng.standard_normal((args.probes, int(fes.ndof))), dtype=np.float64)

    reference_start = time.perf_counter()
    gram.build_hmatrix(
        eps=args.reference_eps,
        leaf=args.reference_leaf_size,
        eta=args.reference_eta)
    reference_build_external_s = time.perf_counter() - reference_start
    reference_action, reference_apply = _timed_apply(
        gram, probes, args.apply_repeats)
    reference_stats = dict(gram.stats())

    def evaluate(
            eps: float, leaf_size: int, eta: float
            ) -> tuple[dict[str, Any], np.ndarray]:
        build_start = time.perf_counter()
        gram.build_hmatrix(eps=eps, leaf=leaf_size, eta=eta)
        build_external_s = time.perf_counter() - build_start
        candidate, apply_timing = _timed_apply(
            gram, probes, args.apply_repeats)
        relative_error = _relative_action_error(candidate, reference_action)
        workload_s = (
            build_external_s
            + args.expected_applies * apply_timing["apply_per_probe_median_s"])
        return {
            "params": {"eps": eps, "leaf_size": leaf_size, "eta": eta},
            "action_relative_error_to_tight_reference": relative_error,
            "build_external_s": build_external_s,
            "apply_timing": _jsonable(apply_timing),
            "hmat_stats": _jsonable(dict(gram.stats())),
            "workload_s": workload_s,
        }, candidate

    current, current_action = evaluate(
        args.current_eps, args.current_leaf_size, args.current_eta)
    current["action_relative_delta_from_current"] = 0.0
    current["accuracy_constraint"] = -args.max_action_relative_delta

    def constraints_func(frozen_trial):
        return [float(frozen_trial.user_attrs.get("accuracy_constraint", math.inf))]

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    if args.sampler == "grid":
        sampler = optuna.samplers.GridSampler(
            {
                "eps": args.eps_values,
                "leaf_size": args.leaf_sizes,
                "eta": args.etas,
            },
            seed=args.seed,
        )
        trial_budget = min(
            args.trials,
            len(args.eps_values) * len(args.leaf_sizes) * len(args.etas))
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", optuna.exceptions.ExperimentalWarning)
            sampler = optuna.samplers.TPESampler(
                seed=args.seed,
                n_startup_trials=min(5, args.trials),
                constraints_func=constraints_func,
            )
        trial_budget = args.trials
    study = optuna.create_study(
        study_name=args.study_name,
        direction="minimize",
        sampler=sampler,
        storage=args.storage,
        load_if_exists=bool(args.storage),
    )
    case_signature = json.dumps({
        "topology": args.topology,
        "order": args.order,
        "curve_order": args.curve_order,
        "cells": args.cells,
        "cyclic_sectors": args.cyclic_sectors,
        "alternating": args.alternating,
        "seed": args.seed,
        "probes": args.probes,
        "apply_repeats": args.apply_repeats,
        "expected_applies": args.expected_applies,
        "max_action_relative_delta": args.max_action_relative_delta,
        "minimum_speedup": args.minimum_speedup,
        "reference_eps": args.reference_eps,
        "reference_leaf_size": args.reference_leaf_size,
        "reference_eta": args.reference_eta,
        "eps_values": args.eps_values,
        "leaf_sizes": args.leaf_sizes,
        "etas": args.etas,
        "kernel_sha256": kernel_sha256,
        "script_sha256": script_sha256,
    }, sort_keys=True)
    stored_signature = study.user_attrs.get("case_signature")
    if stored_signature is not None and stored_signature != case_signature:
        raise SystemExit(
            "the resumed Optuna study belongs to a different mesh/search contract; "
            "use a new --study-name or storage file")
    study.set_user_attr("case_signature", case_signature)
    current_in_search = (
        args.current_eps in args.eps_values
        and args.current_leaf_size in args.leaf_sizes
        and args.current_eta in args.etas)

    def objective(trial) -> float:
        eps = trial.suggest_categorical("eps", args.eps_values)
        leaf_size = trial.suggest_categorical("leaf_size", args.leaf_sizes)
        eta = trial.suggest_categorical("eta", args.etas)
        measured, candidate_action = evaluate(eps, leaf_size, eta)
        relative_delta = _relative_action_error(candidate_action, current_action)
        measured["action_relative_delta_from_current"] = relative_delta
        measured["accuracy_constraint"] = (
            relative_delta - args.max_action_relative_delta)
        for key, value in measured.items():
            if key != "params":
                trial.set_user_attr(key, value)
        return float(measured["workload_s"])

    study.optimize(
        objective, n_trials=trial_budget, gc_after_trial=True,
        catch=(RuntimeError,))
    feasible = [
        trial for trial in study.trials
        if trial.state == optuna.trial.TrialState.COMPLETE
        and float(trial.user_attrs.get("accuracy_constraint", math.inf)) <= 0.0
    ]
    best = min(feasible, key=lambda trial: float(trial.value)) if feasible else None
    recommended: dict[str, Any] | None = None
    if float(current["accuracy_constraint"]) <= 0.0:
        recommended = {"source": "current", **current}
    changed = [
        trial for trial in feasible
        if trial.params != current["params"]
    ]
    best_changed = min(
        changed, key=lambda trial: float(trial.value)) if changed else None
    if best_changed is not None and (
            recommended is None
            or float(best_changed.value) * args.minimum_speedup
            <= float(recommended["workload_s"])):
        recommended = {
            "source": "optuna_trial",
            "trial_number": int(best_changed.number),
            "params": _jsonable(best_changed.params),
            **_jsonable(best_changed.user_attrs),
        }

    payload = {
        "schema": SCHEMA,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "hostname": hostname,
        "local_diagnostic": local_diagnostic,
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "radia_version": radia.__version__,
        "ngsolve_version": ng.__version__,
        "optuna_version": optuna.__version__,
        "ngsolve_threads": int(ng.ngsglobals.numthreads),
        "implementation": {
            "radia_pybind_path": str(kernel_path),
            "radia_pybind_sha256": kernel_sha256,
            "script_sha256": script_sha256,
        },
        "case": {
            "topology": args.topology,
            "order": args.order,
            "curve_order": int(mesh.GetCurveOrder()),
            "cells_per_axis": args.cells,
            "elements": int(mesh.ne),
            "hdiv_dof": int(fes.ndof),
            "charge_dof": int(gram.ndof()),
            "cyclic_sectors": args.cyclic_sectors,
            "alternating": args.alternating,
            "charge_basis_and_operator_shell_s": basis_s,
        },
        "objective": {
            "formula": "build_external_s + expected_applies * apply_per_probe_median_s",
            "expected_applies": args.expected_applies,
            "probe_count": args.probes,
            "apply_repeats": args.apply_repeats,
            "max_action_relative_delta": args.max_action_relative_delta,
            "minimum_speedup": args.minimum_speedup,
            "seed": args.seed,
        },
        "search_space": {
            "sampler": args.sampler,
            "eps": args.eps_values,
            "leaf_size": args.leaf_sizes,
            "eta": args.etas,
            "current_parameters_in_search": current_in_search,
        },
        "reference": {
            "eps": args.reference_eps,
            "leaf_size": args.reference_leaf_size,
            "eta": args.reference_eta,
            "build_external_s": reference_build_external_s,
            "apply_timing": reference_apply,
            "hmat_stats": _jsonable(reference_stats),
        },
        "current_parameters": _jsonable(current),
        "best_feasible_trial": None if best is None else _trial_record(best),
        "recommended_parameters": recommended,
        "trials": [_trial_record(trial) for trial in study.trials],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    recommendation_summary = None if recommended is None else {
        "source": recommended["source"],
        "trial_number": recommended.get("trial_number"),
        "params": recommended["params"],
        "workload_s": recommended["workload_s"],
        "action_relative_delta_from_current": (
            recommended["action_relative_delta_from_current"]),
        "action_relative_error_to_tight_reference": (
            recommended["action_relative_error_to_tight_reference"]),
    }
    print(json.dumps({
        "output": str(args.output),
        "feasible_trials": len(feasible),
        "recommended": recommendation_summary,
    }, indent=2))
    return 0 if recommended is not None else 2


if __name__ == "__main__":
    raise SystemExit(main())
