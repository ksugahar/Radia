"""Certify C-type numerical accuracy by mesh convergence and three routes.

Each mesh level runs in a fresh Python process.  This both isolates NGSolve/C++
state and preserves the public ``run_three_engine.py`` contract.  Acceptance
requires converged nonlinear solves, contracting mesh increments for HDiv-MMM,
HCurl reduced-A, and H1 TOSCA mixed total/reduced Omega, a small fine-mesh three-route spread,
and a conservative combined numerical-uncertainty envelope.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUNNER = HERE / "run_three_engine.py"
ENGINES = ("hdiv_mmm", "reduced_a", "mixed_total_reduced_omega")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _relative_rms(reference: np.ndarray, candidate: np.ndarray) -> float:
    denominator = _rms(reference)
    if denominator <= 0.0:
        raise RuntimeError("zero reference field in convergence analysis")
    return _rms(candidate - reference) / denominator


def _load_level(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "radia.validation.c-type-formulation-comparison.v4":
        raise RuntimeError(f"unexpected level result schema: {path}")
    if set(payload.get("engines", {})) != set(ENGINES):
        raise RuntimeError(f"all three formulations are required: {path}")
    if not payload.get("nonlinear_converged", False):
        raise RuntimeError(f"nonlinear formulation did not converge: {path}")
    return payload


def _run_level(
    mesh_dir: Path,
    output: Path,
    options: argparse.Namespace,
) -> dict:
    command = [
        sys.executable,
        str(RUNNER),
        "--mesh-dir",
        str(mesh_dir),
        "--output",
        str(output),
        "--mode",
        "nonlinear",
        "--hdiv-order",
        str(options.hdiv_order),
        "--hdiv-gram-eps",
        str(options.hdiv_gram_eps),
        "--fem-order",
        str(options.fem_order),
        "--reduced-a-solver",
        options.reduced_a_solver,
        "--reduced-a-relax",
        str(options.reduced_a_relax),
        "--nonlinear-tolerance",
        str(options.nonlinear_tolerance),
        "--nonlinear-maximum-iterations",
        str(options.nonlinear_maximum_iterations),
        "--threads",
        str(options.threads),
        "--gap-core-half-length",
        str(options.gap_core_half_length),
        "--relative-rms-tolerance",
        "1.0",
    ]
    if options.resume:
        command.append("--resume")
    log_path = output.with_suffix(".log")
    print(
        json.dumps(
            {"event": "mesh_level_start", "mesh_dir": str(mesh_dir)},
            sort_keys=True,
        ),
        flush=True,
    )
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            cwd=HERE,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
        returncode = process.wait()
    if returncode != 0 or not output.is_file():
        raise RuntimeError(
            f"three-engine level failed with exit code {returncode}; see {log_path}"
        )
    return _load_level(output)


def _selector(payload: dict) -> np.ndarray:
    points = np.asarray(payload["observation_points_m"], dtype=float)
    half_length = float(payload["gap_core_half_length_m"])
    return np.abs(points[:, 0]) <= half_length + 1e-14


def _engine_convergence(
    payloads: list[dict],
    level_names: list[str],
    engine: str,
    refinement_ratio: float,
) -> dict[str, object]:
    selector = _selector(payloads[-1])
    fields = [
        np.asarray(payload["median_plane_projected_fields_T"][engine], dtype=float)[
            selector
        ]
        for payload in payloads
    ]
    coarse_increment = _rms(fields[-2] - fields[-3])
    fine_increment = _rms(fields[-1] - fields[-2])
    fine_norm = _rms(fields[-1])
    if fine_norm <= 0.0:
        raise RuntimeError(f"zero fine field for {engine}")
    contracting = coarse_increment > 0.0 and fine_increment < coarse_increment
    contraction = (
        fine_increment / coarse_increment if coarse_increment > 0.0 else math.inf
    )
    observed_order = (
        math.log(coarse_increment / fine_increment) / math.log(refinement_ratio)
        if contracting and fine_increment > 0.0
        else None
    )
    richardson = None
    if observed_order is not None:
        denominator = refinement_ratio**observed_order - 1.0
        if denominator > 0.0:
            richardson = fine_increment / denominator / fine_norm
    last_step = fine_increment / fine_norm
    uncertainty = max(last_step, richardson if richardson is not None else math.inf)
    return {
        "convergence_levels": level_names[-3:],
        "penultimate_increment_absolute_rms_T": coarse_increment,
        "last_increment_absolute_rms_T": fine_increment,
        "last_increment_relative_rms": last_step,
        "contraction_ratio": contraction,
        "contracting": contracting,
        "observed_order": observed_order,
        "richardson_relative_uncertainty": richardson,
        "conservative_relative_uncertainty": uncertainty,
    }


def _replicate_metrics(reference: dict, replicate: dict) -> dict[str, object]:
    reference_machine = str(reference.get("machine", "")).strip()
    replicate_machine = str(replicate.get("machine", "")).strip()
    if not reference_machine or not replicate_machine:
        raise RuntimeError("replicate results must record both machine names")
    if reference_machine.casefold() == replicate_machine.casefold():
        raise RuntimeError("accuracy replication must use an independent host")
    if reference["mesh_result_sha256"] != replicate["mesh_result_sha256"]:
        raise RuntimeError("replicate result uses a different fine mesh contract")
    if reference["comparison_contract"] != replicate["comparison_contract"]:
        raise RuntimeError("replicate result uses a different formulation contract")
    if reference.get("radia_version") != replicate.get("radia_version"):
        raise RuntimeError("replicate result uses a different Radia version")
    if reference.get("engine_checkpoint_contracts") != replicate.get(
        "engine_checkpoint_contracts"
    ):
        raise RuntimeError("replicate result uses different implementation inputs")
    selector = _selector(reference)
    rows = {}
    for engine in ENGINES:
        left = np.asarray(
            reference["median_plane_projected_fields_T"][engine], dtype=float
        )[selector]
        right = np.asarray(
            replicate["median_plane_projected_fields_T"][engine], dtype=float
        )[selector]
        rows[engine] = {"relative_rms": _relative_rms(left, right)}
    return {
        "reference_machine": reference_machine,
        "replicate_machine": replicate_machine,
        "engines": rows,
        "maximum_relative_rms": max(row["relative_rms"] for row in rows.values()),
    }


def analyze(
    manifest_path: Path,
    level_paths: list[Path],
    output: Path,
    options: argparse.Namespace,
) -> dict[str, object]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "radia.validation.c-type-cubit-mesh-family.v1":
        raise RuntimeError("unexpected mesh family schema")
    if not manifest.get("passed", False):
        raise RuntimeError("mesh family contract is not passing")
    if len(level_paths) < 3:
        raise RuntimeError("at least three mesh levels are required")
    payloads = [_load_level(path) for path in level_paths]
    points = payloads[0]["observation_points_m"]
    for payload in payloads[1:]:
        if payload["observation_points_m"] != points:
            raise RuntimeError("mesh levels use different physical observations")

    ratio = float(manifest["refinement_ratio"])
    level_names = [str(row["name"]) for row in manifest["levels"]]
    convergence = {
        engine: _engine_convergence(payloads, level_names, engine, ratio)
        for engine in ENGINES
    }
    convergence_passed = all(
        row["contracting"]
        and row["observed_order"] is not None
        and row["observed_order"] >= options.minimum_observed_order
        and row["conservative_relative_uncertainty"]
        <= options.mesh_uncertainty_tolerance
        for row in convergence.values()
    )

    fine = payloads[-1]
    selector = _selector(fine)
    fine_fields = {
        engine: np.asarray(
            fine["median_plane_projected_fields_T"][engine], dtype=float
        )[selector]
        for engine in ENGINES
    }
    consensus = np.mean(np.stack(tuple(fine_fields.values())), axis=0)
    consensus_norm = _rms(consensus)
    if consensus_norm <= 0.0:
        raise RuntimeError("zero fine-mesh consensus field")
    consensus_deviation = {
        engine: _rms(field - consensus) / consensus_norm
        for engine, field in fine_fields.items()
    }
    maximum_consensus_deviation = max(consensus_deviation.values())
    final_pairwise_spread = float(
        fine["maximum_gap_core_pairwise_relative_rms"]
    )
    maximum_discretization_uncertainty = max(
        float(row["conservative_relative_uncertainty"])
        for row in convergence.values()
    )
    combined_uncertainty = (
        maximum_discretization_uncertainty + maximum_consensus_deviation
    )

    replicate = None
    reproducibility_passed = False
    if options.replicate_final_result is not None:
        replicate_payload = _load_level(options.replicate_final_result.resolve())
        replicate = _replicate_metrics(fine, replicate_payload)
        reproducibility_passed = (
            replicate["maximum_relative_rms"]
            <= options.reproducibility_tolerance
        )

    final_spread_passed = (
        final_pairwise_spread <= options.cross_formulation_tolerance
    )
    combined_uncertainty_passed = (
        combined_uncertainty <= options.combined_uncertainty_tolerance
    )
    passed = bool(
        convergence_passed
        and final_spread_passed
        and combined_uncertainty_passed
        and reproducibility_passed
    )
    result = {
        "schema": "radia.validation.c-type-absolute-accuracy-certificate.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "passed": passed,
        "machine": platform.node(),
        "claim": {
            "quantity": "gauge-invariant B in the C-yoke useful gap core",
            "analytic_absolute_truth_claimed": False,
            "certified_statement": (
                "three independent formulations converge toward one common "
                "mesh-refined field within the reported numerical envelope"
            ),
            "scope": "nonlinear order-2 Cubit/ACIS C-type model with Kelvin open boundary",
        },
        "mesh_family": str(manifest_path),
        "mesh_family_sha256": sha256(manifest_path),
        "level_results": [
            {
                "name": row["name"],
                "scale": row["scale"],
                "result": str(path),
                "result_sha256": sha256(path),
                "engines": payload["engines"],
            }
            for row, path, payload in zip(manifest["levels"], level_paths, payloads)
        ],
        "refinement_ratio": ratio,
        "convergence_levels": level_names[-3:],
        "engine_convergence": convergence,
        "fine_mesh_pairwise": fine["pairwise_median_projected_gap_core"],
        "fine_mesh_maximum_pairwise_relative_rms": final_pairwise_spread,
        "fine_mesh_consensus_relative_deviation": consensus_deviation,
        "maximum_consensus_relative_deviation": maximum_consensus_deviation,
        "maximum_discretization_relative_uncertainty": (
            maximum_discretization_uncertainty
        ),
        "combined_relative_numerical_uncertainty": combined_uncertainty,
        "reproducibility": replicate,
        "independent_host_replicate_required": True,
        "thresholds": {
            "minimum_observed_order": options.minimum_observed_order,
            "mesh_uncertainty_tolerance": options.mesh_uncertainty_tolerance,
            "cross_formulation_tolerance": options.cross_formulation_tolerance,
            "combined_uncertainty_tolerance": options.combined_uncertainty_tolerance,
            "reproducibility_tolerance": options.reproducibility_tolerance,
        },
        "checks": {
            "all_three_nonlinear_formulations_converged": True,
            "mesh_convergence_passed": convergence_passed,
            "fine_mesh_cross_formulation_passed": final_spread_passed,
            "combined_uncertainty_passed": combined_uncertainty_passed,
            "independent_host_reproducibility_passed": reproducibility_passed,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not passed:
        raise RuntimeError(f"C-type accuracy certificate failed; see {output}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-family", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hdiv-order", type=int, choices=(1, 2), default=2)
    parser.add_argument("--hdiv-gram-eps", type=float, default=1.0e-14)
    parser.add_argument("--fem-order", type=int, default=2)
    parser.add_argument(
        "--reduced-a-solver",
        choices=("direct", "bddc", "ams", "auto"),
        default="direct",
    )
    parser.add_argument("--reduced-a-relax", type=float, default=0.1)
    parser.add_argument("--nonlinear-tolerance", type=float, default=2.0e-5)
    parser.add_argument("--nonlinear-maximum-iterations", type=int, default=80)
    parser.add_argument("--threads", type=int, default=38)
    parser.add_argument("--gap-core-half-length", type=float, default=0.010)
    parser.add_argument("--minimum-observed-order", type=float, default=0.5)
    parser.add_argument("--mesh-uncertainty-tolerance", type=float, default=0.005)
    parser.add_argument("--cross-formulation-tolerance", type=float, default=0.005)
    parser.add_argument(
        "--combined-uncertainty-tolerance", type=float, default=0.010
    )
    parser.add_argument("--replicate-final-result", type=Path)
    parser.add_argument("--reproducibility-tolerance", type=float, default=1.0e-9)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--analyze-only",
        action="store_true",
        help="consume level JSON files already present beside --output",
    )
    options = parser.parse_args()
    manifest_path = options.mesh_family.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    output = options.output.resolve()
    level_dir = output.parent / "levels"
    level_dir.mkdir(parents=True, exist_ok=True)
    level_paths = []
    for row in manifest["levels"]:
        result_path = level_dir / f"{row['name']}.json"
        if not options.analyze_only:
            mesh_dir = manifest_path.parent / row["relative_directory"]
            _run_level(mesh_dir, result_path, options)
        elif not result_path.is_file():
            raise FileNotFoundError(result_path)
        level_paths.append(result_path)
    analyze(manifest_path, level_paths, output, options)


if __name__ == "__main__":
    main()
