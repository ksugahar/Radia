"""Read the three-engine comparison across the C-type GAP family.

The gap family (``build_mesh_family.py --gap-layers 6 12 24``) refines ONLY
the gap air: N parallel slabs through the gap and an in-gap size of h/N,
with the iron, air and Kelvin sizes held at one base level.  It is not a
scale family and must not be fed to ``run_mesh_convergence.py``:

* HDiv-MMM uses the iron mesh alone, so it is the SAME at every level.  Here
  that is a gate, not a convergence quantity: the air-mesh-free route must
  not move when only the gap air changes.
* The FEM routes' increments between levels measure gap resolution only.
  Their observed order is in the layer count N, and the Richardson estimate
  bounds the gap-resolution error, not the discretisation error of the whole
  model.

Reported per family (one FEM order):

* the HDiv-MMM identity across levels;
* for reduced-A and mixed Omega, the gap-core increments between consecutive
  levels, contraction, observed order and Richardson estimate (the estimator
  of ``run_mesh_convergence.py``);
* each FEM route's distance from HDiv-MMM per level, and whether it shrinks
  as the gap is refined;
* the pairwise spread per level.

The gate is: HDiv-MMM identical, both FEM routes contracting.  Nothing here
is an absolute-accuracy claim, and a passing gate on one family says nothing
about a family built with other base sizes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_mesh_convergence as convergence  # noqa: E402

ENGINES = convergence.ENGINES
FEM_ENGINES = ("reduced_a", "mixed_total_reduced_omega")
LEVEL_SCHEMA = "radia.validation.c-type-formulation-comparison.v4"
FAMILY_SCHEMA = "radia.validation.c-type-cubit-gap-family.v1"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_level(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != LEVEL_SCHEMA:
        raise RuntimeError(f"unexpected level result schema: {path}")
    if set(payload.get("engines", {})) != set(ENGINES):
        raise RuntimeError(f"all three formulations are required: {path}")
    if payload.get("mode") == "nonlinear" and not payload.get("nonlinear_converged", False):
        raise RuntimeError(f"nonlinear formulation did not converge: {path}")
    return payload


_PER_LEVEL_KEYS = ("iron_vol_sha256", "kelvin_domain_vol_sha256",
                   "observation_points", "implementation_sha256")


def _fem_contract(payload: dict) -> dict:
    """What must agree between levels for the family to be one experiment.

    The engine checkpoint contracts carry the level's own mesh hashes; those
    are stripped.  Everything else -- mode, FEM order, the reduced-A linear
    solver, tolerances, the comparison contract, the Radia version -- must be
    identical, or the levels are not one experiment.
    """
    checkpoints = {
        engine: {key: value for key, value in contract.items()
                 if key not in _PER_LEVEL_KEYS}
        for engine, contract in (payload.get("engine_checkpoint_contracts") or {}).items()
    }
    return {
        "mode": payload.get("mode"),
        "comparison_contract": payload.get("comparison_contract"),
        "engine_checkpoint_contracts": checkpoints,
        "radia_version": payload.get("radia_version"),
    }


def analyze(manifest: dict, results: dict[str, dict], result_paths: dict[str, Path],
            *, hdiv_identity_tolerance: float) -> dict:
    if manifest.get("schema") != FAMILY_SCHEMA:
        raise RuntimeError("not a C-type gap family manifest")
    if not manifest.get("passed", False):
        raise RuntimeError("gap family mesh contract is not passing")
    names = [str(row["name"]) for row in manifest["levels"]]
    if len(names) < 3:
        raise RuntimeError("at least three gap levels are required")
    missing = [name for name in names if name not in results]
    if missing:
        raise RuntimeError(f"no result for gap levels {missing}")
    payloads = [results[name] for name in names]

    points = payloads[0]["observation_points_m"]
    for payload in payloads[1:]:
        if payload["observation_points_m"] != points:
            raise RuntimeError("gap levels use different observation points")
    contract = _fem_contract(payloads[0])
    for name, payload in zip(names[1:], payloads[1:]):
        if _fem_contract(payload) != contract:
            raise RuntimeError(f"gap level {name} was run under a different "
                               "formulation contract or Radia version")
    for name, row in zip(names, manifest["levels"]):
        mesh_sha = row["mesh_result_sha256"]
        if results[name]["mesh_result_sha256"] != mesh_sha:
            raise RuntimeError(f"result for {name} was not computed on the "
                               "manifest's mesh contract")

    ratio = float(manifest["gap_refinement_ratio"])
    selector = convergence._selector(payloads[-1])
    fields = {
        engine: [np.asarray(payload["median_plane_projected_fields_T"][engine],
                            dtype=float)[selector] for payload in payloads]
        for engine in ENGINES
    }

    hdiv_reference = fields["hdiv_mmm"][-1]
    hdiv_identity = [convergence._relative_rms(hdiv_reference, field)
                     for field in fields["hdiv_mmm"]]
    hdiv_identical = max(hdiv_identity) <= hdiv_identity_tolerance

    fem_convergence = {
        engine: convergence._engine_convergence(payloads, names, engine, ratio)
        for engine in FEM_ENGINES
    }
    distance_to_hdiv = {
        engine: [convergence._relative_rms(hdiv_reference, field)
                 for field in fields[engine]]
        for engine in FEM_ENGINES
    }
    approaching_hdiv = {
        engine: all(later < earlier for earlier, later in zip(rows, rows[1:]))
        for engine, rows in distance_to_hdiv.items()
    }

    levels = []
    for name, row, payload in zip(names, manifest["levels"], payloads):
        levels.append({
            "name": name,
            "gap_layers": int(row["gap_layers"]),
            "gap_elements": int(row["gap_elements"]),
            "result": str(result_paths[name]),
            "result_sha256": sha256(result_paths[name]),
            "machine": payload.get("machine"),
            "maximum_gap_core_pairwise_relative_rms":
                float(payload["maximum_gap_core_pairwise_relative_rms"]),
            "pairwise_median_projected_gap_core":
                payload["pairwise_median_projected_gap_core"],
            "engines": {engine: {key: payload["engines"][engine].get(key)
                                 for key in ("ndof", "mesh_elements", "runtime_s")}
                        for engine in ENGINES},
        })

    passed = bool(hdiv_identical and all(
        row["contracting"] for row in fem_convergence.values()))
    return {
        "schema": "radia.validation.c-type-gap-family-report.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "passed": passed,
        "machine": platform.node(),
        "claim": {
            "refined": manifest.get("refined"),
            "held_fixed": manifest.get("held_fixed"),
            "statement": ("gap-only refinement: the FEM routes' increments and "
                          "Richardson estimates bound the gap-resolution error "
                          "of THIS base mesh; HDiv-MMM is the fixed air-mesh-free "
                          "route, not a converged truth"),
            "analytic_absolute_truth_claimed": False,
        },
        "mode": contract["mode"],
        "radia_version": contract["radia_version"],
        "fem_contract": contract["engine_checkpoint_contracts"],
        "gap_refinement_ratio": ratio,
        "levels": levels,
        "hdiv_identity": {
            "relative_rms_to_finest": hdiv_identity,
            "tolerance": hdiv_identity_tolerance,
            "identical": hdiv_identical,
        },
        "fem_convergence_in_gap_layers": fem_convergence,
        "distance_to_hdiv_relative_rms": distance_to_hdiv,
        "distance_to_hdiv_shrinks_every_level": approaching_hdiv,
        "checks": {
            "hdiv_mmm_unchanged_by_gap_refinement": hdiv_identical,
            "reduced_a_contracting": fem_convergence["reduced_a"]["contracting"],
            "mixed_omega_contracting":
                fem_convergence["mixed_total_reduced_omega"]["contracting"],
        },
    }


def _parse_result(value: str) -> tuple[str, Path]:
    name, _, raw = value.partition("=")
    if not name or not raw:
        raise argparse.ArgumentTypeError("result must be NAME=PATH")
    return name.strip(), Path(raw.strip())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gap-family", type=Path, required=True)
    parser.add_argument("--result", action="append", type=_parse_result, required=True,
                        help="NAME=PATH, one per gap level, e.g. n06=runs/n06.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hdiv-identity-tolerance", type=float, default=1.0e-9)
    options = parser.parse_args()
    manifest = json.loads(options.gap_family.read_text(encoding="utf-8"))
    result_paths = {name: path.resolve() for name, path in options.result}
    results = {name: load_level(path) for name, path in result_paths.items()}
    report = analyze(manifest, results, result_paths,
                     hdiv_identity_tolerance=options.hdiv_identity_tolerance)
    report["gap_family"] = str(options.gap_family.resolve())
    report["gap_family_sha256"] = sha256(options.gap_family)
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"event": "gap_family_report", "passed": report["passed"],
                      "output": str(options.output)}, sort_keys=True))
    if not report["passed"]:
        raise RuntimeError(f"C-type gap family gate failed; see {options.output}")


if __name__ == "__main__":
    main()
