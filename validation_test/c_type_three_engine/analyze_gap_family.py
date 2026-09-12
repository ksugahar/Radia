"""Read the three-engine comparison across the C-type GAP family.

The separately certified gap family refines ONLY
the gap air: N parallel slabs through the gap and an in-gap size of h/N,
with the iron, air and Kelvin sizes held at one base level.  It is not a
scale family and must not be fed to ``run_mesh_convergence.py``:

* HDiv-MMM uses the iron mesh alone, so it is the SAME at every level.  Here
  that is a gate, not a convergence quantity: the air-mesh-free route must
  not move when only the gap air changes.
* The FEM routes' increments between levels measure gap resolution only.
  Their observed order is in the layer count N. Richardson extrapolation is
  conditional on asymptotic convergence, not a rigorous error bound.

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


_PER_LEVEL_KEYS = ("kelvin_domain_vol_sha256", "observation_points")


def _fem_contract(payload: dict) -> dict:
    """What must agree between levels for the family to be one experiment.

    The Kelvin mesh hash changes with the gap; the iron mesh hash must not.
    Everything else -- mode, FEM order, the reduced-A linear
    solver, tolerances, the comparison contract, the Radia version -- must be
    identical, or the levels are not one experiment.
    """
    contracts = payload.get("engine_checkpoint_contracts") or {}
    if set(contracts) != set(ENGINES) or any(
        not contract.get("implementation_sha256")
        for contract in contracts.values()
    ) or not contracts["hdiv_mmm"].get("iron_vol_sha256"):
        raise RuntimeError("all engine implementation and iron mesh hashes are required")
    checkpoints = {
        engine: {key: value for key, value in contract.items()
                 if key not in _PER_LEVEL_KEYS}
        for engine, contract in contracts.items()
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
    if not np.isfinite(hdiv_identity_tolerance) or hdiv_identity_tolerance < 0:
        raise ValueError("hdiv_identity_tolerance must be finite and nonnegative")
    names = [str(row["name"]) for row in manifest["levels"]]
    if len(names) != len(set(names)):
        raise RuntimeError("gap level names must be unique")
    if len(names) < 3:
        raise RuntimeError("at least three gap levels are required")
    missing = [name for name in names if name not in results]
    if missing:
        raise RuntimeError(f"no result for gap levels {missing}")
    payloads = [results[name] for name in names]
    for payload in payloads:
        if (payload.get("schema") != LEVEL_SCHEMA
                or set(payload.get("engines", {})) != set(ENGINES)):
            raise RuntimeError("each gap level must contain all three formulations")
        if payload.get("mode") not in ("linear", "nonlinear"):
            raise RuntimeError("unknown gap level mode")
        if (payload["mode"] == "nonlinear"
                and not payload.get("nonlinear_converged", False)):
            raise RuntimeError("nonlinear formulation did not converge")

    points = payloads[0]["observation_points_m"]
    point_array = np.asarray(points, dtype=float)
    if (point_array.ndim != 2 or point_array.shape[1] != 3
            or not len(point_array) or not np.isfinite(point_array).all()):
        raise RuntimeError("observation points must be a nonempty finite (N,3) array")
    for payload in payloads[1:]:
        if (payload["observation_points_m"] != points
                or payload["gap_core_half_length_m"]
                != payloads[0]["gap_core_half_length_m"]):
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
    layers = np.asarray([row["gap_layers"] for row in manifest["levels"]], dtype=float)
    if (not np.isfinite(ratio) or ratio <= 1 or not np.isfinite(layers).all()
            or np.any(layers <= 0) or np.any(layers != np.floor(layers))
            or not np.allclose(layers[1:] / layers[:-1], ratio, rtol=1e-12, atol=0)):
        raise RuntimeError("gap layers must follow the stated increasing refinement ratio")
    half_length = float(payloads[0]["gap_core_half_length_m"])
    if not np.isfinite(half_length) or half_length <= 0:
        raise RuntimeError("gap core half length must be positive and finite")
    selector = convergence._selector(payloads[-1])
    if not np.any(selector):
        raise RuntimeError("no observation points in the gap core")
    for payload in payloads:
        for engine in ENGINES:
            field = np.asarray(
                payload["median_plane_projected_fields_T"][engine], dtype=float)
            if field.shape != point_array.shape or not np.isfinite(field).all():
                raise RuntimeError("all fields must be finite and match the observation points")
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
    for index, (name, row, payload) in enumerate(zip(names, manifest["levels"], payloads)):
        pairs = {}
        for left_index, left in enumerate(ENGINES):
            for right in ENGINES[left_index + 1:]:
                a, b = fields[left][index], fields[right][index]
                pairs[f"{left}__vs__{right}"] = {
                    "relative_rms": convergence._relative_rms(a, b),
                    "maximum_absolute_difference_T": float(np.max(np.linalg.norm(a - b, axis=1))),
                }
        levels.append({
            "name": name,
            "gap_layers": int(row["gap_layers"]),
            "gap_elements": int(row["gap_elements"]),
            "result": str(result_paths[name]),
            "result_sha256": sha256(result_paths[name]),
            "machine": payload.get("machine"),
            "maximum_gap_core_pairwise_relative_rms":
                max(pair["relative_rms"] for pair in pairs.values()),
            "pairwise_median_projected_gap_core":
                pairs,
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
                          "conditional Richardson estimates describe refinement "
                          "sensitivity, not a rigorous error bound; "
                          "HDiv-MMM is the fixed air-mesh-free "
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
    name, raw = name.strip(), raw.strip()
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
    names = [name for name, _ in options.result]
    if len(names) != len(set(names)):
        parser.error("duplicate --result level name")
    manifest = json.loads(options.gap_family.read_text(encoding="utf-8"))
    result_paths = {name: path.resolve() for name, path in options.result}
    results = {name: load_level(path) for name, path in result_paths.items()}
    report = analyze(manifest, results, result_paths,
                     hdiv_identity_tolerance=options.hdiv_identity_tolerance)
    report["gap_family"] = str(options.gap_family.resolve())
    report["gap_family_sha256"] = sha256(options.gap_family)
    serialized = json.dumps(report, indent=2, allow_nan=False) + "\n"
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(serialized, encoding="utf-8")
    print(json.dumps({"event": "gap_family_report", "passed": report["passed"],
                      "output": str(options.output)}, sort_keys=True))
    if not report["passed"]:
        raise RuntimeError(f"C-type gap family gate failed; see {options.output}")


if __name__ == "__main__":
    main()
