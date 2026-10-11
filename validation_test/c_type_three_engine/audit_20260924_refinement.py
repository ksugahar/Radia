"""Recompute historical refinement metrics from immutable Git evidence.

This reads JSON only. It never imports or runs the archived solver drivers.
Use a full-history checkout; this is not current-runtime release acceptance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess


ARCHIVE = "75e49aa58b43558c7a4cd3f48b6f24d6ca169d96"
ROOT = Path(__file__).resolve().parents[2]
PREFIX = "validation_test/c_type_three_engine/results/compute-host_20260924_"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def norm(field):
    require(all(math.isfinite(x) for row in field for x in row), "nonfinite field")
    return math.sqrt(math.fsum(x * x for row in field for x in row))


def difference(left, right):
    require(len(left) == len(right), "field length mismatch")
    require(all(len(row) == 3 for row in left + right), "expected vector fields")
    return norm([[a - b for a, b in zip(x, y)] for x, y in zip(left, right)])


def close(actual, expected):
    require(math.isclose(actual, expected, rel_tol=2e-10, abs_tol=1e-14),
            f"archived metric mismatch: {actual} != {expected}")


def project(points, field):
    require(len(points) == len(field), "point/field length mismatch")
    keys = {tuple(round(x, 12) for x in p): i for i, p in enumerate(points)}
    require(len(keys) == len(points), "duplicate observation points")
    result = []
    for p, row in zip(points, field):
        index = keys[(round(p[0], 12), round(p[1], 12), round(-p[2], 12))]
        result.append([(a + s * b) / 2
                       for a, b, s in zip(row, field[index], (-1, -1, 1))])
    return result


def summarize():
    inputs = {}

    def read(name):
        path = PREFIX + name + ".json"
        raw = subprocess.check_output(["git", "show", f"{ARCHIVE}:{path}"], cwd=ROOT)
        inputs[name] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    names = ("cached", "medium", "fine", "fullfiner", "order3")
    campaigns = {name: read(name + "/cached-three") for name in names}
    base = campaigns["cached"]
    points = base["observation_points_m"]
    mask = [abs(p[0]) <= base["gap_core_half_length_m"] + 1e-12 for p in points]

    def core(field):
        return [row for row, selected in zip(field, mask) if selected]

    summaries = {}
    for name, data in campaigns.items():
        require(data["passed"] and data["mode"] == "linear", "historical campaign failed")
        for key in ("observation_points_m", "coil", "gap_core_half_length_m"):
            require(data[key] == base[key], f"changed comparison contract: {key}")
        require(data["comparison_contract"]["cad_authority_sha256"] ==
                base["comparison_contract"]["cad_authority_sha256"], "changed CAD authority")
        projected = data["median_plane_projected_fields_T"]
        for engine, field in data["fields_T"].items():
            close(difference(project(points, field), projected[engine]), 0.0)
        pairs = {}
        expected_pairs = {"hdiv_mmm__vs__reduced_a", "hdiv_mmm__vs__mixed_total_reduced_omega",
                          "reduced_a__vs__mixed_total_reduced_omega"}
        require(set(data["pairwise_median_projected_gap_core"]) == expected_pairs, "missing engine pair")
        for key in sorted(expected_pairs):
            left, right = key.split("__vs__")
            a, b = core(projected[left]), core(projected[right])
            metric = difference(a, b) / norm(a)
            close(metric, data["pairwise_median_projected_gap_core"][key]["relative_rms"])
            pairs[key] = metric
        close(max(pairs.values()), data["maximum_gap_core_pairwise_relative_rms"])
        summaries[name] = {
            "historical_projected_core_gate_passed": max(pairs.values()) <= data["relative_rms_tolerance"],
            "gate": data["relative_rms_tolerance"], "core_pairwise_relative_rms": pairs,
            "engines": {key: {"primary_dofs": value["ndof"], "runtime_s": value["runtime_s"]}
                        for key, value in data["engines"].items()},
        }

    fem = {}
    saved = read("fullfiner/fem_refinement")
    for engine in ("reduced_a", "mixed_total_reduced_omega"):
        levels = [campaigns[key] for key in ("medium", "fine", "fullfiner")]
        contracts = [d["engine_checkpoint_contracts"][engine] for d in levels]
        for key in ("mode", "fem_order", "implementation_sha256"):
            require(all(c[key] == contracts[0][key] for c in contracts), f"changed {engine} {key}")
        fields = [core(d["median_plane_projected_fields_T"][engine]) for d in levels]
        increments = [difference(a, b) / norm(fields[-1]) for a, b in zip(fields, fields[1:])]
        for actual, expected in zip(increments, saved["engines"][engine]["core_increments"]):
            close(actual, expected)
        fem[engine] = {"core_increments": increments, "contraction_ratio": increments[1] / increments[0]}

    bdm1_saved = read("fine/three_level_increments")["engines"]["hdiv_mmm"]
    bdm1_fields = [core(campaigns[key]["median_plane_projected_fields_T"]["hdiv_mmm"])
                   for key in ("cached", "medium", "fine")]
    bdm1_increments = [difference(a, b) / norm(bdm1_fields[-1])
                       for a, b in zip(bdm1_fields, bdm1_fields[1:])]
    for actual, expected in zip(bdm1_increments, bdm1_saved["core_increments_normalized_by_fine"]):
        close(actual, expected)
    bdm1_ratio = bdm1_increments[1] / bdm1_increments[0]
    close(bdm1_ratio, bdm1_saved["core_increment_contraction_ratio"])
    bdm1 = {"core_increments": bdm1_increments, "contraction_ratio": bdm1_ratio,
            "status": "not_contracting", "cause": "unresolved",
            "interpretation": "BDM1 refinement did not demonstrate convergence; not an accuracy certificate."}

    bdm_rows = [next(row for row in read(name)["rows"] if row["rule"] == 2) for name in
                ("bdm2levels/medium_bdm2", "bdm2levels/fine_bdm2", "bdm2finer/finer_bdm2")]
    require(all(row["solve_scalars"]["last_solve_converged"] == 1 for row in bdm_rows),
            "BDM2 solve did not converge")
    fields = [core(project(points, row["field_T"])) for row in bdm_rows]
    increments = [difference(a, b) / norm(fields[-1]) for a, b in zip(fields, fields[1:])]
    saved = read("bdm2finer/increments")
    for actual, expected in zip(increments, saved["core_increments"]):
        close(actual, expected)
    bdm = {"core_increments": increments, "contraction_ratio": increments[1] / increments[0],
           "native_residuals": [row["solve_scalars"]["last_solve_final_relative_residual"] for row in bdm_rows]}

    order_sensitivity = {}
    saved = read("order3/order_sensitivity")
    for engine in ("reduced_a", "mixed_total_reduced_omega"):
        current = core(campaigns["order3"]["median_plane_projected_fields_T"][engine])
        order_sensitivity[engine] = {}
        for label, name in (("coarse_p2", "cached"), ("finer_p2", "fullfiner")):
            reference = core(campaigns[name]["median_plane_projected_fields_T"][engine])
            metric = difference(current, reference) / norm(reference)
            close(metric, saved["comparisons"][engine][label])
            order_sensitivity[engine][label] = metric

    sensitivity = {}
    for name in ("replay/replay", "dense/dense_comparison", "farrule/far_rule_comparison",
                 "nearquad/near_quad_comparison", "sourcequad/source_quadrature_comparison"):
        rows = read(name)["rows"]
        reference = rows[0]["field_T"]
        sensitivity[name.split("/")[0]] = [{
            "relative_vector_difference_from_first_row": difference(row["field_T"], reference) / norm(reference),
            "native_converged": bool(row["solve_scalars"]["last_solve_converged"]),
            "native_residual": row["solve_scalars"]["last_solve_final_relative_residual"],
        } for row in rows]

    vacuum = read("vacuum/vacuum")
    require(vacuum["points_m"] == points, "vacuum stencil changed")
    vacuum_error = difference(vacuum["computed_B_T"], vacuum["target_B_T"]) / norm(vacuum["target_B_T"])
    close(vacuum_error, vacuum["raw_all_relative_rms"])
    provenance = read("order3/cached-three-provenance")
    return {
        "schema": "radia.historical-c-type-refinement-audit.v1",
        "status": "historical_metrics_recomputed_not_current_release_acceptance",
        "archive_commit": ARCHIVE,
        "runtime": {"versions": provenance["versions"], "threads": 8},
        "scope": "Linear mu_r=1000, parity-projected gap core; vacuum diagnostic is separate.",
        "limitations": ["NGSolve 2606 evidence, not a 2607 rerun.",
                        "Pairwise agreement and contracting increments are not absolute error bounds.",
                        "No energy/coenergy acceptance or matched-error performance claim.",
                        "Runtime includes exact source-cache preparation; process RSS is not isolated engine memory.",
                        "These historic volume-load results do not validate current auto/surface-flux defaults."],
        "campaigns": summaries, "fem_mesh_sensitivity": fem, "bdm2_mesh_sensitivity": bdm,
        "bdm1_mesh_sensitivity": bdm1,
        "fem_order_sensitivity": order_sensitivity,
        "hdiv_controlled_sensitivities": sensitivity, "vacuum_raw_relative_rms": vacuum_error,
        "input_json_sha256": inputs,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = summarize()
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Recomputed historical metrics from {len(result['input_json_sha256'])} JSON records.")
