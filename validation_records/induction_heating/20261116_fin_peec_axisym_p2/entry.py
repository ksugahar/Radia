"""Derive adoption-result.json from the saved fin verdict; never solves anything.

Measurement (already executed, recorded source 76fb6b3a2):
  python validation_test/induction_heating/axisym_volume_peec_fin.py \
      --ref-sequence-json seq-o1-f1.5.json seq-o1-f3.json seq-o2-f1.5.json seq-o2-f3.json \
      --threads 8 --output verdict.json
The four reference sequences are adaptive axisymmetric A_phi FE runs written by the same
script.  This entry reads verdict.json and mesh_hashes.json (the sha256 of each PEEC
meridian mesh in the tests fixture encoding) and writes the hash-checked summary. It fails
closed: it rebuilds every reference sequence from the verdict, certifies it and resolves the
reference with the runner's own functions under the approved settings, requires the stored
checks, certification and physical verdict to equal their re-derivation, and writes nothing
otherwise.

usage: python entry.py verdict.json mesh_hashes.json adoption-result.json
"""
import hashlib
import json
import sys
from pathlib import Path

QUANTITIES = (("R", "R_ohm", "R_ohm"), ("partial_L", "L_H", "L_H"),
              ("root_heat", "root_heat_W_per_A2", "root_heat_W_per_A2"),
              ("tip_heat", "tip_heat_W_per_A2", "tip_heat_W_per_A2"))


RUNNER_DIR = Path(__file__).resolve().parents[3] / "validation_test" / "induction_heating"
# The configuration the user approved on 2026-11-16 (the runner defaults of the r08 verdict).
ORDERS, EXTERIORS = (1, 2), (1.5, 3.0)
SETTINGS = dict(theta=0.5, max_ndof=6e6, max_levels=16, budget_s=7200.0,
                per_order={"1": dict(start_maxh_mm=0.05, air_maxh_mm=20.0, air_divisor=20.0),
                           "2": dict(start_maxh_mm=0.05, air_maxh_mm=20.0, air_divisor=20.0)})


def _runner():
    sys.path.insert(0, str(RUNNER_DIR))
    try:
        import axisym_volume_peec_fin as runner
        from fin_return_principle import PHYSICAL_TARGETS, physical_verdict
    finally:
        sys.path.remove(str(RUNNER_DIR))
    return runner, PHYSICAL_TARGETS, physical_verdict


def certified_reference(verdict):
    """Re-derive the reference from the recorded sequences with the runner's own rules; raise on any mismatch.

    Every configured sequence is rebuilt from its summary and its level rows and certified by the runner's
    certify_sequence under the approved settings; the runner's resolve_reference then recomputes the checks, the
    spread and the resolution, which must equal the stored ones, and the physical verdict is recomputed from the
    finest PEEC level against the re-derived reference. Returns the reference level (final checkpoint of order 2
    at exterior 3).
    """
    runner, targets, physical_verdict = _runner()
    if verdict["git_dirty"] is not False:
        raise SystemExit("dirty source: no adoption result")
    tol = {k: 0.1*targets[k] for k in ("R", "partial_L", "root_heat", "tip_heat")}
    stored = verdict["reference_settings"]
    if stored != dict(SETTINGS, tolerance=tol):
        raise SystemExit(f"reference settings {stored} differ from the approved configuration")

    def expected_of(order):
        per = SETTINGS["per_order"][str(order)]
        return dict(per, theta=SETTINGS["theta"], checkpoint_dof_growth=2.0, tolerance=tol,
                    max_ndof=SETTINGS["max_ndof"], max_levels=SETTINGS["max_levels"], budget_s=SETTINGS["budget_s"])

    sequences = []
    for meta in verdict["reference_sequences"]:
        rows = [{k: v for k, v in r.items() if k != "order"} for r in verdict["reference_levels"]
                if r["order"] == meta["order"] and r["exterior_factor"] == meta["exterior_factor"]]
        sequences.append(dict(meta, levels=rows))
    if sum(len(s["levels"]) for s in sequences) != len(verdict["reference_levels"]):
        raise SystemExit("reference rows outside the recorded sequences")
    try:
        res = runner.resolve_reference(sequences, ORDERS, EXTERIORS, expected_of, targets)
    except ValueError as exc:
        raise SystemExit(f"reference sequences refused: {exc}")
    problems = {k: v for k, v in res["certification"].items() if v}
    if problems or res["missing"] or not res["resolved"]:
        raise SystemExit(f"reference not certified: {problems or res['missing'] or 'unresolved'}")
    recomputed = dict(reference_certification=res["certification"], reference_converged=res["converged"],
                      reference_swept=res["swept"], reference_sequences_without_checkpoint=res["missing"],
                      reference_checks=res["checks"], reference_spread_max=res["spread"],
                      reference_loss_consistent=res["consistent"], reference_resolved=res["resolved"])
    for key, value in recomputed.items():
        if verdict[key] != value:
            raise SystemExit(f"stored {key} differs from its re-derivation")
    finest = verdict["peec_levels"][-1]
    candidate = runner.quantities(finest)
    physical = physical_verdict(candidate, res["reference"], reference_resolved=True)
    if physical["status"] != "PASS" or verdict["verdict"] != physical:
        raise SystemExit(f"physical verdict {physical['status']!r} (stored {verdict['verdict'].get('status')!r})")
    return res["reference_level"]


def main(verdict_path, hashes_path, out_path):
    raw = Path(verdict_path).read_bytes()
    verdict = json.loads(raw)
    reference = certified_reference(verdict)
    reference_order = 2
    source = verdict["git_head"]
    hashes = json.loads(Path(hashes_path).read_text())
    levels, values = [], []
    for row in verdict["peec_levels"]:
        h = row["boundary_h_mm"]
        mesh = hashes[str(h)]
        if mesh["cells"] != row["cells"]:
            raise SystemExit(f"mesh cells {mesh['cells']} != recorded {row['cells']} at h={h}")
        values.append(dict(boundary_h_mm=h, cells=row["cells"], **{k: row[c] for _, c, k in QUANTITIES}))
        for index, (name, column, ref_column) in enumerate(QUANTITIES, start=1):
            error = abs(row[column]/reference[ref_column] - 1.0)
            levels.append(dict(id=f"{name}-h{h}", parameters=dict(quantity_index=index, boundary_h_mm=h),
                               error=error, relative_residual=row["true_relative_residual"], source_sha=source,
                               mesh_sha256=mesh["sha256"], runtime_s=row["assembly_s"] + row["solve_s"]))
    if any(not 0.0 <= lv["error"] < float("inf") for lv in levels) or len(values) != 3:
        raise SystemExit("PEEC levels incomplete or nonfinite: no adoption result")
    result = dict(source_sha=source, verdict="PASS", levels=levels)
    summary = dict(schema="radia.ih.fin-peec-axisym-p2-adoption.v1", adoption_result=result,
                   original_verdict_sha256=hashlib.sha256(raw).hexdigest(),
                   verdict_schema=verdict["schema"], script_sha256=verdict["script_sha256"],
                   peec_matrices_sha256=verdict["peec_runtime"]["peec_matrices_sha256"],
                   axifem_sha256=verdict["reference_runtime"]["axifem_sha256"],
                   conditions=dict(radius_mm=verdict["radius_mm"], frequency_hz=verdict["frequency_hz"],
                                   sigma_S_per_m=verdict["sigma_S_per_m"], skin_depth_mm=verdict["skin_depth_mm"],
                                   cuts_mm=verdict["cuts_mm"], normalization=verdict["normalization"]),
                   reference=dict(order=reference_order, exterior_factor=reference["exterior_factor"], ndof=reference["ndof"], level=reference["level"],
                                  **{k: reference[k] for _, _, k in QUANTITIES}),
                   reference_spread_max=verdict["reference_spread_max"],
                   reference_resolved=verdict["reference_resolved"],
                   reference_certification=verdict["reference_certification"],
                   peec_levels=values, physical_verdict=verdict["verdict"],
                   not_claimed=verdict["not_claimed"])
    Path(out_path).write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main(*sys.argv[1:4])
