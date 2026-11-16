"""Derive adoption-result.json from the saved fin verdict; never solves anything.

Measurement (already executed, recorded source 76fb6b3a2):
  python validation_test/induction_heating/axisym_volume_peec_fin.py \
      --ref-sequence-json seq-o1-f1.5.json seq-o1-f3.json seq-o2-f1.5.json seq-o2-f3.json \
      --threads 8 --output verdict.json
The four reference sequences are adaptive axisymmetric A_phi FE runs written by the same
script.  This entry reads verdict.json and mesh_hashes.json (the sha256 of each PEEC
meridian mesh in the tests fixture encoding) and writes the hash-checked summary. It fails
closed: a physical FAIL, a dirty source, an open reference flag or certification item, or a
reference sequence whose final level is not a converged checkpoint writes nothing.

usage: python entry.py verdict.json mesh_hashes.json adoption-result.json
"""
import hashlib
import json
import sys
from pathlib import Path

QUANTITIES = (("R", "R_ohm", "R_ohm"), ("partial_L", "L_H", "L_H"),
              ("root_heat", "root_heat_W_per_A2", "root_heat_W_per_A2"),
              ("tip_heat", "tip_heat_W_per_A2", "tip_heat_W_per_A2"))


REFERENCE_ORDER, REFERENCE_EXTERIOR = 2, 3.0


def certified_reference(verdict):
    """Terminal checkpoint of the converged order 2, exterior 3 sequence; raise on anything else.

    The verdict must pass physically, from a clean source, with every reference flag set and every certification
    list empty. The reference rows must hold the levels 0..n-1 once each with growing DOFs, and the last level must
    be a checkpoint.
    """
    status = verdict["verdict"]["status"] if isinstance(verdict["verdict"], dict) else verdict["verdict"]
    if status != "PASS":
        raise SystemExit(f"physical verdict {status!r}: no adoption result")
    if verdict["git_dirty"]:
        raise SystemExit("dirty source: no adoption result")
    for flag in ("reference_resolved", "reference_converged", "reference_swept", "reference_loss_consistent"):
        if verdict[flag] is not True:
            raise SystemExit(f"{flag} is {verdict[flag]!r}: no adoption result")
    if verdict["reference_sequences_without_checkpoint"]:
        raise SystemExit("a reference sequence has no checkpoint: no adoption result")
    open_items = {k: v for k, v in verdict["reference_certification"].items() if v}
    if open_items or not verdict["reference_certification"]:
        raise SystemExit(f"reference certification incomplete: {open_items}")
    sequences = [s for s in verdict["reference_sequences"]
                 if s["order"] == REFERENCE_ORDER and s["exterior_factor"] == REFERENCE_EXTERIOR]
    if len(sequences) != 1 or sequences[0]["converged"] is not True or sequences[0]["stop"] != "converged":
        raise SystemExit(f"order 2 exterior 3 sequence not uniquely converged: {sequences}")
    rows = [r for r in verdict["reference_levels"]
            if r["order"] == REFERENCE_ORDER and r["exterior_factor"] == REFERENCE_EXTERIOR]
    levels = [r["level"] for r in rows]
    if sorted(levels) != list(range(len(rows))):
        raise SystemExit(f"order 2 exterior 3 levels are not 0..n-1 once each: {levels}")
    rows = sorted(rows, key=lambda r: r["level"])
    if any(b["ndof"] <= a["ndof"] for a, b in zip(rows, rows[1:])):
        raise SystemExit("order 2 exterior 3 DOFs do not grow with the level")
    final = rows[-1]
    if final["checkpoint"] is not True or final["within_budget"] is not True:
        raise SystemExit(f"final level {final['level']} is not a checkpoint inside the budget")
    return final


def main(verdict_path, hashes_path, out_path):
    raw = Path(verdict_path).read_bytes()
    verdict = json.loads(raw)
    reference = certified_reference(verdict)
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
                   reference=dict(order=2, exterior_factor=3.0, ndof=reference["ndof"], level=reference["level"],
                                  **{k: reference[k] for _, _, k in QUANTITIES}),
                   reference_spread_max=verdict["reference_spread_max"],
                   reference_resolved=verdict["reference_resolved"],
                   reference_certification=verdict["reference_certification"],
                   peec_levels=values, physical_verdict=verdict["verdict"],
                   not_claimed=verdict["not_claimed"])
    Path(out_path).write_text(json.dumps(summary, indent=1, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main(*sys.argv[1:4])
