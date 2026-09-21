"""How fine does the BEM-A reference have to be, and what does it cost?

With the fin PEEC at its converged discretisation the comparison is entirely
BEM-A: 437 s of a roughly 8 minute run at `maxh = 0.75 mm`.  The compressed
`hacapk_cocr` solver does not help at this size -- it is 5% slower than dense
COCR on the same system, for the same answer -- so the only remaining lever is
the reference mesh itself.

The existing refinement evidence only went finer than 0.75 mm, on the 24 mm
fixture.  This goes coarser on the accepted 48 mm one and reports, per level,
the reference R and whether the delivery gate still closes against the PEEC.
A coarser reference is only useful if it still decides the gate the same way.
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

from compare_beak_fin_bema import compare_peec, solve_bema  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, default=(
        ROOT / "tests/coil_from_cad/fixtures/beak_fin_48mm.step"))
    parser.add_argument("--frequency", type=float, default=150_000.0)
    parser.add_argument("--maxh-mm", type=float, nargs="+",
                        default=[2.0, 1.5, 1.0, 0.75])
    parser.add_argument("--n-peri", type=int, default=256)
    parser.add_argument("--n-stations", type=int, default=5)
    parser.add_argument("--curve-order", type=int, default=2)
    parser.add_argument("--fes-order", type=int, default=0)
    parser.add_argument("--bema-solver", default="cocr")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    rows = []
    for maxh_mm in args.maxh_mm:
        t0 = time.perf_counter()
        state = solve_bema(args.step, frequency=args.frequency,
                           maxh=maxh_mm * 1e-3, bema_solver=args.bema_solver,
                           curve_order=args.curve_order,
                           fes_order=args.fes_order)
        bema_seconds = time.perf_counter() - t0
        t1 = time.perf_counter()
        result = compare_peec(state, args.step, n_peri=args.n_peri,
                              n_stations=args.n_stations,
                              experimental_peec=True)
        peec_seconds = time.perf_counter() - t1
        profile = result["surface_current_profile"]
        errors = profile["integral_metric_errors"]
        row = {
            "maxh_mm": maxh_mm,
            "n_surface_faces": state["bema"]["n_surface_faces"],
            "n_J": state["bema"]["n_J"],
            "bema_seconds": bema_seconds,
            "peec_seconds": peec_seconds,
            "bema_R_uohm": state["bema"]["R_ohm"] * 1e6,
            "bema_L_nH": state["bema"]["L_H"] * 1e9,
            "peec_R_uohm": result["experimental_peec"]["R_ohm"] * 1e6,
            "errors": {
                "beak_current_fraction_relative":
                    errors["beak_current_fraction_relative"],
                "beak_loss_fraction_relative":
                    errors["beak_loss_fraction_relative"],
                "beak_current_centroid_absolute_mm":
                    errors["beak_current_centroid_absolute_m"] * 1e3,
                "probe_H_relative_L2": profile["probe_H_relative_L2"],
            },
            "accepted": profile["accepted"],
        }
        rows.append(row)
        print(json.dumps({"phase": "level", "maxh_mm": maxh_mm,
                          "faces": row["n_surface_faces"],
                          "bema_s": round(bema_seconds, 1),
                          "peec_s": round(peec_seconds, 2),
                          "accepted": row["accepted"]}), flush=True)

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    finest = rows[-1]
    report = {
        "schema": "radia.beak_fin_reference_mesh_cost.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "host": platform.node(),
        "question": ("with the PEEC side converged and cheap, the comparison "
                     "is the BEM-A reference.  How coarse can that reference "
                     "be and still decide the delivery gate the same way?"),
        "fixture": str(args.step),
        "peec": {"n_peri": args.n_peri, "n_stations": args.n_stations},
        "solver_note": ("hacapk_cocr was measured on the 0.75 mm level at "
                        "459.8 s against 437.3 s for dense cocr, agreeing to "
                        "7.5e-7 in R: the H-matrix does not pay at this size, "
                        "so the mesh is the only lever"),
        "levels": rows,
        "reference_R_uohm_at_finest": finest["bema_R_uohm"],
        "relative_R_to_finest": [
            r["bema_R_uohm"] / finest["bema_R_uohm"] - 1.0 for r in rows],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({"phase": "complete", "output": str(args.output),
                      "accepted": [r["accepted"] for r in rows],
                      "bema_seconds": [round(r["bema_seconds"], 1)
                                       for r in rows]}))


if __name__ == "__main__":
    main()
