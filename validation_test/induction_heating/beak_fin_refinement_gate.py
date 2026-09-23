"""Gate 1: perimeter and axial refinement of the beak-fin PEEC.

The delivery gate closed at 48 mm on one discretisation (n_peri=256,
n_stations=33).  A single passing point is not convergence: it does not say
whether the agreement survives refinement, nor whether it was reached from
above or below.  This driver holds the BEM-A reference fixed and refines the
PEEC side in each direction independently, so a drift in the delivery metrics
can be attributed to the perimeter or to the axial discretisation and not to
both at once.

The BEM-A solve is the whole cost -- 1164 s on the 48 mm fixture -- and it does
not depend on either PEEC parameter, so it is solved once and reused.  That is
what `compare_beak_fin_bema.solve_bema` exists for.

Each level also records the gate-1 geometry rejections: station alignment,
zero-area cells, and a missing tip.

The refinement stops where the experimental PEEC does.  Its matrices are dense
in the branch count, and the accepted point -- n_peri=256 with n_stations=33,
about 8200 branches -- already holds 18.5 GB resident.  Doubling either
parameter quadruples that, so levels beyond `--max-branches` are recorded as
skipped rather than attempted; a finer sequence needs a compressed PEEC, not a
bigger machine.
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

from compare_beak_fin_bema import compare_peec, solve_bema  # noqa: E402


def _level(state, step_path, *, n_peri, n_stations, max_branches):
    # Longitudinal branches alone; the transverse ones only add to this, so
    # the estimate is a lower bound on what the dense assembly would hold.
    estimate = n_peri * (n_stations - 1)
    if estimate > max_branches:
        return {"n_peri": n_peri, "n_stations": n_stations,
                "skipped": True,
                "estimated_longitudinal_branches": int(estimate),
                "max_branches": int(max_branches),
                "reason": ("the experimental PEEC assembles dense matrices in "
                           "the branch count; this level does not fit")}
    result = compare_peec(state, step_path, n_peri=n_peri,
                          n_stations=n_stations, experimental_peec=True)
    profile = result["surface_current_profile"]
    peec = profile["integral_metrics"]["peec"]
    bema = profile["integral_metrics"]["bema"]
    errors = profile["integral_metric_errors"]
    return {
        "n_peri": n_peri,
        "n_stations": n_stations,
        "skipped": False,
        "n_branches": result["peec_topology"]["n_branches"],
        "geometry_checks": result["geometry_checks"],
        "peec_R_uohm": result["experimental_peec"]["R_ohm"] * 1e6,
        "relative_R_error": result["experimental_peec"]["relative_R_error"],
        "relative_L_error": result["experimental_peec"]["relative_L_error"],
        "peec": {
            "beak_current_fraction": peec["beak_current_fraction"],
            "beak_loss_fraction": peec["beak_loss_fraction"],
            "tip_loss_fraction": peec["tip_loss_fraction"],
            "beak_current_centroid_x_mm":
                peec["beak_current_centroid_x_m"] * 1e3,
        },
        "bema": {
            "beak_current_fraction": bema["beak_current_fraction"],
            "beak_loss_fraction": bema["beak_loss_fraction"],
            "tip_loss_fraction": bema["tip_loss_fraction"],
            "beak_current_centroid_x_mm":
                bema["beak_current_centroid_x_m"] * 1e3,
        },
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, default=(
        ROOT / "tests/coil_from_cad/fixtures/beak_fin_48mm.step"))
    parser.add_argument("--frequency", type=float, default=150_000.0)
    parser.add_argument("--maxh", type=float, default=0.00075)
    parser.add_argument("--curve-order", type=int, default=2)
    parser.add_argument("--fes-order", type=int, default=0)
    parser.add_argument("--bema-solver", default="cocr")
    parser.add_argument("--base-n-peri", type=int, default=256)
    parser.add_argument("--base-n-stations", type=int, default=33)
    parser.add_argument("--peri-levels", type=int, nargs="+",
                        default=[32, 64, 128, 256])
    parser.add_argument("--station-levels", type=int, nargs="+",
                        default=[5, 9, 17, 33])
    parser.add_argument("--max-branches", type=int, default=10000,
                        help="levels above this are recorded as skipped")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    print(json.dumps({"phase": "bema", "fixture": str(args.step),
                      "maxh_m": args.maxh}), flush=True)
    state = solve_bema(args.step, frequency=args.frequency, maxh=args.maxh,
                       bema_solver=args.bema_solver,
                       curve_order=args.curve_order,
                       fes_order=args.fes_order)
    print(json.dumps({"phase": "bema_done",
                      "R_uohm": state["bema"]["R_ohm"] * 1e6,
                      "n_J": state["bema"]["n_J"]}), flush=True)

    perimeter, axial = [], []
    for n_peri in args.peri_levels:
        row = _level(state, args.step, n_peri=n_peri,
                     n_stations=args.base_n_stations,
                     max_branches=args.max_branches)
        perimeter.append(row)
        print(json.dumps({"phase": "perimeter", **{
            k: row.get(k) for k in ("n_peri", "n_stations", "accepted", "skipped")}}),
            flush=True)
    for n_stations in args.station_levels:
        row = _level(state, args.step, n_peri=args.base_n_peri,
                     n_stations=n_stations,
                     max_branches=args.max_branches)
        axial.append(row)
        print(json.dumps({"phase": "axial", **{
            k: row.get(k) for k in ("n_peri", "n_stations", "accepted", "skipped")}}),
            flush=True)

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    report = {
        "schema": "radia.beak_fin_refinement_gate.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "host": platform.node(),
        "question": (
            "the delivery gate closed at one discretisation.  Does it survive "
            "refinement of the perimeter and of the axial stations "
            "independently, and do the geometry rejections hold at every "
            "level?"),
        "fixture": str(args.step),
        "frequency_hz": args.frequency,
        "bema_reference": {
            "maxh_m": args.maxh,
            "curve_order": args.curve_order,
            "fes_order": args.fes_order,
            "solver": args.bema_solver,
            "R_uohm": state["bema"]["R_ohm"] * 1e6,
            "L_nH": state["bema"]["L_H"] * 1e9,
            "n_J": state["bema"]["n_J"],
            "n_surface_faces": state["bema"]["n_surface_faces"],
            "note": ("solved once; it depends on the fixture, the mesh and the "
                     "basis, and on neither PEEC refinement parameter"),
        },
        "perimeter_refinement": {"n_stations": args.base_n_stations,
                                 "levels": perimeter},
        "axial_refinement": {"n_peri": args.base_n_peri, "levels": axial},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({"phase": "complete", "output": str(args.output),
                      "perimeter_accepted": [r.get("accepted") for r in perimeter],
                      "axial_accepted": [r.get("accepted") for r in axial]}))


if __name__ == "__main__":
    main()
