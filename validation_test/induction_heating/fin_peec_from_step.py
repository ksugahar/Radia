"""STEP -> fin-graded surface PEEC -> terminal and per-station fin metrics.

General-geometry runner (any swept single-solid conductor).  Exploratory:
the surface-branch PEEC kernel is the same experimental one used by the
straight fixture comparison; no workpiece coupling is applied here.  The
solve itself lives in ``radia.fin_surface_solver`` and is shared with the
``--coil-solver fin-surface`` backend of ``calc_inductance.py``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


def run(step_path: Path, *, n_lanes=64, n_stations=20, n_outline=2048,
        lane_grading="auto", tip_lanes=8, frequency=150_000.0, sigma=5.8e7,
        cad_units_per_meter=1.0, route="auto", probe_points=41):
    from radia.fin_surface_solver import solve_fin_surface

    sol = solve_fin_surface(
        step_path, frequency=frequency, sigma=sigma, n_lanes=n_lanes,
        n_stations=n_stations, n_outline=n_outline, lane_grading=lane_grading,
        tip_lanes=tip_lanes, cad_units_per_meter=cad_units_per_meter,
        route=route, probe_points=probe_points)
    return {
        "step_path": str(step_path), "frequency_hz": frequency,
        "sigma_S_per_m": sigma, "skin_depth_m": sol["skin_depth_m"],
        "sweep": sol["sweep"].meta,
        "peec": {"R_ohm": sol["R_ohm"], "L_external_H": sol["L_external_H"],
                 "n_branches": sol["n_branches"], "n_nodes": sol["n_nodes"]},
        "stations": sol["stations"],
        "status": ("exploratory surface PEEC; fin metrics per station, no "
                   "workpiece coupling, kernel not yet accepted as IH backend"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, required=True)
    parser.add_argument("--n-lanes", type=int, default=64)
    parser.add_argument("--n-stations", type=int, default=20)
    parser.add_argument("--n-outline", type=int, default=2048)
    parser.add_argument("--tip-lanes", type=int, default=8)
    parser.add_argument("--lane-grading", choices=("uniform", "auto"),
                        default="auto")
    parser.add_argument("--route", choices=("auto", "straight_prism",
                                            "section_planes"), default="auto")
    parser.add_argument("--frequency", type=float, default=150_000.0)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--cad-units-per-meter", type=float, default=1.0,
                        help="1000 for a millimetre STEP (auto-checked)")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.step, n_lanes=args.n_lanes, n_stations=args.n_stations,
                 n_outline=args.n_outline, lane_grading=args.lane_grading,
                 tip_lanes=args.tip_lanes, frequency=args.frequency,
                 sigma=args.sigma, cad_units_per_meter=args.cad_units_per_meter,
                 route=args.route)
    text = json.dumps(result, indent=2)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
