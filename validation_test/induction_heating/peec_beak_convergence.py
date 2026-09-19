"""Run one experimental 3-D beak-fin PEEC convergence point."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


def run(step_path: Path, *, n_peri: int, max_axial_m: float,
        frequency=150_000.0, sigma=5.8e7, lane_grading="uniform"):
    from build123d import import_step

    from radia.peec_fin_topology import (
        assemble_experimental_fin_peec,
        build_hybrid_surface_topology_from_straight_prism_step,
    )

    solid = import_step(str(step_path)).solids()[0]
    bbox = solid.bounding_box()
    length_m = (float(bbox.max.Z) - float(bbox.min.Z)) / 1000.0
    n_stations = math.ceil(length_m / max_axial_m) + 1
    graph, cad = build_hybrid_surface_topology_from_straight_prism_step(
        step_path, n_peri=n_peri, n_stations=n_stations,
        lane_grading=lane_grading)
    rings = np.empty((graph.n_stations, graph.n_lanes, 3))
    for station in range(graph.n_stations - 1):
        lo = station * graph.n_lanes
        hi = lo + graph.n_lanes
        rings[station] = graph.branch_xyz[lo:hi, 0]
    rings[-1] = graph.branch_xyz[
        (graph.n_stations - 2) * graph.n_lanes:
        (graph.n_stations - 1) * graph.n_lanes, 1]
    mu0 = 4e-7 * math.pi
    omega = 2 * math.pi * frequency
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    peec, _ = assemble_experimental_fin_peec(
        graph, rings, sigma=sigma, sheet_depth=delta)
    r_branch = np.diag(peec.R_dc) if peec.R_dc.ndim == 2 else peec.R_dc
    current = peec.compute_branch_currents(
        frequency, [1.0], Zs=1j * r_branch)
    resistance = float(np.real(np.vdot(current, r_branch * current)))
    return {"length_m": length_m, "n_peri": n_peri,
            "n_stations": n_stations,
            "axial_step_m": length_m / (n_stations - 1),
            "n_branches": len(graph.branches), "R_ohm": resistance,
            "lane_grading": cad["lane_grading"],
            "fin_analysis": cad["fin_analysis"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, required=True)
    parser.add_argument("--n-peri", type=int, required=True)
    parser.add_argument("--max-axial-mm", type=float, default=1.5)
    parser.add_argument("--lane-grading", choices=("uniform", "auto"),
                        default="uniform")
    args = parser.parse_args()
    print(json.dumps(run(args.step, n_peri=args.n_peri,
                         max_axial_m=args.max_axial_mm * 1e-3,
                         lane_grading=args.lane_grading), indent=2))


if __name__ == "__main__":
    main()
