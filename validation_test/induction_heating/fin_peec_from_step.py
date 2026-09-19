"""STEP -> fin-graded surface PEEC -> terminal and per-station fin metrics.

General-geometry runner (any swept single-solid conductor).  Exploratory:
the surface-branch PEEC kernel is the same experimental one used by the
straight fixture comparison; no workpiece coupling is applied here.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))


def run(step_path: Path, *, n_lanes=64, n_stations=20, n_outline=2048,
        lane_grading="auto", tip_lanes=8, frequency=150_000.0, sigma=5.8e7,
        cad_units_per_meter=1.0, route="auto", probe_points=41):
    from radia.fin_section import feature_panel_weights
    from radia.fin_sweep import fin_graph_from_step, probe_points_3d
    from radia.peec_fin_topology import assemble_experimental_fin_peec
    from radia.peec_proximity import _biot_savart_H

    mu0 = 4e-7 * math.pi
    omega = 2 * math.pi * frequency
    delta = math.sqrt(2 / (omega * mu0 * sigma))
    rs = 1.0 / (sigma * delta)

    sweep = fin_graph_from_step(
        step_path, n_lanes=n_lanes, n_stations=n_stations,
        n_outline=n_outline, lane_grading=lane_grading, tip_lanes=tip_lanes,
        cad_units_per_meter=cad_units_per_meter, route=route)
    graph = sweep.graph
    peec, widths = assemble_experimental_fin_peec(
        graph, sweep.rings, sigma=sigma, sheet_depth=delta)
    r_branch = np.diag(peec.R_dc) if peec.R_dc.ndim == 2 else peec.R_dc
    current = peec.compute_branch_currents(frequency, [1.0], Zs=1j * r_branch)
    r_total = float(np.real(np.vdot(current, r_branch * current)))
    l_external = float(np.real(np.vdot(current, peec.L @ current)))
    widths = np.asarray(widths)

    starts = graph.branch_xyz[:, 0, :]
    ends = graph.branch_xyz[:, 1, :]
    per_station = []
    for feat, st in zip(sweep.features, sweep.stations):
        s = feat.index
        if s >= graph.n_stations - 1:
            break
        lo, hi = s * graph.n_lanes, (s + 1) * graph.n_lanes
        i_lane = current[lo:hi]
        w_lane = widths[lo:hi]
        k_lane = i_lane / w_lane
        loss = rs * np.abs(i_lane) ** 2 / w_lane            # W/m per lane
        entry = {"station": s, "fin": feat.fin is not None,
                 "graded": bool(feat.graded),
                 "lane_current_sum_A": float(abs(np.sum(i_lane))),
                 "lateral_loss_W_per_m": float(np.sum(loss)),
                 "max_absK_over_mean": float(
                     np.max(np.abs(k_lane)) / (1.0 / np.sum(w_lane)))}
        fin = feat.fin
        if fin is not None:
            wb = feature_panel_weights(feat.analysis, fin, feat.lane_uv)
            wt = feature_panel_weights(feat.analysis, fin, feat.lane_uv, tip=True)
            lane_s = feat.analysis.arclength_of(feat.lane_uv)
            axial = fin.project(feat.lane_uv)
            offset = max(5.0 * delta, 4.0 * fin.tip_radius)
            probe = probe_points_3d(feat, st, offset, n=probe_points)
            h = _biot_savart_H(probe, starts, ends, current)
            h_abs = np.sqrt(np.sum(np.abs(h) ** 2, axis=1))
            entry.update({
                "beak_current_fraction": float(abs(np.sum(i_lane * wb)) /
                                               abs(np.sum(i_lane))),
                "beak_loss_fraction": float(np.sum(loss * wb) / np.sum(loss)),
                "tip_loss_fraction": float(np.sum(loss * wt) / np.sum(loss)),
                "beak_current_centroid_axial_m": float(
                    np.sum(np.abs(i_lane) * wb * axial) /
                    np.sum(np.abs(i_lane) * wb)),
                "tip_lanes": tip_lanes,
                "tip_arc_lane_centers": int(np.count_nonzero(
                    fin.on_tip_arc(lane_s, feat.analysis.perimeter))),
                "tip_overlap_majority_panels": int(np.count_nonzero(wt > 0.5)),
                "fin_length_m": fin.length, "tip_radius_m": fin.tip_radius,
                "tip_radius_over_delta": fin.tip_radius / delta,
                "probe_offset_m": offset,
                "probe_center_H_abs_A_per_m": float(h_abs[len(h_abs) // 2]),
                "probe_H_abs_A_per_m": h_abs.tolist(),
                "probe_xyz_m": probe.tolist(),
            })
        per_station.append(entry)

    return {
        "step_path": str(step_path), "frequency_hz": frequency,
        "sigma_S_per_m": sigma, "skin_depth_m": delta,
        "sweep": sweep.meta,
        "peec": {"R_ohm": r_total, "L_external_H": l_external,
                 "n_branches": len(graph.branches),
                 "n_nodes": len(graph.nodes)},
        "stations": per_station,
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
