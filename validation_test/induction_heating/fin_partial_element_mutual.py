"""Gate 2: is the filamentary mutual defensible between adjacent fin bands?

The fin PEEC adds each surface band through `add_connected_segment`, which
creates one segment with no parent.  `MutualInductanceRectBar` -- the
cross-section-averaged kernel -- is reached only for sub-filaments sharing a
parent, so every branch pair in the fin graph takes the filamentary Neumann
formula instead.  Verified behaviourally: the built L matrix reproduces
Grover's equal-parallel-filament closed form to machine precision.

That is the path the kernel's own comment warns about, in its words, for the
"spurious circulating current artifact that filamentary Neumann gives for
close parallel bars".  On this fin the bands are exactly that case: at
n_peri=256 the perimeter spacing is 0.0894 mm, the band width is the same
0.0894 mm, and the sheet depth is the 0.1706 mm skin depth, so neighbours
touch and their separation is smaller than their own depth.

This driver quantifies the error: for each sampled neighbour pair it evaluates
the filamentary mutual the kernel uses and the cross-section-averaged mutual
over the same two bars, with the bands oriented on the actual surface.  It
also reports how much the averaged value depends on that orientation, which
matters only if the averaged path is adopted.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

from fin_partial_element_geometry import (  # noqa: E402
    MU0, bar_mutual, filament_mutual, kernel_frame, surface_frame,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=Path, default=(
        ROOT / "tests/coil_from_cad/fixtures/beak_fin_48mm.step"))
    parser.add_argument("--n-peri", type=int, default=256)
    parser.add_argument("--n-stations", type=int, default=33)
    parser.add_argument("--frequency", type=float, default=150_000.0)
    parser.add_argument("--sigma", type=float, default=5.8e7)
    parser.add_argument("--pairs", type=int, default=32)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    from radia.peec_fin_topology import (
        build_hybrid_surface_topology_from_straight_prism_step,
    )

    delta = math.sqrt(2.0 / (2.0 * math.pi * args.frequency * MU0 * args.sigma))
    graph, cad = build_hybrid_surface_topology_from_straight_prism_step(
        args.step, n_peri=args.n_peri, n_stations=args.n_stations)
    fin = dict(cad)["section_analysis"].primary

    n_lanes = graph.n_lanes
    mid = (graph.n_stations - 1) // 2
    first = mid * n_lanes
    ring = graph.branch_xyz[first:first + n_lanes, 0]
    ring_next = graph.branch_xyz[first:first + n_lanes, 1]
    axis = ring_next[0] - ring[0]
    length = float(np.linalg.norm(axis))
    axis = axis / length
    xy = ring[:, :2]

    spacing = np.empty(n_lanes)
    for k in range(n_lanes):
        spacing[k] = 0.5 * (np.linalg.norm(xy[k] - xy[(k - 1) % n_lanes])
                            + np.linalg.norm(xy[(k + 1) % n_lanes] - xy[k]))

    assumed = kernel_frame(axis)
    tip = fin.tip_mask(xy)
    beak = fin.fin_mask(xy)
    step = max(1, n_lanes // args.pairs)

    rows = []
    for k in range(0, n_lanes, step):
        j = (k + 1) % n_lanes
        frame_i = surface_frame(xy[(k - 1) % n_lanes], xy[k],
                                xy[(k + 1) % n_lanes])
        frame_j = surface_frame(xy[(j - 1) % n_lanes], xy[j],
                                xy[(j + 1) % n_lanes])
        centre_i = 0.5 * (ring[k] + ring_next[k])
        centre_j = 0.5 * (ring[j] + ring_next[j])
        delta_c = centre_j - centre_i
        gap = float(np.linalg.norm(delta_c - np.dot(delta_c, axis) * axis))
        m_filament = filament_mutual(length, gap)
        m_true = bar_mutual(centre_i, centre_j, axis, length, spacing[k],
                            delta, spacing[j], delta, frame_i, frame_j)
        m_assumed = bar_mutual(centre_i, centre_j, axis, length, spacing[k],
                               delta, spacing[j], delta, assumed, assumed)
        rows.append({
            "lane": int(k),
            "region": "tip" if tip[k] else "beak" if beak[k] else "body",
            "centre_gap_m": gap,
            "band_width_m": float(spacing[k]),
            "sheet_depth_m": delta,
            "gap_over_depth": gap / delta,
            "M_filament_H": m_filament,
            "M_cross_section_averaged_H": m_true,
            "filament_relative_error": float(m_filament / m_true - 1.0),
            "orientation_sensitivity": float(m_assumed / m_true - 1.0),
        })

    err = np.array([r["filament_relative_error"] for r in rows])
    orient = np.array([r["orientation_sensitivity"] for r in rows])
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    report = {
        "schema": "radia.fin_partial_element_mutual.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_head": head,
        "host": platform.node(),
        "question": (
            "the fin branches take the filamentary Neumann mutual because "
            "they are not sub-filaments of a shared parent.  How far is that "
            "from the cross-section-averaged value on bands that touch?"),
        "fixture": str(args.step),
        "n_peri": args.n_peri,
        "n_stations": args.n_stations,
        "band_length_m": length,
        "sheet_depth_m": delta,
        "kernel_path": {
            "used": "MutualInductance (filamentary Neumann)",
            "not_used": "MutualInductanceRectBar (cross-section averaged)",
            "reason": ("the averaged kernel requires si.parent_segment == "
                       "sj.parent_segment >= 0, which holds only among the "
                       "sub-filaments of one subdivided segment; "
                       "add_connected_segment creates unsubdivided segments"),
            "verified": ("the built L matrix reproduces Grover's equal, "
                         "aligned, parallel filament formula to machine "
                         "precision at four separations"),
        },
        "filament_relative_error": {
            "min": float(err.min()), "max": float(err.max()),
            "mean": float(err.mean()), "max_abs": float(np.max(np.abs(err))),
        },
        "orientation_sensitivity_of_the_averaged_value": {
            "min": float(orient.min()), "max": float(orient.max()),
            "max_abs": float(np.max(np.abs(orient))),
            "note": ("only relevant if the averaged path is adopted; the "
                     "segment has no way to express its cross-section "
                     "orientation today"),
        },
        "pairs": rows,
        "not_claimed": (
            "this is the partial element in isolation.  It does not say what "
            "the assembled terminal impedance or the delivery metrics do, and "
            "it is not a replacement element."),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "phase": "complete",
        "filament_max_abs_relative_error":
            report["filament_relative_error"]["max_abs"],
        "filament_mean_relative_error":
            report["filament_relative_error"]["mean"],
        "orientation_max_abs":
            report["orientation_sensitivity_of_the_averaged_value"]["max_abs"],
        "output": str(args.output)}))


if __name__ == "__main__":
    main()
