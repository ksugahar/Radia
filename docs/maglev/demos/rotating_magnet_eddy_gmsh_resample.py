"""Resample the rotating-magnet eddy-current record with native Gmsh 5.

Input: the HDF5 written by ``rotating_magnet_eddy_animation_export.py`` (read
only).  Its samples are a 25 x 25 grid on the plate mid-plane, so the record is
a declared 2-D sampling surface, not a volume: Gmsh interpolates the sampled J
with Quad4 cells inside each solver step.  Query points are the cell centres
(non-node points, ids ``i * (n2 - 1) + k``).  The output is a separate
``radia.field-animation.v1`` HDF5 with the original times, solver steps and
magnet poses; ``blender_field_animation`` can render it unchanged.

The script checks the Gmsh result against an independent NumPy evaluation of
the same Quad4 interpolant (the mean of the four cell corners) and writes a
JSON report next to the output.

Run (qualified Gmsh 5 runtime, RADIA_GMSH_PROVENANCE set):
    python rotating_magnet_eddy_gmsh_resample.py --source rotating_magnet_eddy.h5 \
        --out rotating_magnet_eddy_gmsh.h5
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import h5py
import numpy as np

from radia_mcp.gmsh import sampling_surface as ss


def cell_centres(points, ids):
    normal, axes, (u, w, c0), (n1, n2) = ss.grid_association(points, ids)
    uc, wc = 0.5 * (u[1:] + u[:-1]), 0.5 * (w[1:] + w[:-1])
    q = np.empty(((n1 - 1) * (n2 - 1), 3))
    q[:, normal] = c0
    q[:, axes[0]] = np.repeat(uc, n2 - 1)
    q[:, axes[1]] = np.tile(wc, n1 - 1)
    return q, np.arange(len(q), dtype=np.int64), (n1, n2)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--steps", default="all", help="comma list of step indices or 'all'")
    a = ap.parse_args()
    report_path = Path(a.out).with_suffix(".json")
    if Path(a.out).suffix.lower() == ".json":
        raise SystemExit(f"--out {a.out} ends in .json; the report would take that name")
    if report_path.exists():
        raise SystemExit(f"report {report_path} already exists; choose a new --out")
    with h5py.File(a.source, "r") as f:
        pts, ids, J, n = f["samples/point_m"][()], f["samples/id"][()], f["field/J_A_per_m2"][()], len(f["time_s"])
    steps = list(range(n)) if a.steps == "all" else [int(s) for s in a.steps.split(",")]
    q, qid, (n1, n2) = cell_centres(pts, ids)
    t0 = time.perf_counter()
    summary = ss.resample(a.source, a.out, q, qid, steps=steps)
    seconds = time.perf_counter() - t0
    grid = J[:, np.argsort(ids)].reshape(n, n1, n2, 3)[steps]
    oracle = 0.25 * (grid[:, :-1, :-1] + grid[:, 1:, :-1] + grid[:, 1:, 1:] + grid[:, :-1, 1:])
    oracle = oracle.reshape(len(steps), -1, 3)
    with h5py.File(a.out, "r") as g:
        V = g["field/J_A_per_m2"][()]
    scale = float(np.abs(oracle).max())
    err = float(np.abs(V - oracle).max() / scale)
    report = {**summary, "seconds_resample": seconds, "shape": list(V.shape), "unit": "A/m^2",
              "oracle": "NumPy Quad4 cell-centre value (mean of the four corners)",
              "max_abs_J_A_per_m2": scale, "max_rel_error_vs_oracle": err}
    with open(report_path, "x", encoding="utf-8") as stream:  # never replace
        stream.write(json.dumps(report, indent=2))
    print(json.dumps({k: report[k] for k in ("out_sha256", "shape", "max_rel_error_vs_oracle",
                                             "seconds_resample")}))
    if not err < 1e-12:
        raise SystemExit(f"Gmsh result differs from the Quad4 oracle: {err:.3e}")


if __name__ == "__main__":
    main()
