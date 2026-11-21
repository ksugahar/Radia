"""Blender replay of the Gmsh-resampled rotating-magnet record.

Writes the .blend with ``blender_field_animation`` (no video), reopens it in a
separate Blender process and compares, for several nonadjacent steps, the
physical time, sample ids, magnet pose, rendered arrow centres (with the
display offset accounted for) and the signed J stored in the scene against the
HDF5.  It then compares the reopened signed J with an independent NumPy Quad4
evaluation from the *source* record, so the chain source -> Gmsh -> HDF5 ->
.blend is checked end to end.

Needs RADIA_BLENDER_PYTHON (Python 3.13 with bpy) and, outside an editable
checkout, RADIA_BLENDER_WORKER.  Rendering uses at most --threads CPU threads.

Run:
    python rotating_magnet_eddy_gmsh_replay_blender.py --source rotating_magnet_eddy.h5 \
        --resampled rotating_magnet_eddy_gmsh.h5 --out-dir blender_replay
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import h5py
import numpy as np

from radia_mcp.blender.field_animation import field_animation
from radia_mcp.gmsh import sampling_surface as ss


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", required=True)
    ap.add_argument("--resampled", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--steps", default="5,18,36", help="nonadjacent step indices to check")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    steps = [int(s) for s in a.steps.split(",")]
    t0 = time.perf_counter()
    rec = field_animation(a.resampled, a.out_dir,
                          {"video": False, "preview_steps": steps, "samples": 8,
                           "resolution": [640, 360], "threads": a.threads})
    seconds = time.perf_counter() - t0
    with h5py.File(a.source, "r") as f:
        pts, ids, J = f["samples/point_m"][()], f["samples/id"][()], f["field/J_A_per_m2"][()]
        t_src = f["time_s"][()]
    _, _, _, (n1, n2) = ss.grid_association(pts, ids)
    grid = J[:, np.argsort(ids)].reshape(len(t_src), n1, n2, 3)
    oracle = 0.25 * (grid[:, :-1, :-1] + grid[:, 1:, :-1] + grid[:, 1:, 1:] + grid[:, :-1, 1:])
    with h5py.File(a.resampled, "r") as g:
        src_index = g["source_step_index"][()]
        V = g["field/J_A_per_m2"][()]
        t_out = g["time_s"][()]
    per_step = {}
    for k in steps:
        o = oracle[src_index[k]].reshape(-1, 3)
        per_step[str(k)] = {
            "time_s": float(t_out[k]), "source_time_s": float(t_src[src_index[k]]),
            "hdf5_vs_source_oracle_max_rel": float(np.abs(V[k] - o).max() / np.abs(o).max()),
            "blend_readback": rec["readback"]["steps"].get(str(k)),
        }
    report = {"status": rec["status"], "failures": rec["readback"]["failures"],
              "blend": rec["artifacts"]["blend"], "sha256": rec["sha256"],
              "blender": rec["blender"], "dataset_sha256": rec["dataset"]["sha256"],
              "shape": list(V.shape), "unit": "A/m^2", "steps": per_step,
              "tolerances": rec["readback"]["tolerances"], "seconds": seconds}
    out = Path(a.out_dir) / "blender_replay.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "failures": report["failures"], "report": str(out)}))
    if rec["status"] != "ok" or any(s["hdf5_vs_source_oracle_max_rel"] > 1e-12 for s in per_step.values()):
        raise SystemExit("Blender replay mismatch")


if __name__ == "__main__":
    main()
