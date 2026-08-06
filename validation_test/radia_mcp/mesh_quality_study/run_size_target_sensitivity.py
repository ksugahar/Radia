"""Is a mesher's MIN quality a reproducible property? (control for the sweep)

Promoted from C:/temp/mesh_quality_study (2026-08-06) with its committed
results JSON (Data Persistence Policy). Re-run with
`python run_size_target_sensitivity.py` (requires Cubit + netgen + gmsh +
build123d; scratch meshes land in artifacts/, gitignored).

run_nonmonotone_sweep.py counted refinement steps where the finer mesh had
a WORSE worst element. That count only means something if min quality is
otherwise a stable function of the size target. This script measures the
two ways it could fail to be:

  A. REPEATABILITY -- re-mesh the same STEP at the same target N times.
  B. TARGET SENSITIVITY -- mesh across a ~0.15 % window of the size
     target (a change far below any meaningful refinement step).

Both meshers are measured the same way, so neither claim is asymmetric.

Quality-class run (correctness, not timing) -- LAB execution allowed.
"""
import json
import os
import platform
import statistics
import time
from datetime import datetime
from pathlib import Path

from build123d import Sphere, export_step
from radia_mcp.build123d.archetypes import c_core
from radia_mcp.cubit.server import _netgen_mesh_to_msh, _run_batch
from radia_mcp.gmsh.msh_inspect import mesh_quality

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, "artifacts")
os.makedirs(OUT, exist_ok=True)

REPEATS = 5


def summarize(q):
    bts = q.get("by_type") or []
    n = sum(bt["n_elements"] for bt in bts)
    return {"n": n,
            "min": round(min(bt["min_quality"] for bt in bts), 4),
            "mean": round(sum(bt["mean_quality"] * bt["n_elements"]
                              for bt in bts) / n, 4)}


def netgen_once(step, maxh, tag):
    msh = Path(OUT) / f"st_{tag}.msh"
    _netgen_mesh_to_msh(Path(step), maxh, msh)
    return summarize(mesh_quality(msh))


def cubit_once(step, size, tag):
    msh = os.path.join(OUT, f"st_{tag}.msh")
    r = _run_batch(step, ["volume all scheme tetmesh",
                          f"volume all size {size}", "mesh volume all",
                          "block 1 add volume all", 'block 1 name "mesh"',
                          f'export gmsh "{msh.replace(os.sep, "/")}" overwrite'],
                   timeout_s=600)
    assert r["status"] == "ok", r
    return summarize(mesh_quality(msh))


def _stats(runs, key_name, keys):
    mins = [r["min"] for r in runs]
    ns = [r["n"] for r in runs]
    return {key_name: keys,
            "min_values": mins,
            "min_range": [min(mins), max(mins)],
            "min_spread": round(max(mins) - min(mins), 4),
            "min_stdev": round(statistics.pstdev(mins), 4),
            "n_values": ns,
            "n_spread": max(ns) - min(ns)}


def window(center, n=10):
    """~0.15 % window around the size target, denser near the centre."""
    return [round(center * (1.0 + d), 6) for d in
            (0.0, 2.5e-4, 5.0e-4, 5.5e-4, 6.0e-4, 6.5e-4,
             7.0e-4, 7.5e-4, 1.0e-3, 1.25e-3)][:n]


def main():
    t0 = time.time()
    cc = os.path.join(OUT, "c_core.step")
    sp = os.path.join(OUT, "sphere.step")
    export_step(c_core(width=80, height=60, depth=25, leg=15, gap=8), cc)
    export_step(Sphere(1.0), sp)

    results = {"timestamp": datetime.now().isoformat(),
               "hostname": platform.node(),
               "referee": "gmsh minSICN", "order": 1, "repeats": REPEATS,
               "repeatability": {}, "target_sensitivity": {}}

    # ---- A. repeatability at a FIXED target -------------------------
    for name, kind, step, size in (
            ("netgen_c_core_6.5", "netgen", cc, 6.5),
            ("netgen_sphere_0.3", "netgen", sp, 0.3),
            ("cubit_c_core_5.0", "cubit", cc, 5.0),
            ("cubit_sphere_0.35", "cubit", sp, 0.35)):
        runs = [(netgen_once(step, size, f"rp_{name}_{i}") if kind == "netgen"
                 else cubit_once(step, size, f"rp_{name}_{i}"))
                for i in range(REPEATS)]
        s = _stats(runs, "target", size)
        results["repeatability"][name] = s
        print(f"[repeat ] {name:20s} min={s['min_values']} "
              f"spread={s['min_spread']:.4f}", flush=True)

    # ---- B. sensitivity to a TINY change of the target --------------
    for name, kind, step, center in (
            ("netgen_c_core", "netgen", cc, 3.1),
            ("cubit_c_core", "cubit", cc, 4.0)):
        keys = window(center)
        runs = [(netgen_once(step, k, f"ws_{name}_{i}") if kind == "netgen"
                 else cubit_once(step, k, f"ws_{name}_{i}"))
                for i, k in enumerate(keys)]
        s = _stats(runs, "targets", keys)
        s["window_pct"] = round((max(keys) - min(keys)) / min(keys) * 100, 4)
        results["target_sensitivity"][name] = s
        print(f"[window ] {name:20s} over {s['window_pct']:.3f}% of target: "
              f"min {s['min_range'][0]:.4f}-{s['min_range'][1]:.4f} "
              f"(spread {s['min_spread']:.4f}), n {min(s['n_values'])}-"
              f"{max(s['n_values'])}", flush=True)

    out = os.path.join(_HERE, "results_size_target_sensitivity.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    print(f"saved {out} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
