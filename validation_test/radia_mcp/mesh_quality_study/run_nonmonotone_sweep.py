"""Map netgen's min-quality NON-MONOTONICITY under uniform refinement.

Promoted from C:/temp/mesh_quality_study (2026-08-06) with its committed
results JSON (Data Persistence Policy). Re-run with
`python run_nonmonotone_sweep.py` (requires Cubit + netgen + gmsh +
build123d; scratch meshes land in artifacts/, gitignored).

The equal-budget study (run_study.py) found one striking data point:
calibrating c_core from maxh 8.0 to 6.97 made the mesh FINER but dropped
min minSICN from 0.47 to 0.12. This script asks whether that was a fluke
or a systematic property, by sweeping maxh continuously (24 points per
geometry) and counting refinement steps where the finer mesh has a WORSE
worst element. Cubit tetmesh at three sizes provides the reference band.

Quality-class run (correctness, not timing) -- LAB execution allowed.
"""
import json
import os
import platform
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

# a refinement step counts as non-monotone when min quality drops by more
# than this (0.02 absorbs the referee's own reproducibility noise)
DROP_TOL = 0.02
N_POINTS = 24


def summarize(q):
    bts = q.get("by_type") or []
    n = sum(bt["n_elements"] for bt in bts)
    return {"n": n,
            "min": round(min(bt["min_quality"] for bt in bts), 4),
            "mean": round(sum(bt["mean_quality"] * bt["n_elements"]
                              for bt in bts) / n, 4)}


def netgen_point(step, maxh, tag):
    msh = Path(OUT) / f"sw_{tag}.msh"
    _netgen_mesh_to_msh(Path(step), maxh, msh)
    s = summarize(mesh_quality(msh))
    s["maxh"] = round(maxh, 4)
    return s


def cubit_point(step, size, tag):
    msh = os.path.join(OUT, f"sw_{tag}_cubit.msh")
    r = _run_batch(step, ["volume all scheme tetmesh",
                          f"volume all size {size}", "mesh volume all",
                          "block 1 add volume all", 'block 1 name "mesh"',
                          f'export gmsh "{msh.replace(os.sep, "/")}" overwrite'],
                   timeout_s=600)
    assert r["status"] == "ok", r
    s = summarize(mesh_quality(msh))
    s["size"] = size
    return s


def frange(a, b, steps):
    return [a + (b - a) * i / (steps - 1) for i in range(steps)]


def main():
    t0 = time.time()
    cc = os.path.join(OUT, "c_core.step")
    sp = os.path.join(OUT, "sphere.step")
    export_step(c_core(width=80, height=60, depth=25, leg=15, gap=8), cc)
    export_step(Sphere(1.0), sp)

    results = {"timestamp": datetime.now().isoformat(),
               "hostname": platform.node(),
               "referee": "gmsh minSICN", "order": 1,
               "drop_tolerance": DROP_TOL, "n_points": N_POINTS,
               "sweeps": {}}

    for name, step, lo, hi, sizes in (
            ("c_core", cc, 2.5, 10.0, [8.0, 5.0, 3.0]),
            ("sphere", sp, 0.2, 0.7, [0.5, 0.35, 0.25])):
        pts = []
        for i, maxh in enumerate(frange(hi, lo, N_POINTS)):
            s = netgen_point(step, maxh, f"{name}_{i:02d}")
            pts.append(s)
            print(f"{name} maxh={s['maxh']:7.4f}  n={s['n']:6d} "
                  f"min={s['min']:.3f} mean={s['mean']:.3f}", flush=True)
        cub = [cubit_point(step, sz, f"{name}_{sz}") for sz in sizes]
        for s in cub:
            print(f"{name} CUBIT size={s['size']:5.2f} n={s['n']:6d} "
                  f"min={s['min']:.3f} mean={s['mean']:.3f}", flush=True)
        mins = [p["min"] for p in pts]
        # refinement steps where the FINER mesh has a WORSE worst element
        drops = [(pts[i]["maxh"], round(mins[i] - mins[i + 1], 3))
                 for i in range(len(pts) - 1)
                 if mins[i + 1] < mins[i] - DROP_TOL]
        results["sweeps"][name] = {
            "netgen": pts, "cubit_ref": cub,
            "netgen_min_range": [min(mins), max(mins)],
            "cubit_min_range": [min(c["min"] for c in cub),
                                max(c["min"] for c in cub)],
            "non_monotone_drops": drops,
        }
        print(f"{name}: netgen min range {min(mins):.3f}-{max(mins):.3f}, "
              f"cubit band {min(c['min'] for c in cub):.3f}-"
              f"{max(c['min'] for c in cub):.3f}, "
              f"{len(drops)} non-monotone drops\n", flush=True)

    out = os.path.join(_HERE, "results_nonmonotone_sweep.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    print(f"saved {out} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
