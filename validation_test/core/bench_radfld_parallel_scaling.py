"""rad.Fld batch parallel scaling -- regression guard for the allocation-free per-triangle field.

Background (2026-07-04): rad.Fld's obs-point batch (radTApplication::ComputeFieldBatch, an
ngcore::ParallelFor over points) used to SLOW DOWN with threads (0.4x at 38 cores on idle mdx) because
the per-triangle field RadFieldFromTriangleFaceGlobal allocated EIGHT std::vectors + ran a NESTED
ngcore::ParallelFor per call (thousands of calls per observation point) -> Windows CRT heap-lock
contention.  Delegating it to the allocation-free RadFieldFromTriangleFaceWithBasis (the closed form the
HACApK path already used, which scales) fixed it: single-thread ~3x faster AND ~20x thread scaling.
A RegionTaskManager self-wrap in the batch entry makes a BARE rad.Fld (no caller `with TaskManager()`)
run parallel too.

RUN ON mdx (idle) per the Benchmark Policy -- timing must not be measured on the codex-contended LAB.
Isolated env: PYTHONPATH -> a HEAD src/radia copy over the system Python; MKL/OMP=1 so TaskManager owns
the swept parallelism.  Writes results_radfld_parallel_scaling.json next to this file.

  python bench_radfld_parallel_scaling.py
"""
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import json, time, platform
from datetime import datetime
import numpy as np
import radia as rad
from ngsolve import SetNumThreads, TaskManager


def build(nsrc_side=6):
    rad.UtiDelAll()
    d = 1.0 / nsrc_side
    h = 0.40 * d
    M = [0.0, 0.0, 954930.0]
    objs = []
    for i in range(nsrc_side):
        cx = (i + 0.5) * d
        for j in range(nsrc_side):
            cy = (j + 0.5) * d
            for k in range(nsrc_side):
                cz = (k + 0.5) * d
                v = [[cx - h, cy - h, cz - h], [cx + h, cy - h, cz - h],
                     [cx + h, cy + h, cz - h], [cx - h, cy + h, cz - h],
                     [cx - h, cy - h, cz + h], [cx + h, cy - h, cz + h],
                     [cx + h, cy + h, cz + h], [cx - h, cy + h, cz + h]]
                objs.append(rad.ObjHexahedron(v, M))
    return rad.ObjCnt(objs)


def main():
    cores = os.cpu_count()
    cont = build(6)                                  # 216 hexahedron sources
    N = 8000
    pts = (np.random.RandomState(0).rand(N, 3) * 2.0 + 1.5).tolist()   # far obs
    warm = pts[:400]
    thread_list = sorted({t for t in [1, 2, 4, 8, 16, cores] if t <= cores})

    def measure(nt, wrapped):
        SetNumThreads(nt)
        if wrapped:
            with TaskManager():
                rad.Fld(cont, "b", warm)
            best = 1e9
            for _ in range(3):
                with TaskManager():
                    t = time.perf_counter(); rad.Fld(cont, "b", pts); best = min(best, time.perf_counter() - t)
        else:
            rad.Fld(cont, "b", warm)                 # BARE -- no caller TaskManager (self-wrap must parallelise)
            best = 1e9
            for _ in range(3):
                t = time.perf_counter(); rad.Fld(cont, "b", pts); best = min(best, time.perf_counter() - t)
        return best

    results = []
    base_bare = None
    for nt in thread_list:
        tb = measure(nt, wrapped=False)
        tw = measure(nt, wrapped=True)
        if base_bare is None:
            base_bare = tb
        results.append(dict(nthreads=nt, t_bare=tb, t_wrapped=tw, speedup_bare_vs_1t=base_bare / tb))
        print("threads=%2d  bare=%7.3fs (%.1fx)  wrapped=%7.3fs" % (nt, tb, base_bare / tb, tw), flush=True)

    out = dict(timestamp=datetime.now().isoformat(), hostname=platform.node(),
               benchmark="radfld_batch_parallel_scaling",
               problem=dict(cores=cores, N_src=216, N_obs=N, field="b",
                            geometry="216 ObjHexahedron sources, 8000 far obs points",
                            radia_version=rad.__version__, python_version=platform.python_version(),
                            platform=platform.platform(),
                            note="bare = no caller TaskManager (exercises the RegionTaskManager self-wrap)"),
               results=results)
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_radfld_parallel_scaling.json")
    with open(fn, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("Results saved to", fn, flush=True)
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
