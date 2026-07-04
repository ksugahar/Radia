"""Step-2 crossover benchmark: H-matrix-accelerated rad.Fld (_FieldEvalHMatrix) vs direct rad.Fld,
plus _FieldEvalHMatrix thread-scaling.

RUN ON mdx (idle) per the Benchmark Policy -- timing MUST NOT be measured on LAB (codex-contaminated).
Isolated env: PYTHONPATH -> C:\\temp\\radiabench\\src (HEAD build with _FieldEvalHMatrix); system deps.
MKL/OMP forced to 1 thread so TaskManager (SetNumThreads) alone owns the parallelism we sweep.

Writes results_fieldeval_hmatrix.json next to the script (Benchmark Policy schema).
"""
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import sys, json, time, platform
from datetime import datetime
import numpy as np
import psutil
import radia as rad
import radia._radia_pybind as _rp
from ngsolve import TaskManager, SetNumThreads


def peak_mem_mb():
    m = psutil.Process(os.getpid()).memory_info()
    return getattr(m, "peak_wset", m.rss) / (1024.0 * 1024.0)


def build_case(n_side):
    """n_side^3 hexahedron cubes filling [0,L]^3 (gap between cubes); obs = 1 pt just outside each +x
    face -> OUTSIDE every element (per-element superposition == rad.Fld only outside)."""
    rad.UtiDelAll()
    L = 1.0
    d = L / n_side
    h = 0.40 * d                      # half-extent -> 0.2*d gap between neighbour faces
    M = [0.0, 0.0, 954930.0]
    objs, obs = [], []
    for i in range(n_side):
        cx = (i + 0.5) * d
        for j in range(n_side):
            cy = (j + 0.5) * d
            for k in range(n_side):
                cz = (k + 0.5) * d
                v = [[cx - h, cy - h, cz - h], [cx + h, cy - h, cz - h],
                     [cx + h, cy + h, cz - h], [cx - h, cy + h, cz - h],
                     [cx - h, cy - h, cz + h], [cx + h, cy - h, cz + h],
                     [cx + h, cy + h, cz + h], [cx - h, cy + h, cz + h]]
                objs.append(rad.ObjHexahedron(v, M))
                obs.append([cx + h + 0.05 * d, cy, cz])     # outside this cube AND its +x neighbour
    cont = rad.ObjCnt(objs)
    return cont, np.asarray(obs, float)


def run_case(n_side, eps, nthreads, do_direct=True):
    cont, obs = build_case(n_side)
    N = obs.shape[0]
    obs_list = obs.tolist()
    obs_flat = obs.reshape(-1).tolist()
    SetNumThreads(nthreads)
    t_direct = None
    Bdir = None
    with TaskManager():
        if do_direct:
            t0 = time.perf_counter()
            Bdir = np.asarray(rad.Fld(cont, "b", obs_list), float).reshape(-1, 3)
            t_direct = time.perf_counter() - t0
        t0 = time.perf_counter()
        G = _rp._FieldEvalHMatrix(cont, obs_flat, eps=eps, field_type="b")
        t_build = time.perf_counter() - t0
        x = [0.0] * (3 * G.n_obs()) + list(G.src_magnetization())
        t0 = time.perf_counter()
        y = G.matvec(x)
        t_matvec = time.perf_counter() - t0
    Bhm = np.asarray(y[:3 * N], float).reshape(-1, 3)
    err = None
    if Bdir is not None:
        err = float(np.linalg.norm(Bhm - Bdir) / (np.linalg.norm(Bdir) + 1e-300))
    pk = peak_mem_mb()
    rad.UtiDelAll()
    t_total = t_build + t_matvec
    rec = dict(n_side=n_side, N_src=N, N_obs=N, embed_dof=3 * G.n_src() + 3 * G.n_obs(),
               nthreads=nthreads, eps=eps, t_direct=t_direct, t_hm_build=t_build,
               t_hm_matvec=t_matvec, t_hm_total=t_total, err=err, peak_memory_mb=pk,
               speedup_total=(t_direct / t_total) if (t_direct and t_total > 0) else None,
               converged=(err is None or err < 10 * eps + 1e-6))
    return rec


def main():
    cores = os.cpu_count()
    eps = 1e-6
    # NEAR / co-located geometry (obs interleaved with the sources): the WORST case for the symmetric embed
    # -- clusters mix obs+src, so admissible blocks have the checkerboard zero pattern (obs-obs=0, src-src=0)
    # that ACA cannot compress -> ~4-9% error floor AND ~fully-dense build (slower than direct).  Small
    # sweep is enough to show the loss; larger n only deepens it (build is O(N^2) here).  See bench_far for
    # the far-separated regime where the H-matrix WINS and is eps-accurate.
    sweep = [4, 6, 8]
    direct_cap_N = 3500
    crossover = []
    for n in sweep:
        N = n ** 3
        rec = run_case(n, eps, cores, do_direct=(N <= direct_cap_N))
        crossover.append(rec)
        sd = ("%.2f" % rec["speedup_total"]) if rec["speedup_total"] else "  -  "
        td = ("%.3f" % rec["t_direct"]) if rec["t_direct"] is not None else "skip "
        er = ("%.2e" % rec["err"]) if rec["err"] is not None else "  -  "
        print("[near] n=%2d N=%5d  t_direct=%6ss  t_hm(build+mv)=%.3f+%.3f=%.3f  speedup=%s  err=%s"
              % (n, N, td, rec["t_hm_build"], rec["t_hm_matvec"], rec["t_hm_total"], sd, er), flush=True)

    out = dict(
        timestamp=datetime.now().isoformat(),
        hostname=platform.node(),
        benchmark="fieldeval_hmatrix_near_colocated",
        problem=dict(field_type="b", eps=eps, cores=cores, geometry="n^3 hexahedron cubes, obs just "
                     "outside each +x face (obs interleaved with src -> checkerboard embed)",
                     radia_version=rad.__version__, python_version=platform.python_version(),
                     platform=platform.platform()),
        crossover=crossover,
    )
    here = os.path.dirname(os.path.abspath(__file__))
    fn = os.path.join(here, "results_fieldeval_hmatrix.json")
    with open(fn, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("Results saved to", fn, flush=True)


if __name__ == "__main__":
    main()
