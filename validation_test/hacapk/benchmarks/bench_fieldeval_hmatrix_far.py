"""Far-separated crossover: src box and obs box separated by a clean gap so the cluster tree splits them
into DISTINCT obs / src clusters (no checkerboard) -- the path-B far-field-map regime.  Confirms
(a) accuracy recovers to eps-level when obs/src are separable, and (b) where H-matrix beats direct.

Run on mdx (idle).  PYTHONPATH -> C:\\temp\\radiabench\\src.  MKL/OMP=1 -> TaskManager owns parallelism.
Writes results_fieldeval_far.json next to the script.
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


def build_far(n_side, S=0.2, gap_in_S=4.0):
    """src = n^3 cubes filling [0,S]^3; obs = n^3 grid in a box shifted +z by (1+gap_in_S)*S -> clean
    gap (obs/src separable into distinct clusters)."""
    rad.UtiDelAll()
    d = S / n_side
    h = 0.40 * d
    M = [0.0, 0.0, 954930.0]
    objs = []
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
    cont = rad.ObjCnt(objs)
    zoff = (1.0 + gap_in_S) * S
    lin = (np.arange(n_side) + 0.5) * d
    ox, oy, oz = np.meshgrid(lin, lin, lin + zoff, indexing="ij")
    obs = np.column_stack([ox.ravel(), oy.ravel(), oz.ravel()])
    return cont, obs


def run_case(n_side, eps, nthreads, do_direct=True):
    cont, obs = build_far(n_side)
    N = obs.shape[0]
    obs_list = obs.tolist()
    obs_flat = obs.reshape(-1).tolist()
    SetNumThreads(nthreads)
    t_direct, Bdir = None, None
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
    err = None if Bdir is None else float(np.linalg.norm(Bhm - Bdir) / np.linalg.norm(Bdir))
    t_total = t_build + t_matvec
    rad.UtiDelAll()
    return dict(n_side=n_side, N_src=N, N_obs=N, nthreads=nthreads, eps=eps, t_direct=t_direct,
                t_hm_build=t_build, t_hm_matvec=t_matvec, t_hm_total=t_total, err=err,
                peak_memory_mb=peak_mem_mb(),
                speedup_total=(t_direct / t_total) if (t_direct and t_total > 0) else None)


def main():
    cores = os.cpu_count()
    eps = 1e-6
    sweep = [4, 6, 8, 10, 12, 14]
    direct_cap_N = 3500
    results = []
    for n in sweep:
        N = n ** 3
        r = run_case(n, eps, cores, do_direct=(N <= direct_cap_N))
        results.append(r)
        sd = ("%.2f" % r["speedup_total"]) if r["speedup_total"] else "  -  "
        td = ("%.3f" % r["t_direct"]) if r["t_direct"] is not None else "skip "
        er = ("%.2e" % r["err"]) if r["err"] is not None else "  -  "
        print("[far] n=%2d N=%5d  t_direct=%6ss  t_hm=%.3f+%.3f=%.3f  speedup=%s  err=%s"
              % (n, N, td, r["t_hm_build"], r["t_hm_matvec"], r["t_hm_total"], sd, er), flush=True)
    # accuracy vs eps at a fixed multi-cluster far case
    acc = []
    for e in [1e-3, 1e-6, 1e-9]:
        r = run_case(8, e, cores, do_direct=True)
        acc.append(dict(eps=e, err=r["err"], t_hm_total=r["t_hm_total"]))
        print("[far-acc] n=8 eps=%.0e err=%.3e" % (e, r["err"]), flush=True)
    out = dict(timestamp=datetime.now().isoformat(), hostname=platform.node(),
               benchmark="fieldeval_hmatrix_far_separated",
               problem=dict(field_type="b", eps=eps, cores=cores,
                            geometry="src n^3 cubes in [0,S]^3, obs n^3 grid shifted +z by 5S (clean gap)",
                            radia_version=rad.__version__, python_version=platform.python_version(),
                            platform=platform.platform()),
               results=results, accuracy_vs_eps=acc)
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_fieldeval_far.json")
    with open(fn, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("Results saved to", fn, flush=True)


if __name__ == "__main__":
    main()
