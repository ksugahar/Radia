"""hdiv_matvec_scaling.py -- the SCALABILITY evidence for the HDiv-VIM production demag solve.

The production path (radia.vim.hdiv_demag_solve) solves the form-1 soft-iron system
    A m = M_chi h_ext ,   A = M_mass + chi*N   (uniform chi; bit-for-bit the production system),
with N = B^T G B applied as the HACApK analytic charge-Gram H-matrix MATVEC and GMRES
preconditioned by the sparse M_mass^-1.  This is NOT the H-LU factorization path
(hdiv_system_hlu_preconditioner.py) -- it is the cheaper, default, matrix-free solve.

The two scalability claims this benchmark MEASURES across a DOF sweep (unit-cube tet,
n_face ~ 0.4k -> ~50k):

  (1) STORAGE is O(N log N):  the charge-Gram H-matrix `memory_mb` vs the dense O(N^2)
      `dense_memory_mb` -- the compression ratio grows with N (sub-quadratic storage).
      THIS IS LOAD-INDEPENDENT (the H-matrix size is deterministic), so it is the clean
      sub-cubic evidence even on a busy machine.

  (2) ITERATION COUNT is FLAT in N:  the loop-de-Rham modes sit in ker(N) and are carried
      by M_mass (form-1 / no loop-star), so the M_mass^-1-GMRES iteration count does NOT
      grow with N.  Also load-independent (iteration count is deterministic).

Together: one solve = O(1) iters x O(N log N) matvec = O(N log N).  HACApK (the charge-Gram
H-matvec) is the scalable engine; the loop-null form-1 keeps the iteration count bounded.

HONEST on TIMING: the wall-clock columns (t_build_s, t_solve_s) are recorded for reference
but are NOT the headline -- they are unreliable on a loaded machine.  The headline metrics
(memory_mb / compression / iters / relres) are deterministic and valid regardless of load.
Re-run the timing on a QUIET machine for the wall-clock power-law fit.

Run: python hdiv_matvec_scaling.py  (writes hdiv_matvec_scaling.json incrementally)
"""
import inspect
import json
import os
import platform
import sys
import time
from datetime import datetime

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

try:
    import psutil
    def _peak_mb():
        m = psutil.Process(os.getpid()).memory_info()
        return (m.peak_wset if hasattr(m, "peak_wset") else m.rss) / (1024 * 1024)
except Exception:
    def _peak_mb():
        return 0.0

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia"))
import ngsolve as ng
from netgen.csg import CSGeometry, OrthoBrick, Pnt
import radia._radia_pybind as _rp
from radia.vim import _core

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "hdiv_matvec_scaling.json")
CHI = 999.0                 # mu_r = 1000
GMRES_RTOL = 1e-8
GRAM_EPS = 1e-5             # ACA tolerance for the charge-Gram H-matrix
MAXH = [0.30, 0.22, 0.16, 0.12, 0.09, 0.07, 0.055]   # -> n_face ~ 0.4k .. ~50k
_TOL_KW = "rtol" if "rtol" in inspect.signature(spla.gmres).parameters else "tol"


def _one(maxh, geo):
    with ng.TaskManager():
        mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh))
        d = _core.build_demag(mesh, skip_dense_gram=True)
        Mm = d["M_mass"]; Bcsr = d["B_csr"]
        cell_verts = np.asarray(d["cell_verts"], float)
        face_verts = np.asarray(d["face_verts"], float)
        ndof = int(d["ndof"]); n_el = int(d["n_el"]); n_charge = int(d["n_charge"])

        t0 = time.time()
        Hg = _rp._ChargeGramHMatrix(cell_verts=cell_verts.tolist(), face_verts=face_verts.tolist(),
                                    n_el=n_el, eps=GRAM_EPS, leaf=32, eta=2.0, near_factor=2.0)
        t_build = time.time() - t0
        st = Hg.stats()

        def A_apply(v):
            v = np.asarray(v, float)
            return np.asarray(Mm @ v).ravel() + CHI * (
                Bcsr.T @ np.asarray(Hg.matvec((Bcsr @ v).tolist()), float))

        A = spla.LinearOperator((ndof, ndof), matvec=A_apply)
        Mfac = spla.splu(sp.csc_matrix(Mm))
        Minv = spla.LinearOperator((ndof, ndof), matvec=lambda r: Mfac.solve(np.asarray(r, float)))

        rng = np.random.default_rng(0)
        b = rng.standard_normal(ndof)
        it = {"n": 0}
        t0 = time.time()
        x, info = spla.gmres(A, b, M=Minv, restart=400, maxiter=2000,
                             callback=lambda rk: it.__setitem__("n", it["n"] + 1),
                             callback_type="pr_norm", **{_TOL_KW: GMRES_RTOL})
        t_solve = time.time() - t0
        relres = float(np.linalg.norm(A_apply(x) - b) / np.linalg.norm(b))

    return dict(
        maxh=maxh, ndof=ndof, n_el=n_el, n_charge=n_charge,
        # --- LOAD-INDEPENDENT headline metrics ---
        hmat_memory_mb=float(st.get("memory_mb", 0.0)),
        dense_memory_mb=float(st.get("dense_memory_mb", 0.0)),
        compression=float(st.get("compression", 0.0)),
        max_rank=int(st.get("max_rank", 0)),
        n_lowrank=int(st.get("n_lowrank", 0)),
        n_dense=int(st.get("n_dense", 0)),
        iters=int(it["n"]), gmres_info=int(info), relres=relres,
        # --- reference-only (load-contaminated) timing ---
        t_build_s=t_build, t_solve_s=t_solve, peak_memory_mb=_peak_mb(),
    )


def run():
    ng.SetNumThreads(4)   # be a good citizen under load
    geo = CSGeometry(); geo.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)))
    results = []
    print(f"form-1 A=M_mass+{CHI:.0f}*N, charge-Gram H-matvec, M_mass^-1 GMRES (rtol={GMRES_RTOL:.0e})")
    print(f"{'ndof':>6} {'n_chg':>6} | {'Hmat_MB':>9} {'dense_MB':>10} {'compr':>7} {'rank':>5} "
          f"| {'iters':>5} {'relres':>10} | {'t_build':>8} {'t_solve':>8}")
    for maxh in MAXH:
        try:
            r = _one(maxh, geo)
        except Exception as e:
            print(f"maxh={maxh}: FAILED {type(e).__name__}: {e}", flush=True)
            break
        results.append(r)
        print(f"{r['ndof']:>6} {r['n_charge']:>6} | {r['hmat_memory_mb']:>9.2f} "
              f"{r['dense_memory_mb']:>10.1f} {r['compression']:>7.3f} {r['max_rank']:>5} "
              f"| {r['iters']:>5} {r['relres']:>10.2e} | {r['t_build_s']:>8.1f} {r['t_solve_s']:>8.2f}",
              flush=True)
        _save(results)   # incremental: a partial sweep still yields committed data
    return results


def _save(results):
    out = dict(
        timestamp=datetime.now().isoformat(), hostname=platform.node(),
        benchmark="hdiv_matvec_scaling",
        problem=dict(geometry="unit cube (tet)", chi=CHI, mu_r=CHI + 1.0,
                     gmres_rtol=GMRES_RTOL, gram_eps=GRAM_EPS, nthreads=4),
        note=("Scalability of the PRODUCTION matrix-free demag solve: charge-Gram H-matrix storage "
              "(memory_mb vs dense_memory_mb) is O(N log N) and the loop-null M_mass^-1-GMRES iteration "
              "count is flat in N -> O(N log N) per solve.  memory_mb/compression/iters/relres are "
              "load-INDEPENDENT (deterministic); t_build_s/t_solve_s are reference-only (re-time quiet)."),
        results=results,
    )
    with open(OUT, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    run()
    print(f"Results saved to {OUT}")
