"""hdiv_system_hlu_preconditioner.py -- H-LU(A) as a factor-once preconditioner for
the HDiv-VIM soft-iron material system A v = b,  A = M_mass + chi*N  (form-1 uniform
chi).  Compares GMRES iteration counts vs problem size for two preconditioners:

  - M_mass^-1            : the current production preconditioner.  Cheap to apply,
                           but iteration count GROWS with N (M_mass does not capture
                           the non-local demag operator chi*N).
  - H-LU(A)^-1           : the materialize-free system-A H-LU (_HDivVimTetSolver,
                           factored once).  H-LU(A)^-1 A ~ I, so iterations collapse
                           to O(1) and stay BOUNDED in N.

This reproduces the Phase-2 result recorded in the commit log / memory ("iters ~40
-> 1-14").  It is the evidence behind that claim, kept in committed code so it does
not live only in scratch (Data Persistence Policy).

HONEST SCOPE (do not over-read this demo): H-LU(A)^-1 is a *capability*, NOT the
default demag solver.  At the realistic VIM demag sizes (<=~30k DOF) the M_mass^-1
GMRES path is already cheap (~40 flat iterations, sub-second), so the H-LU factor
cost does not pay for itself there; large-N open-boundary demag is the FEM-Kelvin-MG
regime.  What this demo proves is narrow and real: the system-A H-LU is a CORRECT,
materialize-free, factor-once preconditioner that makes the iteration count
N-independent -- useful when many RHS / nonlinear steps amortise one factor.

Run: python hdiv_system_hlu_preconditioner.py  (writes hdiv_system_hlu_preconditioner.json)
"""
import json
import os
import platform
import sys
import time
from datetime import datetime

import inspect
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia"))
import ngsolve as ng
from netgen.csg import CSGeometry, OrthoBrick, Pnt
import radia._radia_pybind as _rp
from radia.vim import _core

HERE = os.path.dirname(os.path.abspath(__file__))
CHI = 999.0  # mu_r = 1000
GMRES_RTOL = 1e-8
_TOL_KW = "rtol" if "rtol" in inspect.signature(spla.gmres).parameters else "tol"


def _face_charge_map(d):
    """per-face (<=2 support) charge map + face centroids from build_demag output."""
    Bc = d["B_csr"].tocsc()
    cent = np.asarray(d["cent"], float)
    ndof = int(d["ndof"])
    fc = np.full(ndof * 2, -1, np.int32)
    fco = np.zeros(ndof * 2)
    fcent = np.zeros((ndof, 3))
    for i in range(ndof):
        s, e = Bc.indptr[i], Bc.indptr[i + 1]
        ids = Bc.indices[s:e][:2]
        vals = Bc.data[s:e][:2]
        for p, (a, v) in enumerate(zip(ids, vals)):
            fc[i * 2 + p] = a
            fco[i * 2 + p] = v
        fcent[i] = cent[ids].mean(0) if len(ids) else 0.0
    return fc, fco, fcent


def _gmres_iters(A, b, Minv):
    it = {"n": 0}
    x, _ = spla.gmres(A, b, M=Minv, restart=400, maxiter=2000,
                      callback=lambda rk: it.__setitem__("n", it["n"] + 1),
                      callback_type="pr_norm", **{_TOL_KW: GMRES_RTOL})
    return x, it["n"]


def run():
    ng.SetNumThreads(6)
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)))
    results = []
    print(f"A v = b,  A = M_mass + {CHI:.0f}*N  (mu_r=1000).  GMRES rtol={GMRES_RTOL:.0e}")
    print(f"{'maxh':>5} {'ndof':>6} | {'iters(M_mass^-1)':>16} {'iters(H-LU)':>12} "
          f"{'relres_HLU':>11} {'t_factor[s]':>11}")
    for maxh in (0.30, 0.22, 0.16, 0.12):
        with ng.TaskManager():
            mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh))
            d = _core.build_demag(mesh)
            Mm = d["M_mass"]
            Bcsr = d["B_csr"]
            cell_verts = np.asarray(d["cell_verts"], float)
            face_verts = np.asarray(d["face_verts"], float)
            ndof = int(d["ndof"])
            n_el = int(d["n_el"])

            # production N apply: N v = B^T G (B v)
            Hg = _rp._ChargeGramHMatrix(cell_verts=cell_verts.tolist(),
                                        face_verts=face_verts.tolist(),
                                        n_el=n_el, eps=1e-5, leaf=32, eta=2.0, near_factor=2.0)

            def A_apply(v):
                v = np.asarray(v, float)
                return np.asarray(Mm @ v).ravel() + CHI * (
                    Bcsr.T @ np.asarray(Hg.matvec((Bcsr @ v).tolist()), float))

            A = spla.LinearOperator((ndof, ndof), matvec=A_apply)
            rng = np.random.default_rng(0)
            b = rng.standard_normal(ndof)

            # preconditioner 1: M_mass^-1
            Mfac = spla.splu(sp.csc_matrix(Mm))
            P1 = spla.LinearOperator((ndof, ndof), matvec=lambda r: Mfac.solve(np.asarray(r, float)))
            _, it_mmass = _gmres_iters(A, b, P1)

            # preconditioner 2: H-LU(A)^-1 (factor once)
            fc, fco, fcent = _face_charge_map(d)
            Mc = sp.coo_matrix(Mm)
            t0 = time.time()
            S = _rp._HDivVimTetSolver(
                face_centroids=fcent.ravel().tolist(), chi=CHI,
                face_charge=fc.tolist(), face_coef=fco.tolist(),
                mI=Mc.row.tolist(), mJ=Mc.col.tolist(), mV=Mc.data.tolist(),
                cell_verts=cell_verts.tolist(), face_verts=face_verts.tolist(), n_el=n_el,
                eps=1e-5, leaf=32, eta=2.0, trunc_tol=1e-6, gram_near_factor=2.0)
            t_factor = time.time() - t0
            P2 = spla.LinearOperator(
                (ndof, ndof),
                matvec=lambda r: np.asarray(S.solve(np.asarray(r, float).tolist()), float))
            x2, it_hlu = _gmres_iters(A, b, P2)
            relres = float(np.linalg.norm(A_apply(x2) - b) / np.linalg.norm(b))

        print(f"{maxh:>5} {ndof:>6} | {it_mmass:>16} {it_hlu:>12} {relres:>11.2e} {t_factor:>11.2f}",
              flush=True)
        results.append(dict(maxh=maxh, ndof=ndof, n_el=n_el,
                            iters_mmass_inv=it_mmass, iters_hlu_inv=it_hlu,
                            relres_hlu=relres, t_factor_s=t_factor))
    return results


if __name__ == "__main__":
    results = run()
    out = dict(
        timestamp=datetime.now().isoformat(),
        hostname=platform.node(),
        benchmark="hdiv_system_hlu_preconditioner",
        problem=dict(geometry="unit cube (tet)", chi=CHI, mu_r=CHI + 1.0, gmres_rtol=GMRES_RTOL),
        note=("H-LU(A)^-1 is a factor-once preconditioner capability, NOT the default "
              "demag solver; iters bounded in N vs the growth of M_mass^-1."),
        results=results,
    )
    path = os.path.join(HERE, "hdiv_system_hlu_preconditioner.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"Results saved to {path}")
