"""
hdiv_demag_speedup.py -- root-cause finish: Calderon operator preconditioning of the SYMMETRIC HDiv demag op.

PART A (root cause / operator preconditioning, the RIGHT direction):
  Riesz-map R^{-1}=(mass+a*divdiv)^{-1} FAILED because it INVERTS -> it DAMPS the high-freq
  curl-free modes, but those are exactly the near-null modes of A=(1/chi)M_mass-N (smallest demag
  eigenvalue, single-layer accumulation) that must be AMPLIFIED.  Calderon operator preconditioning
  MULTIPLIES by the complementary-order operator (D.S ~ I), so the right test is M^{-1} = a*divdiv
  (+ mass) as a MATVEC, which is LARGE on high-freq (where A is near-null).  alpha-sweep per N;
  pulse = best iters BOUNDED across N (vs point-Jacobi's ~linear growth).

(A H-LDL^T direct-factorization prototype was removed 2026-06-08: the HDiv-VIM is mu_r-independent, so
 the iterative solve is the production path; the symmetric H-LDL^T H-factorization was deleted.)
"""
import json, os
import numpy as np
from scipy.sparse.linalg import gmres, LinearOperator
from hdiv_demag_quad_self import build

HERE = os.path.dirname(os.path.abspath(__file__))
MUR = 1e4; CHI = MUR - 1.0

def gmres_iters(A, Minv):
    n = A.shape[0]; it = {"n": 0}
    gmres(A, np.ones(n), M=Minv, rtol=1e-8, restart=n, maxiter=2,
          callback=lambda xk: it.__setitem__("n", it["n"]+1), callback_type="pr_norm")
    return it["n"]

# ---------------- PART A: Calderon-direction operator preconditioning ----------------
print(f"=== PART A: Calderon-direction op-precond  M^-1 = a*divdiv + mass  (MATVEC), mu_r={MUR:.0e} ===")
print(f"  (multiply by the complementary operator, NOT invert it)")
print(f"  {'n':>2} {'rankQ':>5} {'pointJac':>9} {'best Calderon (alpha*)':>24}")
cache = {}; rowsA = []
for n in (3, 4, 5, 6):
    d = build(0.12, 4, 6, n=n); cache[n] = d; S = d["Sstar"]; rq = d["rankQ"]
    Mst = S.T @ d["mass_hdiv"] @ S; Dst = S.T @ d["divdiv_hdiv"] @ S
    A = (1.0/CHI)*Mst - S.T @ d["N"] @ S
    dpt = np.diag(A).copy(); dpt[np.abs(dpt) < 1e-30] = 1.0
    it_pt = gmres_iters(A, LinearOperator((rq, rq), matvec=lambda x: x/dpt))
    best = (10**9, None)
    for a in (1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4):
        M = a*Dst + Mst
        it = gmres_iters(A, LinearOperator((rq, rq), matvec=lambda x, MM=M: MM @ x))
        if it < best[0]: best = (it, a)
    print(f"  {n:>2} {rq:>5} {it_pt:>9} {best[0]:>14} (a={best[1]:.0e})")
    rowsA.append({"n": n, "rankQ": int(rq), "pointJac": it_pt, "calderon_best": int(best[0]), "alpha": best[1]})
pj = [r["pointJac"] for r in rowsA]; cb = [r["calderon_best"] for r in rowsA]
print(f"\n  point-Jacobi   : {pj[0]} -> {pj[-1]} ({pj[-1]/pj[0]:.1f}x)")
print(f"  Calderon best  : {cb[0]} -> {cb[-1]} ({cb[-1]/cb[0]:.1f}x)   "
      f"{'<- PULSE (bounded/sublinear)!' if cb[-1]/cb[0] < pj[-1]/pj[0]*0.6 else '<- still grows like point-Jacobi'}")

with open(os.path.join(HERE, "hdiv_demag_speedup.json"), "w") as f:
    json.dump({"partA_calderon": rowsA}, f, indent=2)
print("\nsaved", os.path.join(HERE, "hdiv_demag_speedup.json"))
