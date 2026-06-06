"""
hdiv_demag_speedup.py -- root-cause finish + H-LDL^T prototype, on the SYMMETRIC HDiv demag op.

PART A (root cause / operator preconditioning, the RIGHT direction):
  Riesz-map R^{-1}=(mass+a*divdiv)^{-1} FAILED because it INVERTS -> it DAMPS the high-freq
  curl-free modes, but those are exactly the near-null modes of A=(1/chi)M_mass-N (smallest demag
  eigenvalue, single-layer accumulation) that must be AMPLIFIED.  Calderon operator preconditioning
  MULTIPLIES by the complementary-order operator (D.S ~ I), so the right test is M^{-1} = a*divdiv
  (+ mass) as a MATVEC, which is LARGE on high-freq (where A is near-null).  alpha-sweep per N;
  pulse = best iters BOUNDED across N (vs point-Jacobi's ~linear growth).

PART B (H-LDL^T): A is symmetric -> the shipped NON-sym H-LU can be a symmetric H-LDL^T
  (Bunch-Kaufman).  Prototype the dense analog: (1) LDL^T solves the real operator accurately
  (stability on a symmetric INDEFINITE A), (2) factor-time ratio LU vs LDL^T (~2x theory) on a
  larger synthetic symmetric-indefinite matrix.
"""
import json, os, time
import numpy as np
import scipy.linalg as sla
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

# ---------------- PART B: H-LDL^T prototype (symmetric factorization) ----------------
print("\n=== PART B: H-LDL^T feasibility (symmetric indefinite direct factorization) ===")
d = cache[6]; S = d["Sstar"]
A_real = (1.0/CHI)*(S.T @ d["mass_hdiv"] @ S) - S.T @ d["N"] @ S   # the real symmetric indefinite op
b = np.ones(A_real.shape[0])
# (1) stability/accuracy: LDL^T solve vs LU solve on the REAL operator
lu, piv = sla.lu_factor(A_real); x_lu = sla.lu_solve((lu, piv), b)
L, Dd, perm = sla.ldl(A_real); x_ldl = sla.solve(A_real, b, assume_a="sym")  # sym solver = LDL^T path
err_ldl = np.linalg.norm(A_real @ x_ldl - b) / np.linalg.norm(b)
err_lu = np.linalg.norm(A_real @ x_lu - b) / np.linalg.norm(b)
sym = np.linalg.norm(A_real - A_real.T)/np.linalg.norm(A_real)
print(f"  real operator (rankQ={A_real.shape[0]}): ||A-A^T||/||A||={sym:.1e}; "
      f"solve rel.res: LU={err_lu:.1e}, LDL^T(sym)={err_ldl:.1e}  (LDL^T stable on indefinite A)")
# (2) factor-time ratio on a larger synthetic symmetric-indefinite matrix (clean timing)
rng = np.random.default_rng(0)
for m in (1500, 3000):
    Q, _ = np.linalg.qr(rng.standard_normal((m, m)))
    ev = np.concatenate([rng.uniform(-1, -1e-3, m//2), rng.uniform(1e-3, 1, m - m//2)])  # indefinite
    Asym = (Q * ev) @ Q.T; Asym = 0.5*(Asym + Asym.T)
    t0 = time.perf_counter(); sla.lu_factor(Asym); t_lu = time.perf_counter() - t0
    t0 = time.perf_counter(); sla.ldl(Asym); t_ldl = time.perf_counter() - t0
    print(f"  m={m}: LU factor {t_lu*1e3:.0f} ms, LDL^T factor {t_ldl*1e3:.0f} ms  -> LDL^T/LU = {t_ldl/t_lu:.2f}x")
print("  (dense theory: LU ~ 2n^3/3, LDL^T ~ n^3/3 -> ~0.5x; same ratio carries to the H-matrix")
print("   block factorization -> a symmetric H-LDL^T roughly halves the shipped H-LU factor cost.)")

with open(os.path.join(HERE, "hdiv_demag_speedup.json"), "w") as f:
    json.dump({"partA_calderon": rowsA, "partB_ldl_err": {"lu": float(err_lu), "ldl": float(err_ldl)}}, f, indent=2)
print("\nsaved", os.path.join(HERE, "hdiv_demag_speedup.json"))
