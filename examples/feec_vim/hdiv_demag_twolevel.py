"""
hdiv_demag_twolevel.py -- go/no-go for "H-AMS / coarse-space preconditioning" on the SYMMETRIC
HDiv-type demag operator.

hdiv_demag_blockjacobi.py showed a NAIVE ADDITIVE coarse correction (D^-1 + Q) does NOT bound the
N-growth.  But additive is the wrong design for indefinite A.  Here we test the PROPER SYMMETRIC
TWO-LEVEL (deflation-corrected smoother), with the ORACLE coarse space W = k smallest-|lambda|
eigvecs of A (the best possible coarse space):

  Q   = W (W^T A W)^{-1} W^T                              (coarse solve)
  S   = signed-diag(A)^{-1}                               (point smoother)
  M^{-1} = Q + (I - Q A) S (I - A Q)                      (symmetric hybrid two-level)

This is the textbook two-level preconditioner.  Question: with the IDEAL coarse space in the
CORRECT form, do GMRES iters stay BOUNDED as N grows (fixed k, and k~ndof)?
  - BOUNDED with fixed/modest k -> a coarse space DOES fix the N-scaling => H-AMS-style
    (operator-preconditioning: dense H-matrix apply + a SPARSE scalar-potential coarse space)
    is worth building.  The HDiv SYMMETRY is what makes this work (it failed for non-sym MSC).
  - NOT bounded even with the oracle coarse space -> coarse-space methods won't help; H-ILU stays
    the answer.  (Decisive: if the BEST coarse space can't do it, no cheap surrogate will.)

RESULT (2026-06-07, DECISIVE -- H-AMS / coarse-space preconditioning is NOT promising here):
  The oracle two-level (deflate the 24 smallest-|lambda| eigvecs of A, correct symmetric hybrid
  form) does NOT improve conditioning: cond(M2.A)=1.92e4 == cond(Mpoint.A)=1.93e4, preconditioned
  spectrum min stays 2.36e-4, GMRES iters identical to point-Jacobi (70 @ n=3).  Deflating MORE
  modes (k~N/6) made it WORSE (485) -- the smallest-|lambda| eigvecs of an INDEFINITE A are a bad
  coarse basis.  Physical reason: the demag operator's spectrum is a SPREAD INDEFINITE CONTINUUM,
  not a few isolated near-null modes.  Its ONE genuine kernel (the loops) is already handled by
  construction (field-null).  The remaining near-field-null STAR modes GROW IN NUMBER with N, so
  no FIXED low-dim coarse space keeps up -> AMS-style coarse preconditioning cannot bound the
  N-growth.  => the right tool is a HIERARCHICAL factorization (H-ILU/H-LU) that captures the FULL
  operator, NOT a coarse-space auxiliary -- not H-AMS.  (This is NOT a sparse-vs-dense issue: AMS's
  sparse auxiliary operators exist, but the spectral STRUCTURE AMS exploits is absent.)  NOTE: the
  production HDiv-VIM path is the mu_r-independent ITERATIVE solve (MINRES/CG/Picard); a symmetric
  H-LDL^T direct-factorization variant was prototyped but REMOVED 2026-06-08 -- N-scalability, if it
  proves necessary, goes via the kept H-LU/H-ILU or a better preconditioner (the open M0 question).
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

def two_level_M(A, k):
    w, V = np.linalg.eigh(A)
    W = V[:, np.argsort(np.abs(w))[:k]]                 # oracle coarse space (smallest |lambda|)
    AW = A @ W; Ac_inv = np.linalg.inv(W.T @ AW)
    d = np.diag(A).copy(); d[np.abs(d) < 1e-30] = 1.0   # signed-diag smoother
    def apply(v):
        qv = W @ (Ac_inv @ (W.T @ v))                   # Q v
        r  = v - A @ qv                                 # (I - A Q) v
        sr = r / d                                      # S r
        return qv + (sr - W @ (Ac_inv @ (W.T @ (A @ sr))))   # Q v + (I - Q A) S (I - A Q) v
    return LinearOperator(A.shape, matvec=apply)

def point_M(A):
    d = np.diag(A).copy(); d[np.abs(d) < 1e-30] = 1.0
    return LinearOperator(A.shape, matvec=lambda x: x/d)

res = {"mu_r": MUR, "rows": []}
print(f"=== PROPER symmetric two-level (oracle coarse) vs point-Jacobi, distorted 0.12, mu_r={MUR:.0e} ===")
print(f"  {'n':>2} {'ndof':>5} {'pointJac':>9} {'2lvl k=12':>10} {'2lvl k=24':>10} {'2lvl k~ndof/6':>14}")
for n in (2, 3, 4, 5, 6):
    d = build(0.12, 4, 6, n=n); A = (1.0/CHI)*np.eye(d["ndof"]) - d["N"]; nd = d["ndof"]
    it_pt = gmres_iters(A, point_M(A))
    it_12 = gmres_iters(A, two_level_M(A, min(12, nd-1)))
    it_24 = gmres_iters(A, two_level_M(A, min(24, nd-1)))
    it_pp = gmres_iters(A, two_level_M(A, max(2, nd//6)))
    print(f"  {n:>2} {nd:>5} {it_pt:>9} {it_12:>10} {it_24:>10} {it_pp:>14}")
    res["rows"].append({"n": n, "ndof": int(nd), "pointJac": it_pt,
                        "twolvl_k12": it_12, "twolvl_k24": it_24, "twolvl_kprop": it_pp})

r = res["rows"]
def grow(key): return f"{r[0][key]} -> {r[-1][key]} ({r[-1][key]/max(r[0][key],1):.1f}x)"
print(f"\n  ndof {r[0]['ndof']} -> {r[-1]['ndof']} ({r[-1]['ndof']/r[0]['ndof']:.1f}x):")
print(f"  point-Jacobi   : {grow('pointJac')}")
print(f"  two-level k=12 : {grow('twolvl_k12')}   (fixed small coarse)")
print(f"  two-level k=24 : {grow('twolvl_k24')}   (fixed larger coarse)")
print(f"  two-level k~N/6: {grow('twolvl_kprop')}   (coarse grows with N)")
print("\n  Verdict: if a FIXED-k two-level stays bounded as N grows -> a small coarse space fixes")
print("  the N-scaling -> H-AMS (H-matrix apply + sparse scalar-potential coarse) is worth building.")
print("  If only k~N bounds it -> the coarse space must grow -> closer to a full multilevel / H-LU.")
with open(os.path.join(HERE, "hdiv_demag_twolevel.json"), "w") as f:
    json.dump(res, f, indent=2)
print("\nsaved", os.path.join(HERE, "hdiv_demag_twolevel.json"))
