"""
hdiv_demag_blockjacobi.py -- can a CHEAP block-Jacobi rescue the "Jacobi de sumu" hope at scale?

hdiv_demag_scaling.py showed plain POINT-Jacobi-MINRES iters GROW with N (distorted 0.12:
20->61->130->254->431 over ndof 36->756; ~linear).  mu_r-independence is already won (field-null
loops).  Question: is the N-growth killed by a cheap BLOCK-Jacobi -- additive Schwarz with
ELEMENT patches (each hex's 6 RT0 face-dofs), local 6x6 solves -- WITHOUT any H-LU factorization?
If yes, "Jacobi kurai de" genuinely holds (the remaining work is a cheap block precond, not H-ILU).

M^{-1} = sum_e R_e^T (R_e A R_e^T)^{-1} R_e   (overlapping element patches).  Compared head-to-head
with point-Jacobi (and a naive additive spectral-coarse correction) across n=2..6, full GMRES.

RESULT (2026-06-07, honest -- the hoped-for "cheap block Jacobi rescues scale" did NOT pan out):
  ndof 36->756 (distorted 0.12, mu_r=1e4):
    point-Jacobi    26 -> 283   (grows ~linearly)
    element-block   27 -> 449   (WORSE than point-Jacobi -- local patches do not touch the
                                 long-range coupling that drives the growth)
    naive coarse    30 -> 304   (additive smallest-|lambda| correction: no bound)
  => the N-growth is LONG-RANGE / far-field dominated; a cheap LOCAL block does not help, and a
     naive additive coarse correction does not either.  A proper two-level/deflated/GenEO design
     (or the shipped H-ILU far-field low-rank) is the OPEN scalable preconditioner -- but this is
     a CONDITIONING/N-scaling matter, ORTHOGONAL to the mu_r breakdown (which the field-null HDiv
     loops already removed: plain Jacobi is mu_r-independent).
"""
import json, os
import numpy as np
from scipy.sparse.linalg import gmres, LinearOperator
from hdiv_demag_quad_self import build

HERE = os.path.dirname(os.path.abspath(__file__))
MUR = 1e4; CHI = MUR - 1.0

def point_jacobi_M(A):
    d = np.abs(np.diag(A)).copy(); d[d < 1e-30] = 1.0
    return LinearOperator(A.shape, matvec=lambda x: x/d)

def block_jacobi_M(A, patches):
    """Additive Schwarz with the given dof patches (overlapping).  Precompute each local inverse."""
    inv = []
    for p in patches:
        Ap = A[np.ix_(p, p)]
        inv.append((p, np.linalg.inv(Ap)))
    n = A.shape[0]
    def apply(x):
        y = np.zeros(n)
        for p, Ai in inv:
            y[p] += Ai @ x[p]
        return y
    return LinearOperator((n, n), matvec=apply)

def iters(A, M):
    """FULL GMRES (restart=ndof, no restart) -> honest Krylov count; handles indefinite precond
    (block-Jacobi of an indefinite A has indefinite local blocks, which MINRES would reject)."""
    n = A.shape[0]; it = {"n": 0}
    gmres(A, np.ones(n), M=M, rtol=1e-8, restart=n, maxiter=2,
          callback=lambda xk: it.__setitem__("n", it["n"]+1), callback_type="pr_norm")
    return it["n"]

def coarse_deflated_M(A, k):
    """point-Jacobi + an IDEAL coarse correction: W = the k smallest-|lambda| eigvecs of A,
    M^{-1} = D^{-1} + W (W^T A W)^{-1} W^T.  Tests whether a GLOBAL coarse space (unlike the LOCAL
    element block) bounds the N-growth -> confirms the missing ingredient is a coarse space."""
    w, V = np.linalg.eigh(A)
    idx = np.argsort(np.abs(w))[:k]; W = V[:, idx]
    Ac_inv = np.linalg.inv(W.T @ A @ W)
    d = np.abs(np.diag(A)).copy(); d[d < 1e-30] = 1.0
    return LinearOperator(A.shape, matvec=lambda x: x/d + W @ (Ac_inv @ (W.T @ x)))

res = {"mu_r": MUR, "distorted_0.12": []}
print(f"=== distorted 0.12, mu_r={MUR:.0e}: point-Jac vs element-block-Jac vs coarse-deflated (full GMRES) ===")
print(f"  {'n':>2} {'ndof':>5} {'point-Jac':>10} {'block-Jac':>10} {'coarse k=12':>12} {'coarse k~ndof/8':>15}")
for n in (2, 3, 4, 5, 6):
    d = build(0.12, 4, 6, n=n)
    A = (1.0/CHI)*np.eye(d["ndof"]) - d["N"]
    it_pt = iters(A, point_jacobi_M(A))
    it_bl = iters(A, block_jacobi_M(A, d["el_dofs"]))
    it_c1 = iters(A, coarse_deflated_M(A, min(12, d["ndof"]-1)))
    it_c2 = iters(A, coarse_deflated_M(A, max(1, d["ndof"]//8)))
    print(f"  {n:>2} {d['ndof']:>5} {it_pt:>10} {it_bl:>10} {it_c1:>12} {it_c2:>15}")
    res["distorted_0.12"].append({"n": n, "ndof": int(d["ndof"]), "point_jacobi": it_pt,
                                  "block_jacobi": it_bl, "coarse_k12": it_c1, "coarse_kprop": it_c2})

pj = [r["point_jacobi"] for r in res["distorted_0.12"]]
bj = [r["block_jacobi"] for r in res["distorted_0.12"]]
c1 = [r["coarse_k12"] for r in res["distorted_0.12"]]
c2 = [r["coarse_kprop"] for r in res["distorted_0.12"]]
nd = [r["ndof"] for r in res["distorted_0.12"]]
print(f"\n  ndof {nd[0]}->{nd[-1]} ({nd[-1]/nd[0]:.1f}x):")
print(f"  point-Jacobi   : {pj[0]:>4} -> {pj[-1]:>4} ({pj[-1]/pj[0]:.1f}x)  -- grows ~linearly")
print(f"  block-Jacobi   : {bj[0]:>4} -> {bj[-1]:>4} ({bj[-1]/bj[0]:.1f}x)  -- LOCAL patch: NO help (worse)")
print(f"  coarse k=12    : {c1[0]:>4} -> {c1[-1]:>4} ({c1[-1]/c1[0]:.1f}x)  -- naive additive coarse")
print(f"  coarse k~ndof/8: {c2[0]:>4} -> {c2[-1]:>4} ({c2[-1]/c2[0]:.1f}x)  -- naive additive coarse")
print("  HONEST verdict (verify-first; a hoped-for 'coarse fixes it' was NOT borne out):")
print("   - LOCAL element-block Schwarz does NOT help -> the N-growth is LONG-RANGE / far-field")
print("     dominated (robust finding: local patches cannot fix a dense integral operator).")
print("   - a NAIVE additive smallest-|lambda| coarse correction also does NOT bound it (the")
print("     additive D^-1 + W(W^TAW)^-1 W^T design is poor for indefinite A; a proper two-level /")
print("     deflated-GMRES / GenEO design is the OPEN task -- same conclusion as the MSC A_SS arc).")
print("   => plain Jacobi wins the mu_r axis (breakdown gone); the N-scaling preconditioner is the")
print("      open part.  The shipped H-ILU is the empirical answer; a cheap clean one is unfound.")
with open(os.path.join(HERE, "hdiv_demag_blockjacobi.json"), "w") as f:
    json.dump(res, f, indent=2)
print("\nsaved", os.path.join(HERE, "hdiv_demag_blockjacobi.json"))
