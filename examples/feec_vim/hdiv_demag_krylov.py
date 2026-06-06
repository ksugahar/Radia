"""
hdiv_demag_krylov.py -- BiCGSTAB vs GMRES vs MINRES on the HDiv-type demag operator, and the
SYMMETRY that decides the answer.

Key structural fact: the HDiv-type VIM operator is a GALERKIN energy form N = B^T G B with G the
symmetric Coulomb Gram (G_ab = G_ba, 1/r is symmetric) -> N is SYMMETRIC (PSD), so A=(1/chi)I-N
is SYMMETRIC INDEFINITE.  For a symmetric indefinite system the OPTIMAL Krylov method is MINRES
(3-term recurrence, O(1) memory, no restart, minimizes ||residual||_2).  GMRES is for GENERAL
matrices (long recurrence, needs restart -> the earlier 80000-iter artifact); BiCGSTAB is for
general matrices too (short recurrence but non-optimal, prone to breakdown on indefinite systems).

Contrast: the Radia MSC A_SS operator is COLLOCATION (EIEM2 point eval, asymmetric test/trial) ->
NON-symmetric (memory: ||A-A^T||/||A|| ~ 6.6e-2) -> MINRES is NOT applicable there, which is why
the shipped MSC solver uses BiCGSTAB (stable only with a strong H-LU precond) / GMRES.  The
HDiv-type's symmetry is therefore an ADVANTAGE: it unlocks MINRES (and the symmetric coarse-space
/ Calderon / AMS theory, most of which assumes symmetry).
"""
import json, os
import numpy as np
from scipy.sparse.linalg import minres, gmres, bicgstab, LinearOperator
from hdiv_demag_quad_self import build

HERE = os.path.dirname(os.path.abspath(__file__))

def count(solver, A, **kw):
    it = {"n": 0}
    cb = (lambda xk: it.__setitem__("n", it["n"]+1))
    n = A.shape[0]; b = np.ones(n)
    if solver == "minres":
        d = np.abs(np.diag(A)).copy(); d[d < 1e-30] = 1.0
        x, info = minres(A, b, M=LinearOperator((n, n), matvec=lambda y: y/d), rtol=1e-8, maxiter=4000, callback=cb)
    elif solver == "gmres":
        d = np.diag(A).copy(); d[np.abs(d) < 1e-30] = 1.0
        x, info = gmres(A, b, M=LinearOperator((n, n), matvec=lambda y: y/d), rtol=1e-8,
                        restart=n, maxiter=2, callback=cb, callback_type="pr_norm")
    else:  # bicgstab
        d = np.diag(A).copy(); d[np.abs(d) < 1e-30] = 1.0
        x, info = bicgstab(A, b, M=LinearOperator((n, n), matvec=lambda y: y/d), rtol=1e-8, maxiter=4000, callback=cb)
    rel = np.linalg.norm(b - A @ x) / np.linalg.norm(b)
    return it["n"], rel, info

print("=== HDiv demag: operator symmetry + Krylov method comparison (distorted 0.12) ===")
d = build(0.12, 4, 6, n=4)
N = d["N"]; ndof = d["ndof"]
asym = np.linalg.norm(N - N.T) / np.linalg.norm(N)
print(f"  ndof={ndof}: ||N - N^T||/||N|| = {asym:.2e}  (HDiv Galerkin -> SYMMETRIC; cf MSC collocation ~6.6e-2)")
res = {"ndof": int(ndof), "asym": float(asym), "by_mur": []}
print(f"\n  {'mu_r':>6} | {'MINRES':>16} | {'full-GMRES':>16} | {'BiCGSTAB':>16}   (iters / final rel.res)")
for mur in (1e2, 1e3, 1e4, 1e5):
    chi = mur - 1.0
    A = (1.0/chi)*np.eye(ndof) - N
    rows = {}
    for s in ("minres", "gmres", "bicgstab"):
        it, rel, info = count(s, A)
        rows[s] = (it, rel, info)
    print(f"  {mur:>6.0e} | {rows['minres'][0]:>5}  {rows['minres'][1]:>8.1e} | "
          f"{rows['gmres'][0]:>5}  {rows['gmres'][1]:>8.1e} | "
          f"{rows['bicgstab'][0]:>5}  {rows['bicgstab'][1]:>8.1e}")
    res["by_mur"].append({"mu_r": mur,
                          "minres": [int(rows['minres'][0]), float(rows['minres'][1])],
                          "gmres": [int(rows['gmres'][0]), float(rows['gmres'][1])],
                          "bicgstab": [int(rows['bicgstab'][0]), float(rows['bicgstab'][1])]})
print("\n  (BiCGSTAB 'iters' count matvecs differently; compare convergence + final residual.")
print("   For a SYMMETRIC indefinite A, MINRES is the principled choice -- short recurrence, O(1)")
print("   memory, no restart, no breakdown.  GMRES matches it but costs O(k) memory; BiCGSTAB can")
print("   stall/break on the harder mu_r.)")
with open(os.path.join(HERE, "hdiv_demag_krylov.json"), "w") as f:
    json.dump(res, f, indent=2)
print("\nsaved", os.path.join(HERE, "hdiv_demag_krylov.json"))
