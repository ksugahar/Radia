"""
hdiv_demag_opprecond.py -- does OPERATOR PRECONDITIONING (the real BEM/AMS mechanism) give a
pulse for the HDiv-type demag operator?  (Re-test, because hdiv_demag_twolevel.py tested the
WRONG mechanism -- low-dim coarse DEFLATION -- which is for operators with an isolated near-null
cluster, not for a single-layer-type operator.)

DIAGNOSIS first: the demag star eigenvalues accumulate at 0 (mu_min DROPS with N, mu_max ~const)
-- the SINGLE-LAYER signature (a smoothing, order -1 pseudodifferential operator).  For such an
operator the principled preconditioner is OPERATOR PRECONDITIONING (Hiptmair 2006 / Steinbach-
Wendland / Calderon): precondition with a SPARSE operator of the COMPLEMENTARY order.  In the
HDiv/charge picture the Coulomb Gram G ~ (-Laplacian)^{-1}, so the complementary operator is a
sparse FEM "Laplacian" on HDiv = the H(div) inner product  R = mass + div-div  (the Riesz map).

We deflate the loops EXACTLY (field-null, known) and precondition the STAR block
   A_star = (1/chi) I - Sstar^T N Sstar          (rankQ x rankQ, loops removed)
with the star-projected HDiv Riesz map
   R_star = Sstar^T (mass_hdiv + alpha*divdiv_hdiv) Sstar      (sparse-derived, SPD)
used as M^{-1} = R_star^{-1} (a cheap sparse SPD solve -- AMG-able in production).  GMRES iters
vs N (mu_r=1e4).  Pulse criterion: iters BOUNDED in N (vs the ~linear growth of point-Jacobi).
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

print("=== diagnosis: demag star eigenvalue distribution vs N (single-layer accumulation-at-0?) ===")
print(f"  {'n':>2} {'rankQ':>5} {'mu_max':>9} {'mu_min':>9} {'#(mu<1e-3)':>11} {'#(mu<1e-4)':>11}")
diag = []
cache = {}
for n in (3, 4, 5, 6):
    d = build(0.12, 4, 6, n=n); cache[n] = d
    mu = np.sort(np.abs(np.linalg.eigvals(d["Nsb"]).real))
    print(f"  {n:>2} {d['rankQ']:>5} {mu[-1]:>9.3e} {mu[0]:>9.3e} {int(np.sum(mu<1e-3)):>11} {int(np.sum(mu<1e-4)):>11}")
    diag.append({"n": n, "rankQ": int(d["rankQ"]), "mu_max": float(mu[-1]), "mu_min": float(mu[0]),
                 "n_below_1e-3": int(np.sum(mu < 1e-3))})

print(f"\n=== operator preconditioning on the CONSISTENT Galerkin operator A=(1/chi)M_mass - N ===")
print(f"    (the material term is the MASS matrix, not I -- I was the collocation/MSC form; a")
print(f"     Galerkin N=B^T G B must pair with (1/chi)M_mass).  star block, mu_r={MUR:.0e}.")
print(f"  {'n':>2} {'rankQ':>5} {'pointJac':>9} {'Riesz a=1':>10} {'Riesz a=10':>11} {'Riesz a=0.1':>12}")
rows = []
for n in (3, 4, 5, 6):
    d = cache[n]; S = d["Sstar"]; rq = d["rankQ"]
    Mst = S.T @ d["mass_hdiv"] @ S
    A_star = (1.0/CHI)*Mst - S.T @ d["N"] @ S          # CONSISTENT mass-weighted material term
    dpt = np.diag(A_star).copy(); dpt[np.abs(dpt) < 1e-30] = 1.0
    it_pt = gmres_iters(A_star, LinearOperator((rq, rq), matvec=lambda x: x/dpt))
    its = {}
    for a in (1.0, 10.0, 0.1):
        R = S.T @ (d["mass_hdiv"] + a*d["divdiv_hdiv"]) @ S
        Rinv = np.linalg.inv(R)
        its[a] = gmres_iters(A_star, LinearOperator((rq, rq), matvec=lambda x, Ri=Rinv: Ri @ x))
    print(f"  {n:>2} {rq:>5} {it_pt:>9} {its[1.0]:>10} {its[10.0]:>11} {its[0.1]:>12}")
    rows.append({"n": n, "rankQ": int(rq), "pointJac": it_pt,
                 "riesz_a1": its[1.0], "riesz_a10": its[10.0], "riesz_a0.1": its[0.1]})

def gr(key): return f"{rows[0][key]} -> {rows[-1][key]} ({rows[-1][key]/max(rows[0][key],1):.1f}x)"
print(f"\n  ndof(star) {rows[0]['rankQ']} -> {rows[-1]['rankQ']} ({rows[-1]['rankQ']/rows[0]['rankQ']:.1f}x):")
print(f"  point-Jacobi : {gr('pointJac')}")
print(f"  Riesz a=1    : {gr('riesz_a1')}")
print(f"  Riesz a=10   : {gr('riesz_a10')}")
print(f"  Riesz a=0.1  : {gr('riesz_a0.1')}")
print("  PULSE if a Riesz-map (sparse SPD, AMG-able) growth << point-Jacobi growth -> operator")
print("  preconditioning works -> H-AMS (H-matrix apply + sparse Riesz/AMS auxiliary) is viable.")
with open(os.path.join(HERE, "hdiv_demag_opprecond.json"), "w") as f:
    json.dump({"diag": diag, "opprecond": rows, "mu_r": MUR}, f, indent=2)
print("\nsaved", os.path.join(HERE, "hdiv_demag_opprecond.json"))
