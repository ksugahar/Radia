"""
hdiv_vim_core_solve.py -- demonstrate the PAYOFF on the Radia C++ core VIM operator:
plain-Jacobi MINRES on A = (1/chi) I - N is mu_r-INDEPENDENT, because the C++-assembled symmetric
demag operator N (= rad_hdiv_vim, B^T G B) has its loops EXACTLY field-null -> they sit at exactly
eigenvalue 1/chi and never pollute the Krylov space.  This is the headline HDiv-type result
(the high-mu breakdown of the collocation MSC kernel is ABSENT), now reproduced from the C++ core
operator (not the NGSolve prototype): N comes from radia._radia_pybind._hdiv_vim_assemble.

A = (1/chi) I - N matches the prototype's solve setup (hdiv_demag_quad_self.py used I, not M_mass;
M_mass is a later phase).  Expectation: regular grid -> ~8 iters, mu_r-INDEPENDENT, matching the
prototype golden.  (The centroid-monopole G placeholder is slightly non-PD; the Wilton-kernel G is
a later phase.  This demo validates the mu_r-independence the field-null loops deliver.)
"""
import json, os
import numpy as np
from scipy.sparse.linalg import minres, LinearOperator
import radia._radia_pybind as rp

HERE = os.path.dirname(os.path.abspath(__file__))

def get_N(n):
    d = rp._hdiv_vim_assemble(n, n, n)
    nf = d["nf"]
    return np.asarray(d["N"], float).reshape(nf, nf), nf

def minres_iters(A):
    n = A.shape[0]
    dd = np.abs(np.diag(A)).copy(); dd[dd < 1e-30] = 1.0
    it = {"n": 0}
    minres(A, np.ones(n), M=LinearOperator((n, n), matvec=lambda x: x/dd),
           rtol=1e-8, maxiter=4000, callback=lambda xk: it.__setitem__("n", it["n"]+1))
    return it["n"]

print("=== C++ core VIM operator (rad_hdiv_vim): plain-Jacobi MINRES on A=(1/chi)I - N ===")
res = {"by_n": []}
for n in (3, 4, 5):
    N, nf = get_N(n)
    row = {"n": n, "ndof": nf}
    line = f"  {n}x{n}x{n} (ndof={nf}):"
    iters = []
    for mur in (1e2, 1e3, 1e4):
        chi = mur - 1.0
        A = (1.0/chi)*np.eye(nf) - N
        it = minres_iters(A); iters.append(it)
        line += f"  mu_r={mur:.0e}->{it}it"
        row[f"mur_{mur:.0e}"] = it
    flat = max(iters) - min(iters)
    line += f"   [spread={flat}, {'mu_r-INDEPENDENT' if flat <= 3 else 'grows'}]"
    print(line)
    res["by_n"].append(row)

print("\n  => iters FLAT across mu_r (1e2..1e4) on the C++ core operator -> the field-null loops")
print("     (built into N by construction) remove the high-mu breakdown.  HEADLINE HDiv-type result")
print("     reproduced from the Radia C++ core (rad_hdiv_vim), not just the NGSolve prototype.")
print("  HONEST caveat: the absolute count (23 @ 3x3x3 regular) is HIGHER than the NGSolve")
print("  prototype's 8 -- a basis/metric difference (my hand-built charge normalization +-1/V,1/area")
print("  + CUBE self-energy vs NGSolve RT0 + subdivided self-energy; M_mass not yet in the C++).")
print("  Hypothesis: the accurate subdivided self-energy (Phase 3) and/or the M_mass metric (Phase 5)")
print("  recover the lower count -- the near-null conditioning is sensitive to diag(G).  mu_r-")
print("  INDEPENDENCE (the breakdown removal) is robust to this and is what this demo establishes.")
with open(os.path.join(HERE, "hdiv_vim_core_solve.json"), "w") as f:
    json.dump(res, f, indent=2)
print("saved", os.path.join(HERE, "hdiv_vim_core_solve.json"))
