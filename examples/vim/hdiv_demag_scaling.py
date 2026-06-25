"""
hdiv_demag_scaling.py -- the SWITCH-GATE scalability question for the HDiv-type demag operator.

hdiv_demag_quad_self.py established (3x3x3): loops EXACTLY field-null, plain Jacobi-MINRES is
mu_r-INDEPENDENT (regular 8, distorted ~75), distortion-sensitive.  The remaining question for
adopting HDiv-type over six-face surface-charge is SCALABILITY: how do plain-Jacobi iters grow with N (mesh
refinement) at fixed mu_r?  An integral (dense) operator's point-Jacobi iters often grow with N
-> that is exactly where a better preconditioner (block-Jacobi / H-ILU / AMS) earns its keep.

This measures Jacobi-MINRES iters vs ndof for n = 2..6 hexes/side, REGULAR and DISTORTED(0.12),
mu_r=1e4 fixed.  Verdict logic:
  - iters BOUNDED in N  -> plain Jacobi is N-scalable; no H-ILU needed even at scale.
  - iters GROW in N     -> plain Jacobi insufficient at scale; the H-ILU/AMS preconditioner is
                           justified -- but for CONDITIONING/scaling, NOT for the mu_r breakdown
                           (which the field-null loops already removed).
"""
import json, os
import numpy as np
from hdiv_demag_quad_self import build, jac_iters

HERE = os.path.dirname(os.path.abspath(__file__))
MUR = 1e4; CHI = MUR - 1.0

def measure(n, distort, nsub_hex=4, nsub_quad=6):
    d = build(distort, nsub_hex, nsub_quad, n=n)
    A = (1.0/CHI)*np.eye(d["ndof"]) - d["N"]
    itm, itg = jac_iters(A)
    return dict(n=n, ndof=int(d["ndof"]), n_loop=int(d["n_loop"]), loop_res=float(d["loop_res"]),
                iters=int(itm), cond=float(np.linalg.cond((1.0/CHI)*np.eye(d["rankQ"]) - d["Nsb"])),
                mu_min=float(d["mu"][0]))

res = {"mu_r": MUR, "regular": [], "distorted_0.12": []}
for tag, da in (("regular", 0.0), ("distorted_0.12", 0.12)):
    print(f"\n=== {tag}: Jacobi-MINRES iters vs N (mu_r={MUR:.0e}) ===")
    print(f"  {'n':>2} {'ndof':>5} {'loops':>5} {'loop_res':>9} {'iters':>6} {'cond(star)':>11} {'mu_min':>9}")
    for n in (2, 3, 4, 5, 6):
        r = measure(n, da)
        print(f"  {r['n']:>2} {r['ndof']:>5} {r['n_loop']:>5} {r['loop_res']:>9.1e} "
              f"{r['iters']:>6} {r['cond']:>11.2e} {r['mu_min']:>9.2e}")
        res[tag].append(r)

# verdict: is the iter count bounded in N?
for tag in ("regular", "distorted_0.12"):
    its = [r["iters"] for r in res[tag]]
    ndofs = [r["ndof"] for r in res[tag]]
    growth = its[-1] / max(its[0], 1)
    ndof_growth = ndofs[-1] / ndofs[0]
    print(f"\n[{tag}] iters {its[0]} -> {its[-1]} ({growth:.1f}x) as ndof {ndofs[0]} -> {ndofs[-1]} "
          f"({ndof_growth:.1f}x): {'BOUNDED -> plain Jacobi N-scalable' if growth < 1.6 else 'GROWS -> needs better precond at scale'}")

with open(os.path.join(HERE, "hdiv_demag_scaling.json"), "w") as f:
    json.dump(res, f, indent=2)
print("\nsaved", os.path.join(HERE, "hdiv_demag_scaling.json"))
