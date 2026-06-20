# -*- coding: utf-8 -*-
"""Twist check: is the Kelvin closure a p-method (exact-terminating at order=n)
rather than an h-method? Compare error-vs-DOF along the p-path (raise element
order on a fixed coarse mesh) vs the h-path (refine mesh at fixed order 1)."""
import sys
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_dtn_eigenvalue


def run(degree, order, maxh):
    r = kelvin_dtn_eigenvalue(R=1.0, degree=degree, maxh=maxh, order=order, dim=3)
    nd = None
    for k in ("ndof", "kelvin_ndof", "ndofs", "ndof_vol", "n_dof"):
        if k in r:
            nd = r[k]; break
    return r["rel_err"], nd, r


# probe return keys once
r0 = kelvin_dtn_eigenvalue(R=1.0, degree=2, maxh=0.5, order=1, dim=3)
print("KEYS:", list(r0.keys()))

for n in (2, 3):
    name = {2: "QUADRUPOLE n=2 (image = quadratic)", 3: "OCTUPOLE n=3 (image = cubic)"}[n]
    print(f"\n=== {name} ===")
    print("  p-path  (fixed coarse maxh=0.5, raise order):")
    for o in (1, 2, 3, 4):
        e, nd, _ = run(n, o, 0.5)
        print(f"    order={o}  maxh=0.50  ndof={nd}  rel_err={e:.3e}")
    print("  h-path  (fixed order=1, refine mesh):")
    for h in (0.6, 0.45, 0.32, 0.22, 0.16):
        e, nd, _ = run(n, 1, h)
        print(f"    order=1  maxh={h:.2f}  ndof={nd}  rel_err={e:.3e}")
print("\nDONE")
