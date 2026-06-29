"""Report condition number distribution of P2 Vandermonde across sphere mesh.

For each triangle in the dedup'd sphere mesh, compute cond(V) of the 6x6
{1, s, z, s^2, sz, z^2} Vandermonde. Identify triangles where cond > 1e12
(threshold for FP64 inversion accuracy). These are the ones likely producing
garbage element matrices in the C++ NGSolve assembly (which uses CalcInverse
without re-scaling).
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
sys.path.insert(0, ".")
from test_p2_dedup_mesh import build_dedup_mesh


def vandermonde(rn, zn):
    sn = rn * rn
    V = np.zeros((6, 6))
    V[:, 0] = 1.0
    V[:, 1] = sn
    V[:, 2] = zn
    V[:, 3] = sn * sn
    V[:, 4] = sn * zn
    V[:, 5] = zn * zn
    return V


def main():
    nodes, triangles, is_cond, _ = build_dedup_mesh(4, 4, 12)
    print(f"Sphere (4, 4, 12) dedup'd: {len(triangles)} triangles")

    conds = []
    bad = []
    for k, tri in enumerate(triangles):
        rn = nodes[tri, 0]
        zn = nodes[tri, 1]
        V = vandermonde(rn, zn)
        c = np.linalg.cond(V)
        conds.append(c)
        if c > 1e12 or not np.isfinite(c):
            bad.append((k, c, tri.tolist(), rn.tolist(), zn.tolist()))

    conds = np.array(conds)
    print(f"\nVandermonde cond stats:")
    print(f"  min:  {np.nanmin(conds):.2e}")
    print(f"  median: {np.nanmedian(conds):.2e}")
    print(f"  90%:  {np.nanpercentile(conds, 90):.2e}")
    print(f"  99%:  {np.nanpercentile(conds, 99):.2e}")
    print(f"  max:  {np.nanmax(conds):.2e}")
    print(f"\n  Cond > 1e8:   {np.sum(conds > 1e8)} triangles")
    print(f"  Cond > 1e10:  {np.sum(conds > 1e10)} triangles")
    print(f"  Cond > 1e12:  {np.sum(conds > 1e12)} triangles")
    print(f"  Cond > 1e14:  {np.sum(conds > 1e14)} triangles")
    print(f"  Cond > 1e16:  {np.sum(conds > 1e16)} triangles (FP64 broken)")

    print(f"\nWorst {min(5, len(bad))} bad triangles (cond > 1e12):")
    bad.sort(key=lambda x: -x[1])
    for k, c, tri, rn, zn in bad[:5]:
        print(f"  tri[{k}] cond={c:.2e}")
        print(f"    rn={rn}")
        print(f"    zn={zn}")


if __name__ == "__main__":
    main()
