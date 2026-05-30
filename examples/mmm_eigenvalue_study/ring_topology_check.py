#!/usr/bin/env python
"""Multiply-connected magnet: how much harder?  Build a RING (frame with a
central hole, first Betti number b1 = 1) and compare the true null dimension
(SVD) with what the LOCAL short-cycle (plaquette/edge-ring) basis spans.

Expectation (homology): null = (local contractible loops, im d2) (+) (b1 GLOBAL
non-contractible loops, H1). Local short-cycles span only the contractible part,
so they UNDER-span by exactly b1 -- the global "belt" loop(s) around the hole.

Usage: python ring_topology_check.py [outer] [thick] [nz] [mu_r]
  builds an (outer x outer) frame, 'thick' elements wide, nz layers (hole in middle).
"""
import sys
import numpy as np
import radia as rad

EDGE = 0.01
HEXF = [(0,1,2,3),(4,5,6,7),(0,1,5,4),(3,2,6,7),(0,3,7,4),(1,2,6,5)]


def build_frame(outer, thick, nz, mu_r):
    """Square frame: outer x outer grid minus the central (outer-2*thick)^2 hole."""
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs, cents = [], []
    lo, hi = thick, outer - thick  # hole = [lo,hi) x [lo,hi)
    for iz in range(nz):
        for iy in range(outer):
            for ix in range(outer):
                in_hole = (lo <= ix < hi) and (lo <= iy < hi)
                if in_hole:
                    continue
                x0, y0, z0 = ix*EDGE, iy*EDGE, iz*EDGE
                v = [[x0,y0,z0],[x0+EDGE,y0,z0],[x0+EDGE,y0+EDGE,z0],[x0,y0+EDGE,z0],
                     [x0,y0,z0+EDGE],[x0+EDGE,y0,z0+EDGE],[x0+EDGE,y0+EDGE,z0+EDGE],[x0,y0+EDGE,z0+EDGE]]
                o = rad.ObjHexahedron(v, [0,0,0]); rad.MatApl(o, mat); objs.append(o)
                cents.append(np.array(v))
    return rad.ObjCnt(objs), cents


def main():
    outer, thick, nz, mu_r = 6, 2, 1, 1e5
    a = sys.argv[1:]
    if len(a) >= 1: outer = int(a[0])
    if len(a) >= 2: thick = int(a[1])
    if len(a) >= 3: nz = int(a[2])
    if len(a) >= 4: mu_r = float(a[3])

    grp, verts = build_frame(outer, thick, nz, mu_r)
    ne = len(verts)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(grp)); N = np.asarray(N, float)
    null_dim = int(np.sum(np.linalg.svd(N, compute_uv=False) < 1e-9 * np.linalg.svd(N, compute_uv=False)[0]))

    # adjacency (face-centroid coincidence) + shared DOF (column cosine, axis-aligned OK)
    tol = 0.25 * EDGE
    def rkey(c): return tuple(np.round(c / tol).astype(int))
    fcent = [[verts[e][list(f)].mean(0) for f in HEXF] for e in range(ne)]
    bucket = {}
    for e in range(ne):
        for fc in fcent[e]:
            bucket.setdefault(rkey(fc), []).append(e)
    nbr = {e: set() for e in range(ne)}
    for k, es in bucket.items():
        if len(es) == 2:
            nbr[es[0]].add(es[1]); nbr[es[1]].add(es[0])
    # element-graph cycle rank = E - V + 1 (per connected comp)
    E = sum(len(v) for v in nbr.values()) // 2
    graph_cycle_rank = E - ne + 1

    cols = N / (np.linalg.norm(N, axis=0, keepdims=True) + 1e-30)
    sdc = {}
    def sdof(A, B):
        if (A, B) in sdc: return sdc[(A, B)]
        best, bi, bj = -2, -1, -1
        for i in range(6*A, 6*A+6):
            for j in range(6*B, 6*B+6):
                c = abs(float(cols[:, i] @ cols[:, j]))
                if c > best: best, bi, bj = c, i, j
        sdc[(A, B)] = (bi, bj); sdc[(B, A)] = (bj, bi); return bi, bj

    # short cycles (3 + 4)
    seen, Lcols = set(), []
    def add(order):
        kk = tuple(sorted(order))
        if kk in seen: return
        seen.add(kk)
        v = np.zeros(dof); k = len(order)
        for t in range(k):
            di, dj = sdof(order[t], order[(t+1) % k]); v[di]+=1; v[dj]-=1
        if np.linalg.norm(v) > 1e-9: Lcols.append(v)
    for e in range(ne):
        Ne = sorted(nbr[e])
        for ai in range(len(Ne)):
            for bi_ in range(ai+1, len(Ne)):
                A_, B_ = Ne[ai], Ne[bi_]
                if B_ in nbr[A_]: add((e, A_, B_))
                for d in (nbr[A_] & nbr[B_]) - {e}:
                    add((e, A_, d, B_))
    rankL = int(np.linalg.matrix_rank(np.column_stack(Lcols), tol=1e-9)) if Lcols else 0

    print(f"FRAME outer={outer} thick={thick} nz={nz}: {ne} elem, dof={dof}, mu_r={mu_r:.0e}")
    print(f"  graph cycle rank (F_int - N_elem + 1) = {graph_cycle_rank}")
    print(f"  SVD null_dim                          = {null_dim}")
    print(f"  short-cycle (local) span rank         = {rankL}")
    print(f"  DEFICIT (global belt loops, ~ b1)     = {null_dim - rankL}")
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
