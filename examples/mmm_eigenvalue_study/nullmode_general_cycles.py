#!/usr/bin/env python
"""General (mesh-topology-free) cycle basis: build the loop null modes from the
element-adjacency graph's SHORT CYCLES (3- and 4-cycles), WITHOUT assuming a
grid index.  Adjacency = elements sharing a face (face-centroid coincidence);
shared-face DOF = N-column parallelism.  This is the algorithm to port to C++
(where adjacency + shared DOF come directly from the face geometry).

Verify on a structured grid AND a sheared (distorted) grid: the short-cycle
loop set must (a) be exact null modes (||N L||~0) and (b) span the null space
(rank(L) == null_dim).

Usage: python nullmode_general_cycles.py [nx ny nz] [mu_r] [shear]
"""
import sys
import numpy as np
import radia as rad

EDGE = 0.01
# hex face vertex groups for my vertex order
HEXF = [(0,1,2,3),(4,5,6,7),(0,1,5,4),(3,2,6,7),(0,3,7,4),(1,2,6,5)]


def build(nx, ny, nz, mu_r, shear=0.0):
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs, verts_all = [], []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                x0, y0, z0 = ix*EDGE, iy*EDGE, iz*EDGE
                base = [(x0,y0,z0),(x0+EDGE,y0,z0),(x0+EDGE,y0+EDGE,z0),(x0,y0+EDGE,z0),
                        (x0,y0,z0+EDGE),(x0+EDGE,y0,z0+EDGE),(x0+EDGE,y0+EDGE,z0+EDGE),(x0,y0+EDGE,z0+EDGE)]
                # shear: x += shear * z  (parallelepiped, non-axis-aligned faces)
                v = [[vx + shear*vz, vy, vz] for (vx, vy, vz) in base]
                o = rad.ObjHexahedron(v, [0,0,0]); rad.MatApl(o, mat); objs.append(o)
                verts_all.append(np.array(v))
    return rad.ObjCnt(objs), verts_all


def face_centroids(verts):
    return [verts[list(f)].mean(0) for f in HEXF]


def main():
    nx, ny, nz = 4, 4, 2
    mu_r, shear = 1e5, 0.0
    a = sys.argv[1:]
    if len(a) >= 3: nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4: mu_r = float(a[3])
    if len(a) >= 5: shear = float(a[4])

    grp, verts_all = build(nx, ny, nz, mu_r, shear)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(grp)); N = np.asarray(N, float)
    ne = len(verts_all)
    s = np.linalg.svd(N, compute_uv=False)
    null_dim = int(np.sum(s < 1e-9*s[0]))
    print(f"{nx}x{ny}x{nz}={ne} elem dof={dof} mu_r={mu_r:.0e} shear={shear} null_dim={null_dim}")

    # adjacency by face-centroid coincidence (geometric, no grid index)
    key = {}
    tol = 0.25 * EDGE
    def rkey(c): return tuple(np.round(c / tol).astype(int))
    adj = {}            # (ea,eb) -> True
    fcent = [face_centroids(v) for v in verts_all]
    bucket = {}
    for e in range(ne):
        for fc in fcent[e]:
            bucket.setdefault(rkey(fc), []).append(e)
    for k, elems in bucket.items():
        if len(elems) == 2:
            a_, b_ = elems
            adj[(a_, b_)] = True; adj[(b_, a_)] = True
    # neighbor list
    nbr = {e: set() for e in range(ne)}
    for (a_, b_) in adj:
        nbr[a_].add(b_)

    # shared-face DOF via N-column parallelism
    cols = N / (np.linalg.norm(N, axis=0, keepdims=True) + 1e-30)
    sdof_cache = {}
    def sdof(A, B):
        if (A, B) in sdof_cache: return sdof_cache[(A, B)]
        best, bi, bj = -2.0, -1, -1
        for i in range(6*A, 6*A+6):
            for j in range(6*B, 6*B+6):
                c = abs(float(cols[:, i] @ cols[:, j]))
                if c > best: best, bi, bj = c, i, j
        sdof_cache[(A, B)] = (bi, bj); sdof_cache[(B, A)] = (bj, bi)
        return bi, bj

    # enumerate short cycles (3- and 4-cycles) in the adjacency graph
    cycles = set()
    for e in range(ne):
        Ne = sorted(nbr[e])
        for ai in range(len(Ne)):
            a_ = Ne[ai]
            for bi_ in range(ai+1, len(Ne)):
                b_ = Ne[bi_]
                if b_ in nbr[a_]:                     # 3-cycle e-a-b
                    cycles.add(tuple(sorted((e, a_, b_))) + ("3",))
                # 4-cycle e-a-d-b: common neighbor d of a_ and b_, d!=e
                common = (nbr[a_] & nbr[b_]) - {e}
                for d in common:
                    cyc = (e, a_, d, b_)
                    cycles.add(tuple(sorted(cyc)) + ("4",))

    # build loop vectors (need an ORDERED cycle; recover order from adjacency)
    def order_cycle(nodes):
        nodes = list(nodes)
        # greedy walk through adjacency
        path = [nodes[0]]; rem = set(nodes[1:])
        while rem:
            nxt = next((x for x in rem if x in nbr[path[-1]]), None)
            if nxt is None: return None
            path.append(nxt); rem.discard(nxt)
        if path[0] not in nbr[path[-1]]: return None     # not closed
        return path

    Lcols = []
    for c in cycles:
        nodes = c[:-1]
        order = order_cycle(nodes)
        if order is None: continue
        v = np.zeros(dof); k = len(order)
        for t in range(k):
            A, B = order[t], order[(t+1) % k]
            di, dj = sdof(A, B); v[di] += 1.0; v[dj] -= 1.0
        if np.linalg.norm(v) > 1e-9: Lcols.append(v)
    L = np.column_stack(Lcols)

    resid = np.linalg.norm(N @ L, axis=0) / (np.linalg.norm(L, axis=0)+1e-30)
    rankL = int(np.linalg.matrix_rank(L, tol=1e-9))
    print(f"short cycles: {len(Lcols)} loops (3- & 4-cycles)")
    print(f"null-mode residual ||N*loop||/||loop||: max={resid.max():.2e}")
    print(f"rank(L)={rankL}  null_dim={null_dim}  SPAN={'OK' if rankL==null_dim else 'FAIL'}")
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
