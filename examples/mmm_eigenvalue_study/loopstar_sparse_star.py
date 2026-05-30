#!/usr/bin/env python
"""Verify the SPARSE star basis (what the C++ HACApK sandwich will use) spans the
correct complement of the loop null space -- i.e. gives the SAME field-exact +
well-conditioned restriction as the dense SVD star (loopstar_prototype.py), but
from a SPARSE, locally-constructed basis.

Sparse star basis (MSC sigma-DOF structure; dof = 6*e + f for hexes):
  (1) BOUNDARY face DOFs (not in any shared pair)        -> unit vector e_i
  (2) SYMMETRIC internal:  sigma_A + sigma_B per shared face -> (+1@A, +1@B)
  (3) DIVERGENCE (cut) per element: sum over its shared faces of (+1@own, -1@nbr)
      -- the anti-symmetric cut space; drop one per connected comp (redundancy).
Dimension check: (1)+(2)+(3) = F_b + E + (N_elem-1) = 7*N_elem - E - 1 = star dim.

We QR the sparse columns to an orthonormal span, then run the SAME restriction
test as the dense prototype. If field-dev ~ eps and cond(A_star) is small for
the LINEAR case, the sparse construction is correct and ready to port to C++.

Cube only (axis-aligned), where the column-cosine shared-DOF id is reliable.
Usage: python loopstar_sparse_star.py [n]   (cube n^3, default 4)
"""
import sys
import numpy as np
import radia as rad

EDGE = 0.01
HEXF = [(0,1,2,3),(4,5,6,7),(0,1,5,4),(3,2,6,7),(0,3,7,4),(1,2,6,5)]


def build_cube(n):
    rad.UtiDelAll()
    mat = rad.MatLin(1000.0)
    objs, verts = [], []
    for iz in range(n):
        for iy in range(n):
            for ix in range(n):
                x0, y0, z0 = ix*EDGE, iy*EDGE, iz*EDGE
                v = [[x0,y0,z0],[x0+EDGE,y0,z0],[x0+EDGE,y0+EDGE,z0],[x0,y0+EDGE,z0],
                     [x0,y0,z0+EDGE],[x0+EDGE,y0,z0+EDGE],[x0+EDGE,y0+EDGE,z0+EDGE],[x0,y0+EDGE,z0+EDGE]]
                o = rad.ObjHexahedron(v, [0,0,0]); rad.MatApl(o, mat)
                objs.append(o); verts.append(np.array(v, float))
    return rad.ObjCnt(objs), verts


def adjacency(verts, N):
    """Element adjacency by coincident face centroids; shared-DOF pair by max
    column-cosine of N (robust for axis-aligned)."""
    ne = len(verts)
    tol = 0.25*EDGE
    def rkey(c): return tuple(np.round(c/tol).astype(int))
    fcent = [[verts[e][list(f)].mean(0) for f in HEXF] for e in range(ne)]
    bucket = {}
    for e in range(ne):
        for fc in fcent[e]:
            bucket.setdefault(rkey(fc), []).append(e)
    nbr = {e: set() for e in range(ne)}
    for es in bucket.values():
        if len(es) == 2:
            nbr[es[0]].add(es[1]); nbr[es[1]].add(es[0])
    cols = N / (np.linalg.norm(N, axis=0, keepdims=True) + 1e-30)
    sdof = {}
    for e in range(ne):
        for e2 in nbr[e]:
            if (e, e2) in sdof: continue
            best, bi, bj = -2, -1, -1
            for i in range(6*e, 6*e+6):
                for j in range(6*e2, 6*e2+6):
                    c = abs(float(cols[:, i] @ cols[:, j]))
                    if c > best: best, bi, bj = c, i, j
            sdof[(e, e2)] = (bi, bj); sdof[(e2, e)] = (bj, bi)
    return nbr, sdof


def build_sparse_star(ne, dof, nbr, sdof):
    cols = []
    internal = set()
    for (a, b), (da, db) in sdof.items():
        internal.add(da); internal.add(db)
    # (1) boundary unit vectors
    for d in range(dof):
        if d not in internal:
            v = np.zeros(dof); v[d] = 1.0; cols.append(v)
    # (2) symmetric per internal face (count each undirected face once)
    seen = set()
    for (a, b), (da, db) in sdof.items():
        key = tuple(sorted((da, db)))
        if key in seen: continue
        seen.add(key)
        v = np.zeros(dof); v[da] = 1.0; v[db] = 1.0; cols.append(v)
    # (3) divergence per element (drop element 0 to kill the 1 redundancy)
    for e in range(1, ne):
        v = np.zeros(dof)
        for e2 in nbr[e]:
            da, db = sdof[(e, e2)]
            v[da] += 1.0; v[db] -= 1.0
        if np.linalg.norm(v) > 1e-12: cols.append(v)
    return np.column_stack(cols) if cols else np.zeros((dof, 0))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    grp, verts = build_cube(n)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(grp)); N = np.asarray(N, float)
    s = np.linalg.svd(N, compute_uv=False); null_dim = int((s < 1e-9*s[0]).sum())
    star_dim = dof - null_dim

    nbr, sdof = adjacency(verts, N)
    T = build_sparse_star(len(verts), dof, nbr, sdof)
    rankT = int(np.linalg.matrix_rank(T, tol=1e-9))
    # orthonormal span of the sparse columns
    Q, _ = np.linalg.qr(T)
    U_s = Q[:, :rankT]

    # how well does span(T) avoid the loops?  (component of loops inside span(T))
    Vt = np.linalg.svd(N)[2]; L = Vt.T[:, s < 1e-9*s[0]]
    leak = np.linalg.norm(U_s.T @ L)  # ~0 if span(T) is a clean complement-ish

    print(f"cube {n}^3: dof={dof}, null={null_dim}, star_dim={star_dim}")
    print(f"sparse T: {T.shape[1]} cols, rank={rankT}  (want rank == star_dim={star_dim})")
    print(f"loop leakage into span(T)  ||U_s^T L|| = {leak:.2e}  (small = clean complement)")

    rng = np.random.default_rng(0); b = rng.standard_normal(dof)
    chi = np.full(dof, 1e5 - 1.0)
    A = np.diag(1.0/chi) - N
    sig_full = np.linalg.solve(A, b)
    Ar = U_s.T @ A @ U_s
    sig_r = U_s @ np.linalg.solve(Ar, U_s.T @ b)
    fd = np.linalg.norm(N@sig_r - N@sig_full)/(np.linalg.norm(N@sig_full)+1e-30)
    print(f"\nLINEAR mu=1e5 (sparse star):  cond(A)={np.linalg.cond(A):.2e}  "
          f"cond(A_star)={np.linalg.cond(Ar):.2e}  field dev={fd:.2e}")
    print("PASS if rank==star_dim, leakage small, field dev ~ eps, cond(A_star) small.")
    rad.UtiDelAll()


if __name__ == '__main__':
    main()
