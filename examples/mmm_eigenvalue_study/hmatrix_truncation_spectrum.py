#!/usr/bin/env python
"""Phase A: does low-rank (H-matrix-like) truncation regularize the near-zero
spectrum of the MSC system matrix?

Hypothesis (Sugahara / Ida, 2026-05): the H-matrix ACA approximation perturbs
the exact-zero eigenvalues of the geometric operator G, pushing the smallest
eigenvalues of A = diag(1/(mu_r-1)) - G away from zero, which would explain why
HACApK BiCGSTAB needs fewer iterations and avoids the loop-mode contamination.
(ksugahar himself flagged it might instead get WORSE -- we report honestly.)

This is a CONTROLLED proxy for HACApK done purely in Python (no rebuild):
  1. extract exact dense G from Radia.
  2. geometric-bisection cluster tree over ELEMENT centroids (each element owns
     6 consecutive face DOFs).
  3. standard admissibility min(diam_a,diam_b) <= eta * dist(box_a,box_b);
     admissible far blocks -> truncated SVD (rel. tol eps); near blocks exact.
  4. densify -> G_approx; A_* = diag(1/(mu_r-1)) - G_*.
  5. compare smallest |eig| and condition number vs eps.

Note: sign of G does not affect the near-zero-mode conclusion (G v = 0 for the
null modes, so A v = (1/(mu_r-1)) v regardless).

Usage: python hmatrix_truncation_spectrum.py [n] [mu_r]   (default n=6, mu_r=5000)
"""
import sys
import numpy as np
import radia as rad

EDGE = 0.01
ETA = 2.0
LEAF_ELEMS = 8


def build_cube(n, mu_r, edge=EDGE):
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs, cents = [], []
    for iz in range(n):
        for iy in range(n):
            for ix in range(n):
                x0, y0, z0 = ix * edge, iy * edge, iz * edge
                verts = [
                    [x0, y0, z0], [x0 + edge, y0, z0],
                    [x0 + edge, y0 + edge, z0], [x0, y0 + edge, z0],
                    [x0, y0, z0 + edge], [x0 + edge, y0, z0 + edge],
                    [x0 + edge, y0 + edge, z0 + edge], [x0, y0 + edge, z0 + edge],
                ]
                o = rad.ObjHexahedron(verts, [0, 0, 0])
                rad.MatApl(o, mat)
                objs.append(o)
                cents.append([(ix + 0.5) * edge, (iy + 0.5) * edge, (iz + 0.5) * edge])
    return rad.ObjCnt(objs), np.asarray(cents)


class Node:
    __slots__ = ("idx", "lo", "hi", "child")

    def __init__(self, idx, cents):
        self.idx = idx
        p = cents[idx]
        self.lo = p.min(0)
        self.hi = p.max(0)
        self.child = None


def build_tree(idx, cents, leaf):
    node = Node(idx, cents)
    if len(idx) <= leaf:
        return node
    p = cents[idx]
    ax = int(np.argmax(node.hi - node.lo))
    med = np.median(p[:, ax])
    left = idx[p[:, ax] <= med]
    right = idx[p[:, ax] > med]
    if len(left) == 0 or len(right) == 0:
        order = idx[np.argsort(p[:, ax])]
        half = len(idx) // 2
        left, right = order[:half], order[half:]
    node.child = [build_tree(left, cents, leaf), build_tree(right, cents, leaf)]
    return node


def diam(node):
    return float(np.linalg.norm(node.hi - node.lo))


def box_dist(a, b):
    d = np.maximum(0.0, np.maximum(a.lo - b.hi, b.lo - a.hi))
    return float(np.linalg.norm(d))


def collect(a, b, eta, near, far):
    if box_dist(a, b) > 0 and min(diam(a), diam(b)) <= eta * box_dist(a, b):
        far.append((a, b))
        return
    if a.child is None and b.child is None:
        near.append((a, b))
        return
    A = a.child if a.child else [a]
    B = b.child if b.child else [b]
    for ca in A:
        for cb in B:
            collect(ca, cb, eta, near, far)


def dof_rows(idx):
    return np.concatenate([np.arange(6 * e, 6 * e + 6) for e in idx])


def truncate(G, near, far, eps):
    Ga = np.zeros_like(G)
    stored = 0
    for (a, b) in near:
        ri, cj = dof_rows(a.idx), dof_rows(b.idx)
        blk = G[np.ix_(ri, cj)]
        Ga[np.ix_(ri, cj)] = blk
        stored += blk.size
    for (a, b) in far:
        ri, cj = dof_rows(a.idx), dof_rows(b.idx)
        blk = G[np.ix_(ri, cj)]
        U, s, Vt = np.linalg.svd(blk, full_matrices=False)
        k = int(np.sum(s > eps * s[0])) if s[0] > 0 else 0
        k = max(1, k)
        Ga[np.ix_(ri, cj)] = (U[:, :k] * s[:k]) @ Vt[:k]
        stored += k * (blk.shape[0] + blk.shape[1])
    return Ga, stored / G.size


def smallest(A, k=6):
    w = np.linalg.eigvals(A)
    return np.sort(np.abs(w))[:k], w


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    mu_r = float(sys.argv[2]) if len(sys.argv) > 2 else 5000.0
    inv_chi = 1.0 / (mu_r - 1.0)

    block, cents = build_cube(n, mu_r)
    n_elem = len(cents)
    G, dof = rad.GetInteractMatrix(rad.BuildMatrix(block))
    G = np.asarray(G, dtype=float)
    print(f"cube {n}x{n}x{n} = {n_elem} elem, dof={dof}, mu_r={mu_r}, "
          f"1/(mu_r-1)={inv_chi:.4e}")

    root = build_tree(np.arange(n_elem), cents, LEAF_ELEMS)
    near, far = [], []
    collect(root, root, ETA, near, far)
    print(f"cluster tree: leaf<= {LEAF_ELEMS} elem, eta={ETA}; "
          f"near blocks={len(near)}, far(admissible) blocks={len(far)}")

    A_exact = np.diag(np.full(dof, inv_chi)) - G
    s_ex, w_ex = smallest(A_exact)
    cond_ex = np.abs(w_ex).max() / np.abs(w_ex).min()
    print("\nEXACT A:")
    print("  smallest |lambda|:", np.array2string(s_ex, precision=4))
    print(f"  |lambda|_max={np.abs(w_ex).max():.4e}  cond={cond_ex:.4e}")

    print("\nTRUNCATED A (vary eps):")
    print(f"  {'eps':>8} {'comp':>7} {'min|lam|':>12} {'max|lam|':>12} {'cond':>11} "
          f"{'min/floor':>10}")
    for eps in (1e-2, 1e-3, 1e-4, 1e-6):
        Ga, comp = truncate(G, near, far, eps)
        A_ap = np.diag(np.full(dof, inv_chi)) - Ga
        s_ap, w_ap = smallest(A_ap)
        cond_ap = np.abs(w_ap).max() / np.abs(w_ap).min()
        print(f"  {eps:8.0e} {comp:7.3f} {s_ap[0]:12.4e} {np.abs(w_ap).max():12.4e} "
              f"{cond_ap:11.4e} {s_ap[0]/inv_chi:10.4f}")

    rad.UtiDelAll()


if __name__ == "__main__":
    main()
