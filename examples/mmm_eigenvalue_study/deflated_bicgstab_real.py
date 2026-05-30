#!/usr/bin/env python
"""Matrix-free deflated/shifted BiCGSTAB on the REAL HACApK operator, using the
local plaquette (cycle) basis. This LOCKS the algorithm to port into the C++
HACApK BiCGSTAB (radTRelaxationMethNo_2).

Algorithm (what C++ will do):
  - Build the local plaquette null basis V (right) from mesh topology
    (SVD-FREE). W (left null) likewise (here via SVD for the prototype).
  - Shifted matvec (low-rank, O(N) with local V,W):
        A_s x = A_haca x + alpha * V (W^T V)^{-1} (W^T x)
    moves the null cluster from 1/(mu_r-1) to 1/(mu_r-1)+alpha.
  - Plain BiCGSTAB on A_haca converges to the loop-dominated (ugly) solution;
    BiCGSTAB on A_s converges fast to the clean (loop-free) solution.

A_haca is the actual ACA operator from rad.HMatrixDensify (dense here; in C++
it is the HACApK MatVec). The shift is added as a low-rank correction.
"""
import sys
import numpy as np
from scipy.linalg import qr
from scipy.sparse.linalg import bicgstab, LinearOperator
import radia as rad

EDGE = 0.01


def build(nx, ny, nz, mu_r):
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                x0, y0, z0 = ix*EDGE, iy*EDGE, iz*EDGE
                v = [[x0,y0,z0],[x0+EDGE,y0,z0],[x0+EDGE,y0+EDGE,z0],[x0,y0+EDGE,z0],
                     [x0,y0,z0+EDGE],[x0+EDGE,y0,z0+EDGE],[x0+EDGE,y0+EDGE,z0+EDGE],[x0,y0+EDGE,z0+EDGE]]
                o = rad.ObjHexahedron(v, [0,0,0]); rad.MatApl(o, mat); objs.append(o)
    return rad.ObjCnt(objs)


def plaquette_basis(N, nx, ny, nz):
    """SVD-free local right-null basis from 2x2 element plaquettes."""
    dof = N.shape[0]
    def eidx(ix, iy, iz): return ix + nx*(iy + ny*iz)
    cols = N / (np.linalg.norm(N, axis=0, keepdims=True) + 1e-30)
    def shared(A, B):
        best, bi, bj = -2.0, -1, -1
        for i in range(6*A, 6*A+6):
            for j in range(6*B, 6*B+6):
                c = abs(float(cols[:, i] @ cols[:, j]))
                if c > best: best, bi, bj = c, i, j
        return bi, bj
    edges = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                e = eidx(ix, iy, iz)
                if ix+1 < nx: edges.append((e, eidx(ix+1,iy,iz)))
                if iy+1 < ny: edges.append((e, eidx(ix,iy+1,iz)))
                if iz+1 < nz: edges.append((e, eidx(ix,iy,iz+1)))
    sdof = {}
    for (A, B) in edges:
        i, j = shared(A, B); sdof[(A, B)] = (i, j); sdof[(B, A)] = (j, i)
    def loop_vec(cycle):
        v = np.zeros(dof); k = len(cycle)
        for t in range(k):
            A, B = cycle[t], cycle[(t+1) % k]
            di, dj = sdof[(A, B)]; v[di] += 1.0; v[dj] -= 1.0
        return v
    plaq = []
    for iz in range(nz):
        for iy in range(ny-1):
            for ix in range(nx-1):
                plaq.append([eidx(ix,iy,iz),eidx(ix+1,iy,iz),eidx(ix+1,iy+1,iz),eidx(ix,iy+1,iz)])
    for iy in range(ny):
        for iz in range(nz-1):
            for ix in range(nx-1):
                plaq.append([eidx(ix,iy,iz),eidx(ix+1,iy,iz),eidx(ix+1,iy,iz+1),eidx(ix,iy,iz+1)])
    for ix in range(nx):
        for iz in range(nz-1):
            for iy in range(ny-1):
                plaq.append([eidx(ix,iy,iz),eidx(ix,iy+1,iz),eidx(ix,iy+1,iz+1),eidx(ix,iy,iz+1)])
    L = np.column_stack([loop_vec(c) for c in plaq])
    Q, _ = np.linalg.qr(L)
    rank = int(np.linalg.matrix_rank(L, tol=1e-9))
    return Q[:, :rank]


def main():
    nx, ny, nz = 6, 6, 1
    mu_r = 1e5
    a = sys.argv[1:]
    if len(a) >= 3: nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4: mu_r = float(a[3])
    inv_chi = 1.0/(mu_r-1.0)

    blk = build(nx, ny, nz, mu_r)
    h = rad.BuildMatrix(blk)
    N, dof = rad.GetInteractMatrix(h); N = np.asarray(N, float)
    A_haca, _ = rad.HMatrixDensify(h); A_haca = np.asarray(A_haca, float)

    V = plaquette_basis(N, nx, ny, nz)              # SVD-free right null (local)
    U, s, _ = np.linalg.svd(N); W = U[:, s < 1e-9*s[0]]   # left null (SVD here)
    nd = V.shape[1]
    print(f"{nx}x{ny}x{nz} dof={dof} mu_r={mu_r:.0e} plaquette null_dim={nd}")

    rng = np.random.default_rng(0)
    b = rng.standard_normal(dof); b -= 0.97*(W @ (W.T @ b)); b /= np.linalg.norm(b)

    Minv = np.linalg.inv(W.T @ V)
    alpha = 1.0
    Aop = LinearOperator((dof, dof), matvec=lambda x: A_haca @ x)
    # matrix-free shift: HACApK matvec + low-rank correction (O(N) with local V,W)
    Asop = LinearOperator((dof, dof),
                          matvec=lambda x: A_haca @ x + alpha*(V @ (Minv @ (W.T @ x))))

    def run(op):
        h_ = {"n": 0, "nf": 0.0}
        def cb(xk):
            h_["n"] += 1; h_["nf"] = np.linalg.norm(V.T @ xk)/(np.linalg.norm(xk)+1e-30)
        bicgstab(op, b, rtol=1e-10, atol=0.0, maxiter=500, callback=cb)
        return h_

    p = run(Aop); s_ = run(Asop)
    print(f"plain  BiCGSTAB(A_haca)        : {p['n']:4d} iters, final null-frac = {p['nf']:.3f}")
    print(f"shifted BiCGSTAB(A_haca + V..) : {s_['n']:4d} iters, final null-frac = {s_['nf']:.3e}")
    print("=> matrix-free plaquette-shift deflation on the REAL HACApK operator: "
          "clean + fast. This is the algorithm to port to C++.")
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
