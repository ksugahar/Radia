#!/usr/bin/env python
"""C++-viable shift: A + alpha * L L^T with the RAW (non-orthonormal) plaquette
basis L -- NO orthonormalization (O(N^3)), NO (L^T L)^-1 inverse.

Each matvec adds the low-rank term L(L^T x): L^T x is O(nnz)=O(N) (8 entries
per plaquette), then L(...) is O(N). L L^T acts ONLY on span(L) = the null
space (L L^T x = 0 for x perpendicular to every plaquette), so it lifts the
null eigenvalues of A away from zero without touching the physical complement
(empirically; A is non-normal so we verify numerically).

This is the exact operation to put in RadHACApKMSCManager::MatVec.
"""
import sys
import numpy as np
import radia as rad
from deflated_bicgstab_real import build


def raw_plaquettes(N, nx, ny, nz):
    """RAW plaquette loop vectors L (dof x n_plaq), non-orthonormal, +/-1 entries."""
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
    def loop_vec(cyc):
        v = np.zeros(dof); k = len(cyc)
        for t in range(k):
            A, B = cyc[t], cyc[(t+1) % k]
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
    return np.column_stack([loop_vec(c) for c in plaq])


def main():
    nx, ny, nz = 4, 4, 2
    mu_r = 1e5
    a = sys.argv[1:]
    if len(a) >= 3: nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4: mu_r = float(a[3])
    inv_chi = 1.0/(mu_r-1.0)

    blk = build(nx, ny, nz, mu_r)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(blk)); N = np.asarray(N, float)
    A = np.diag(np.full(dof, inv_chi)) - N
    L = raw_plaquettes(N, nx, ny, nz)
    Vo, _ = np.linalg.qr(L); Vo = Vo[:, :np.linalg.matrix_rank(L, tol=1e-9)]   # for MEASURING only
    nf = lambda x: np.linalg.norm(Vo.T @ x)/(np.linalg.norm(x)+1e-30)
    print(f"{nx}x{ny}x{nz} dof={dof} mu_r={mu_r:.0e}  L shape={L.shape} (raw plaquettes), "
          f"||N L||/||L||={np.linalg.norm(N@L)/np.linalg.norm(L):.1e}")

    rng = np.random.default_rng(0)
    b = rng.standard_normal(dof); b /= np.linalg.norm(b)
    sig_ugly = np.linalg.solve(A, b)
    sig_ref = sig_ugly - Vo @ (Vo.T @ sig_ugly)   # projection-clean reference

    print(f"\nA + alpha*L L^T  (RAW L, no orthonormalization, no inverse):")
    print(f"  {'alpha':>8} {'cond(A_s)':>11} {'null-frac':>10} {'reldiff vs clean':>17}")
    for alpha in (1e-2, 1e-1, 1.0, 10.0, 100.0):
        A_s = A + alpha * (L @ L.T)
        sig = np.linalg.solve(A_s, b)
        rd = np.linalg.norm(sig - sig_ref)/(np.linalg.norm(sig_ref)+1e-30)
        print(f"  {alpha:8.0e} {np.linalg.cond(A_s):11.3e} {nf(sig):10.3e} {rd:17.3e}")
    print(f"  (cond(A) without shift = {np.linalg.cond(A):.3e})")
    print("\n=> if some alpha reconditions (cond ~O(10-100)) AND null-frac drops,")
    print("   the RAW-L shift is C++-viable: matvec += alpha*L(L^T x), O(N), no SVD.")
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
