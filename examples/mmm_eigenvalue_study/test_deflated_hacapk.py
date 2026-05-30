#!/usr/bin/env python
"""End-to-end test of the C++ matrix-free loop-mode deflation in HACApK.

Builds the local plaquette (cycle) basis L (CSR), passes it via
rad.SetHACApKDeflation, and solves a high-mu_r block with method=2 (HACApK)
WITHOUT vs WITH deflation. Without deflation the converged solution is the
loop-dominated "ugly" one (low alignment); with deflation it is the clean
"beautiful" one (high alignment).

Usage: python test_deflated_hacapk.py [nx ny nz] [mu_r] [alpha]
"""
import sys
import numpy as np
import radia as rad

EDGE = 0.01
D = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
B0 = 0.1


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
    cube = rad.ObjCnt(objs)
    ext = rad.ObjBckg(lambda p: [B0*D[0], B0*D[1], B0*D[2]])
    return rad.ObjCnt([cube, ext]), cube


def plaquette_csr(N, nx, ny, nz):
    """Build CSR (offsets, dofs, signs) of the 2x2 plaquette loop basis."""
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
    offs, dofs, signs = [0], [], []
    for cyc in plaq:
        k = len(cyc)
        for t in range(k):
            A, B = cyc[t], cyc[(t+1) % k]
            di, dj = sdof[(A, B)]
            dofs += [di, dj]; signs += [1.0, -1.0]
        offs.append(len(dofs))
    return (np.asarray(offs, np.int32), np.asarray(dofs, np.int32),
            np.asarray(signs, np.float64), len(plaq))


def align(cube):
    M = np.array([m[1] for m in rad.ObjM(cube)], float)
    nrm = np.linalg.norm(M, axis=1); g = nrm > 1e-12
    return float(np.mean((M[g] @ D) / nrm[g]))


def main():
    nx, ny, nz = 4, 4, 2
    mu_r, alpha = 1e5, 1.0
    a = sys.argv[1:]
    if len(a) >= 3: nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4: mu_r = float(a[3])
    if len(a) >= 5: alpha = float(a[4])

    print("has SetHACApKDeflation:", hasattr(rad, "SetHACApKDeflation"))

    # plaquette basis from a probe build
    grp, cube = build(nx, ny, nz, mu_r)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(grp))
    offs, dofs, signs, nplaq = plaquette_csr(np.asarray(N, float), nx, ny, nz)
    print(f"{nx}x{ny}x{nz} dof={dof} mu_r={mu_r:.0e}: {nplaq} plaquettes, nnz={len(dofs)}")

    rad.SolverConfig(bicgstab_tol=1e-10, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

    # undeflated HACApK
    rad.SetHACApKDeflation(np.zeros(0, np.int32), np.zeros(0, np.int32), np.zeros(0), 0.0)
    grp, cube = build(nx, ny, nz, mu_r)
    rad.Solve(grp, 0.001, 2000, 2)
    a_undef = align(cube)

    # deflated HACApK
    rad.SetHACApKDeflation(offs, dofs, signs, alpha)
    grp, cube = build(nx, ny, nz, mu_r)
    rad.Solve(grp, 0.001, 2000, 2)
    a_def = align(cube)

    print(f"\nHACApK undeflated : mean alignment = {a_undef:.3f}  (loop-dominated = ugly)")
    print(f"HACApK deflated   : mean alignment = {a_def:.3f}  (clean = beautiful)")
    print("=> deflation works in C++" if a_def > a_undef + 0.2 else "=> CHECK: little change")
    rad.SetHACApKDeflation(np.zeros(0, np.int32), np.zeros(0, np.int32), np.zeros(0), 0.0)
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
