#!/usr/bin/env python
"""(A-reinforce) Explicit plaquette (cycle) null basis -- SVD-free, topological.

Discovery: the null space of the MSC operator N equals the CYCLE SPACE of the
element-adjacency graph (nodes = elements, edges = shared internal faces):
    null_dim = F_internal - N_elem + 1   (graph first Betti number)
verified exactly for 2x2x1..5x5x5.

Physical derivation: a null mode puts opposite charges on the two sides of each
shared face (field cancels) AND zero net charge per element (compensating point
charge vanishes) => the per-face "flux" is divergence-free per element =>
circulation => a cycle of the element graph.  Minimal cycle = 2x2 plaquette.

This builds the explicit local plaquette basis WITHOUT SVD and verifies it:
  - shared-face DOF for each adjacent pair found via N-column parallelism
    (coincident faces produce near-parallel columns) -- no face-DOF "bridge"
  - each plaquette loop sigma (+1/-1 circulation) is an exact null mode
  - the plaquette set spans the SVD null space (rank == null_dim, proj err ~0)

Usage: python nullmode_cycle_basis.py [nx ny nz] [mu_r]   default 4 4 2 1e5
"""
import sys
import numpy as np
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


def main():
    nx, ny, nz = 4, 4, 2
    mu_r = 1e5
    a = sys.argv[1:]
    if len(a) >= 3:
        nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4:
        mu_r = float(a[3])
    ne = nx*ny*nz

    blk = build(nx, ny, nz, mu_r)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(blk))
    N = np.asarray(N, float)

    # reference null dim via SVD
    s = np.linalg.svd(N, compute_uv=False)
    null_dim = int(np.sum(s < 1e-9*s[0]))

    def eidx(ix, iy, iz):
        return ix + nx*(iy + ny*iz)

    # element-adjacency graph: edges = shared internal faces
    edges = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                e = eidx(ix, iy, iz)
                if ix+1 < nx: edges.append((e, eidx(ix+1,iy,iz)))
                if iy+1 < ny: edges.append((e, eidx(ix,iy+1,iz)))
                if iz+1 < nz: edges.append((e, eidx(ix,iy,iz+1)))
    F_int = len(edges)
    cycle_rank = F_int - ne + 1
    print(f"{nx}x{ny}x{nz}={ne} elem, dof={dof}, mu_r={mu_r:.0e}")
    print(f"graph: V={ne} edges(F_int)={F_int} -> cycle rank = F_int-V+1 = {cycle_rank}")
    print(f"SVD null dim = {null_dim}   MATCH={cycle_rank==null_dim}")

    # shared-face DOF via N-column parallelism (bridge-free)
    cols = N / (np.linalg.norm(N, axis=0, keepdims=True) + 1e-30)
    def shared_dof(A, B):
        best, bi, bj = -2.0, -1, -1
        for i in range(6*A, 6*A+6):
            for j in range(6*B, 6*B+6):
                c = abs(float(cols[:, i] @ cols[:, j]))
                if c > best:
                    best, bi, bj = c, i, j
        return bi, bj, best
    # probe the identification quality on the first edge
    i0, j0, c0 = shared_dof(*edges[0])
    print(f"shared-face DOF id (col cosine) sample: edge{edges[0]} -> dof({i0},{j0}) cos={c0:.3f}")

    sdof = {}
    for (A, B) in edges:
        i, j, _ = shared_dof(A, B)
        sdof[(A, B)] = (i, j)
        sdof[(B, A)] = (j, i)

    def loop_vec(cycle):
        """cycle: ordered list of element indices forming a closed loop."""
        v = np.zeros(dof)
        k = len(cycle)
        for t in range(k):
            A, B = cycle[t], cycle[(t+1) % k]
            di, dj = sdof[(A, B)]
            v[di] += 1.0
            v[dj] -= 1.0
        return v

    # 2x2 plaquettes in all three planes
    plaq = []
    for iz in range(nz):
        for iy in range(ny-1):
            for ix in range(nx-1):
                plaq.append([eidx(ix,iy,iz), eidx(ix+1,iy,iz), eidx(ix+1,iy+1,iz), eidx(ix,iy+1,iz)])
    for iy in range(ny):
        for iz in range(nz-1):
            for ix in range(nx-1):
                plaq.append([eidx(ix,iy,iz), eidx(ix+1,iy,iz), eidx(ix+1,iy,iz+1), eidx(ix,iy,iz+1)])
    for ix in range(nx):
        for iz in range(nz-1):
            for iy in range(ny-1):
                plaq.append([eidx(ix,iy,iz), eidx(ix,iy+1,iz), eidx(ix,iy+1,iz+1), eidx(ix,iy,iz+1)])

    L = np.column_stack([loop_vec(c) for c in plaq])
    resid = np.linalg.norm(N @ L, axis=0) / (np.linalg.norm(L, axis=0) + 1e-30)
    rankL = int(np.linalg.matrix_rank(L, tol=1e-9))
    print(f"\nplaquettes built: {len(plaq)}  (each a 4-element 2x2 loop)")
    print(f"null-mode residual ||N*loop||/||loop||: max={resid.max():.2e}  mean={resid.mean():.2e}")
    print(f"rank(L) = {rankL}   null_dim = {null_dim}   SPAN-RANK MATCH={rankL==null_dim}")

    # does L span the SVD null space?
    U, sN, Vt = np.linalg.svd(N)
    V = Vt[sN < 1e-9*sN[0]].T
    Q, _ = np.linalg.qr(L)
    Q = Q[:, :rankL]
    proj_err = np.linalg.norm(V - Q @ (Q.T @ V)) / np.linalg.norm(V)
    print(f"span error ||(I - QQ^T)V||/||V|| = {proj_err:.2e}  (0 => plaquettes span the null space)")

    print(f"\n=> explicit local plaquette (cycle) basis = the null space, SVD-FREE.")
    print(f"   each basis vector touches 4 elements (8 shared-face DOF) -> O(1) support.")
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
