#!/usr/bin/env python
"""Null-space dimension of the MSC geometric operator N vs block arrangement.

Follow-up to block_spectrum.py: the 4x4x2 block showed 33 exact-zero
eigenvalues of N.  Here we characterize how the null-space dimension scales
with the number / arrangement of cubes.  Uses SVD (robust rank for the
non-symmetric N).  null_dim = # singular values < tol * sigma_max.
"""
import numpy as np
import radia as rad

EDGE = 0.01
TOL = 1e-9  # relative singular-value threshold for "zero"


def build_block(nx, ny, nz, mu_r=1000.0, edge=EDGE):
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
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
    return rad.ObjCnt(objs), len(objs)


def null_dim(nx, ny, nz):
    block, n = build_block(nx, ny, nz)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(block))
    N = np.asarray(N, dtype=float)
    s = np.linalg.svd(N, compute_uv=False)
    nd = int(np.sum(s < TOL * s[0]))
    return n, dof, nd, s


print(f"{'block':>10} {'n_elem':>7} {'dof':>5} {'null_dim':>9} "
      f"{'null/elem':>10} {'sigma_min/max':>14}")
sizes = [(1, 1, 1), (2, 1, 1), (2, 2, 1), (2, 2, 2),
         (3, 3, 3), (4, 4, 2), (4, 4, 4), (2, 2, 8), (5, 5, 5)]
for (nx, ny, nz) in sizes:
    n, dof, nd, s = null_dim(nx, ny, nz)
    ratio = nd / n if n else 0.0
    print(f"{f'{nx}x{ny}x{nz}':>10} {n:7d} {dof:5d} {nd:9d} "
          f"{ratio:10.3f} {s[-1] / s[0]:14.3e}")

rad.UtiDelAll()
