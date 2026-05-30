#!/usr/bin/env python
"""Eigenvalue spectrum of the MMM/MSC interaction operator (BLOCK study).

Background (CEFC 2026 full paper, eigenvalue study):
  Yano observed that BiCGSTAB from a zero initial guess passes a 'beautiful'
  (all-aligned, physical) solution at intermediate iterations and then twists
  into an 'ugly' (loop/vortex) solution = the LU solution.  Hypothesis
  (Sugahara / Ida): the exact solution carries a near-zero-eigenvalue mode of
  the MSC operator -- a magnetization-loop pattern that produces almost no
  external field.  This script tests step 1: do such near-zero eigenvalues
  exist in the geometric interaction operator?

What this computes:
  N   = rad.GetInteractMatrix(rad.BuildMatrix(block))  -- mu_r-independent
        geometric interaction operator (CLAUDE.md: stored +N, non-symmetric
        MSC collocation).  dof = 6 * n_elem.
  A   = diag(1/(mu_r-1)) - N                            -- system matrix that
        LU / BiCGSTAB / HACApK actually solve (sign per CLAUDE.md ComputeEntry;
        verified separately by solver calibration).
  eig(N), eig(A) via numpy.linalg.eig (general, complex spectra).

Usage:
  python block_spectrum.py [nx ny nz] [mu_r]
  default: 4 4 2 (=32 cubes), mu_r=5000  (Yano's BLOCK48 description)
"""
import sys
import os
import numpy as np
import radia as rad

EDGE = 0.01  # cube edge length [m]


def build_block(nx, ny, nz, mu_r, edge=EDGE):
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


def main():
    nx, ny, nz = 4, 4, 2
    mu_r = 5000.0
    args = sys.argv[1:]
    if len(args) >= 3:
        nx, ny, nz = int(args[0]), int(args[1]), int(args[2])
    if len(args) >= 4:
        mu_r = float(args[3])

    n_elem = nx * ny * nz
    print(f"BLOCK {nx}x{ny}x{nz} = {n_elem} cubes, mu_r={mu_r}, dof={6*n_elem}")

    block, n = build_block(nx, ny, nz, mu_r)
    h = rad.BuildMatrix(block)
    N, dof = rad.GetInteractMatrix(h)
    N = np.asarray(N, dtype=float)
    print(f"N shape={N.shape}  symmetric={np.allclose(N, N.T, atol=1e-9)}"
          f"  maxasym={np.max(np.abs(N - N.T)):.3e}")

    inv_chi = 1.0 / (mu_r - 1.0)
    A = np.diag(np.full(dof, inv_chi)) - N   # system matrix (sign per CLAUDE.md)

    for name, Mtx in (("N (geometric)", N), ("A = diag(1/(mu_r-1)) - N", A)):
        w = np.linalg.eigvals(Mtx)
        idx = np.argsort(np.abs(w))
        w = w[idx]
        print(f"\n=== eig({name}) ===")
        print(f"  |lambda| range: [{np.abs(w[0]):.3e}, {np.abs(w[-1]):.3e}]")
        print(f"  smallest 12 |lambda|:")
        for k in range(min(12, dof)):
            print(f"    {k:3d}: {w[k].real:+.6e} {w[k].imag:+.6e}j  |.|={np.abs(w[k]):.6e}")
        thr = 10.0 * inv_chi
        n_near = int(np.sum(np.abs(w) < thr))
        print(f"  count |lambda| < 10/(mu_r-1) (={thr:.2e}): {n_near}")
        print(f"  count |lambda| < 1e-3: {int(np.sum(np.abs(w) < 1e-3))}")

    # Save full spectra + eigenvectors of N for later mode visualization
    wN, vN = np.linalg.eig(N)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       f"spectrum_{nx}x{ny}x{nz}_mu{int(mu_r)}.npz")
    np.savez(out, N=N, A=A, eigval_N=wN, eigvec_N=vN,
             mu_r=mu_r, dims=np.array([nx, ny, nz]), inv_chi=inv_chi)
    print(f"\nsaved: {out}")
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
