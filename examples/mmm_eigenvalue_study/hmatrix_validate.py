#!/usr/bin/env python
"""(B) Validate the ACTUAL HACApK (ACA+) operator vs the exact dense operator.

Uses the new rad.HMatrixDensify(handle) binding (densifies the real HACApK
system matrix A = -N + diag(1/chi) via MatVec on unit vectors) and compares it
with the exact A = diag(1/chi) - N from rad.GetInteractMatrix.

Answers the original "Phase B" question (does the real ACA+ shift the near-zero
eigenvalues?) and cross-checks the A = diag(1/chi) - N sign convention used in
the eigenvalue study.

Usage: python hmatrix_validate.py [n] [mu_r] [eps]   default 6 1e5 1e-4
"""
import sys
import numpy as np
import radia as rad

EDGE = 0.01


def build(n, mu_r):
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs = []
    for iz in range(n):
        for iy in range(n):
            for ix in range(n):
                x0, y0, z0 = ix*EDGE, iy*EDGE, iz*EDGE
                v = [[x0,y0,z0],[x0+EDGE,y0,z0],[x0+EDGE,y0+EDGE,z0],[x0,y0+EDGE,z0],
                     [x0,y0,z0+EDGE],[x0+EDGE,y0,z0+EDGE],[x0+EDGE,y0+EDGE,z0+EDGE],[x0,y0+EDGE,z0+EDGE]]
                o = rad.ObjHexahedron(v, [0,0,0]); rad.MatApl(o, mat); objs.append(o)
    return rad.ObjCnt(objs)


def smallest(M, k=6):
    w = np.linalg.eigvals(M)
    aw = np.sort(np.abs(w))
    return aw[:k], aw[-1]/aw[0]


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    mu_r = float(sys.argv[2]) if len(sys.argv) > 2 else 1e5
    eps = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-4
    chi = mu_r - 1.0
    inv_chi = 1.0/chi

    rad.SolverConfig(hacapk_eps=eps, hacapk_leaf=10, hacapk_eta=2.0)
    blk = build(n, mu_r)
    h = rad.BuildMatrix(blk)

    N, dof = rad.GetInteractMatrix(h)
    N = np.asarray(N, float)
    A_exact = np.diag(np.full(dof, inv_chi)) - N

    A_haca, dof2 = rad.HMatrixDensify(h)
    A_haca = np.asarray(A_haca, float)
    print(f"n={n}^3 dof={dof} (haca dof={dof2}) mu_r={mu_r:.0e} eps={eps:.0e}")
    print(f"A_haca shape={A_haca.shape}")

    # convention + ACA error: A_haca should equal diag(1/chi)-N up to ACA error
    rel = np.linalg.norm(A_haca - A_exact) / np.linalg.norm(A_exact)
    diag_ok = np.allclose(np.diag(A_haca), np.diag(A_exact), rtol=1e-6, atol=1e-8)
    print(f"\n||A_haca - A_exact||/||A_exact|| = {rel:.3e}   (= ACA far-block error)")
    print(f"diagonals match (sign convention OK): {diag_ok}")

    se, ce = smallest(A_exact)
    sh, ch = smallest(A_haca)
    print(f"\nexact  : smallest|lam|={np.array2string(se, precision=3)}  cond={ce:.3e}")
    print(f"hacapk : smallest|lam|={np.array2string(sh, precision=3)}  cond={ch:.3e}")
    print(f"1/(mu_r-1) = {inv_chi:.3e}")
    print(f"\n=> does real ACA+ shift the near-zero eigenvalues away from 1/(mu_r-1)? "
          f"min|lam| haca/exact = {sh[0]/se[0]:.4f}")

    rad.UtiDelAll()


if __name__ == "__main__":
    main()
