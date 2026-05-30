#!/usr/bin/env python
"""Remove the MSC near-null (loop) modes so the converged solution is clean.

Two textbook techniques, prototyped on the dense system (small model) to
validate the concept before any HACApK/C++ integration:

  (1) right-projection DEFLATION of the solution
        sigma_clean = sigma - V (V^T sigma)
      removes the loop content from a given solution. Because N V = 0, the
      collocation field N*sigma is preserved EXACTLY.

  (2) spectral EIGENVALUE SHIFT (Ida's suggestion, "固有値シフト付き")
        A_s = A + alpha * V (W^T V)^{-1} W^T
      moves the cluster of near-zero eigenvalues (the null modes of N, which
      sit at exactly 1/(mu_r-1) in A) up to 1/(mu_r-1)+alpha, so a tight
      converged solve no longer amplifies them -> the converged solution is
      loop-free without relying on lucky early stopping.

V = right null space of N (SVD right-singular vecs with sigma~0)
W = left  null space of N (SVD left-singular  vecs with sigma~0)

Honest scalability note printed at the end: dim(null) ~ 20% of DOF, so
forming V explicitly is O(N * 0.2N) = O(N^2) and defeats HACApK. The modes
are LOCAL loops (minimal mode = 2x2 plaquette), so the scalable path is a
local cycle/loop basis (loop-star style), not dense SVD.

Usage: python nullmode_removal.py [nx ny nz] [mu_r]   default 4 4 2 1e5
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
                x0, y0, z0 = ix * EDGE, iy * EDGE, iz * EDGE
                v = [[x0, y0, z0], [x0 + EDGE, y0, z0], [x0 + EDGE, y0 + EDGE, z0], [x0, y0 + EDGE, z0],
                     [x0, y0, z0 + EDGE], [x0 + EDGE, y0, z0 + EDGE], [x0 + EDGE, y0 + EDGE, z0 + EDGE], [x0, y0 + EDGE, z0 + EDGE]]
                o = rad.ObjHexahedron(v, [0, 0, 0])
                rad.MatApl(o, mat)
                objs.append(o)
    return rad.ObjCnt(objs)


def main():
    nx, ny, nz = 4, 4, 2
    mu_r = 1e5
    a = sys.argv[1:]
    if len(a) >= 3:
        nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4:
        mu_r = float(a[3])
    chi = mu_r - 1.0
    inv_chi = 1.0 / chi

    block = build(nx, ny, nz, mu_r)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(block))
    N = np.asarray(N, dtype=float)
    A = np.diag(np.full(dof, inv_chi)) - N
    print(f"block {nx}x{ny}x{nz}={nx*ny*nz} elem, dof={dof}, mu_r={mu_r:.0e}, 1/chi={inv_chi:.3e}")

    # --- null space of N via SVD ---
    U, s, Vt = np.linalg.svd(N)
    tol = 1e-9 * s[0]
    null = s < tol
    nd = int(null.sum())
    V = Vt[null].T.copy()          # right null space  (dof x nd)
    W = U[:, null].copy()          # left  null space  (dof x nd)
    print(f"null-space dim = {nd}  ({100.0*nd/dof:.1f}% of DOF)")

    condA = np.linalg.cond(A)
    print(f"cond(A) = {condA:.3e}   (min|eig|~1/chi={inv_chi:.2e})")

    # --- physical-ish RHS (deterministic); mechanism is robust to b ---
    rng = np.random.default_rng(0)
    b = rng.standard_normal(dof)
    b /= np.linalg.norm(b)

    sigma = np.linalg.solve(A, b)            # "ugly" converged solution
    def null_frac(x):
        return np.linalg.norm(V.T @ x) / (np.linalg.norm(x) + 1e-30)
    print(f"\nconverged solve:   ||sigma||={np.linalg.norm(sigma):.3e}  "
          f"null-fraction={null_frac(sigma):.3f}")

    # --- (1) projection deflation ---
    sigma_def = sigma - V @ (V.T @ sigma)
    field_change = np.linalg.norm(N @ sigma_def - N @ sigma) / (np.linalg.norm(N @ sigma) + 1e-30)
    print(f"(1) projection:    null-fraction={null_frac(sigma_def):.3e}  "
          f"rel field change ||N(d-sigma)||={field_change:.2e}  (N V=0 -> ~0)")

    # --- (2) spectral eigenvalue shift A_s = A + alpha V (W^T V)^-1 W^T ---
    WtV = W.T @ V
    print(f"\n(2) eigenvalue shift  (cluster currently at 1/chi={inv_chi:.2e}):")
    print(f"     {'alpha':>10} {'cond(A_s)':>12} {'null-frac(sigma_s)':>18} {'||sig_s-sig_def||/||sig_def||':>30}")
    Minv = np.linalg.inv(WtV)
    shift_op = V @ Minv @ W.T
    for alpha in (1e-3, 1e-2, 1e-1, 1.0):
        A_s = A + alpha * shift_op
        sigma_s = np.linalg.solve(A_s, b)
        rel = np.linalg.norm(sigma_s - sigma_def) / (np.linalg.norm(sigma_def) + 1e-30)
        print(f"     {alpha:10.0e} {np.linalg.cond(A_s):12.3e} {null_frac(sigma_s):18.3e} {rel:30.3e}")

    print("\nSCALABILITY: full null basis V is dof x nd with nd ~ "
          f"{100.0*nd/dof:.0f}% of DOF -> O(N^2) storage, defeats HACApK.")
    print("The modes are LOCAL loops (minimal = 2x2 plaquette) -> a local")
    print("cycle/loop basis (loop-star style) is the scalable path for the")
    print("HACApK BiCGSTAB solver (build the loop basis from mesh topology,")
    print("deflate matrix-free). Dense SVD here is a concept prototype only.")
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
