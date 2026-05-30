#!/usr/bin/env python
"""Can the C++ deflation use ONLY the plaquette right-null basis V (no SVD, no
separate left-null W)?  Two questions:

  (1) Is the plaquette basis V ALSO a LEFT null basis (||N^T V|| ~ 0)?
      If yes, W = V and the oblique shift V(W^T V)^-1 W^T = V(V^T V)^-1 V^T
      needs only V.
  (2) Does the V-only SYMMETRIC shift A_s = A + alpha V V^T recondition and
      give the clean (loop-free) solution despite A being non-symmetric?

Reference clean solution: sigma_clean = sigma_ugly - V(V^T V)^-1 V^T sigma_ugly
(projection removes the loops; field N*sigma preserved since N V = 0).
"""
import sys
import numpy as np
import radia as rad
from deflated_bicgstab_real import build, plaquette_basis  # reuse


def main():
    nx, ny, nz = 4, 4, 2
    mu_r = 1e5
    a = sys.argv[1:]
    if len(a) >= 3: nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4: mu_r = float(a[3])
    inv_chi = 1.0/(mu_r-1.0)

    blk = build(nx, ny, nz, mu_r)
    h = rad.BuildMatrix(blk)
    N, dof = rad.GetInteractMatrix(h); N = np.asarray(N, float)
    A = np.diag(np.full(dof, inv_chi)) - N
    V = plaquette_basis(N, nx, ny, nz)      # orthonormal plaquette right-null
    nd = V.shape[1]
    print(f"{nx}x{ny}x{nz} dof={dof} mu_r={mu_r:.0e} plaquette null_dim={nd}")

    # (1) is V also a LEFT null basis?
    rnull = np.linalg.norm(N @ V) / np.linalg.norm(V)
    lnull = np.linalg.norm(N.T @ V) / np.linalg.norm(V)
    print(f"\n(1) ||N V||/||V|| = {rnull:.2e} (right null)   "
          f"||N^T V||/||V|| = {lnull:.2e} (left null?)")
    print(f"    => V is {'ALSO a left-null basis' if lnull < 1e-8 else 'NOT a left-null basis'}")

    # reference clean (projection) solution
    rng = np.random.default_rng(0)
    b = rng.standard_normal(dof); b /= np.linalg.norm(b)
    Pp = V @ V.T
    sig_ugly = np.linalg.solve(A, b)
    sig_ref = sig_ugly - Pp @ sig_ugly
    nf = lambda x: np.linalg.norm(V.T @ x)/(np.linalg.norm(x)+1e-30)

    # (2) V-only symmetric shift
    alpha = 1.0
    A_sym = A + alpha * Pp
    sig_sym = np.linalg.solve(A_sym, b)
    cond_sym = np.linalg.cond(A_sym)
    agree_sym = np.linalg.norm(sig_sym - sig_ref) / (np.linalg.norm(sig_ref)+1e-30)
    print(f"\n(2) V-only symmetric shift A + {alpha}*V V^T:")
    print(f"    cond {np.linalg.cond(A):.2e} -> {cond_sym:.2e}, null-frac={nf(sig_sym):.3e}, "
          f"rel diff vs projection-clean = {agree_sym:.3e}")

    # V+W oblique shift (W from SVD) for comparison
    U, s, _ = np.linalg.svd(N); W = U[:, s < 1e-9*s[0]]
    A_obl = A + alpha * (V @ np.linalg.inv(W.T @ V) @ W.T)
    sig_obl = np.linalg.solve(A_obl, b)
    print(f"    [ref] V+W oblique shift: cond {np.linalg.cond(A_obl):.2e}, "
          f"null-frac={nf(sig_obl):.3e}, rel diff vs clean={np.linalg.norm(sig_obl-sig_ref)/np.linalg.norm(sig_ref):.3e}")

    print("\nConclusion: C++ needs only the plaquette V if (1) lnull~0 OR (2) the "
          "V-only shift matches the projection-clean solution.")
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
