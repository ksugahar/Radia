#!/usr/bin/env python
"""Loop-star gauge prototype (dense, small) -- verify the ALPHA-FREE cure BEFORE
the C++ HACApK sandwich.

Idea: the loop null space (ker N) is what makes A = diag(1/chi) - N ill-conditioned
(loop eigenvalues at 1/chi). Instead of the fragile alpha-shift, RESTRICT the solve
to the STAR (non-loop) subspace: solve  (U_s^T A U_s) y = U_s^T b,  sigma = U_s y,
where U_s spans the complement of ker N. No alpha, no cliff.

This dense prototype (SVD to get ker N) checks, on a small cube:
  LINEAR (uniform chi): restriction is FIELD-EXACT (N*sigma matches the full solve)
    and WELL-CONDITIONED (cond(U_s^T A U_s) << cond(A)) -- the clean win.
  NONLINEAR proxy (non-uniform chi): restriction DRIFTS (loops couple to star via
    the non-scalar diag(1/chi)) -- documents the caveat.

The C++ version will replace the dense U_s/SVD with the SPARSE loop basis + a
sandwich around the unchanged HACApK MatVec; this prototype only validates the math.

Usage: python loopstar_prototype.py [n]   (cube n^3, default 4)
"""
import sys
import numpy as np
import radia as rad

EDGE = 0.01


def build_cube(n):
    rad.UtiDelAll()
    mat = rad.MatLin(1000.0)   # value irrelevant; we rebuild A by hand below
    objs = []
    for iz in range(n):
        for iy in range(n):
            for ix in range(n):
                x0, y0, z0 = ix*EDGE, iy*EDGE, iz*EDGE
                v = [[x0,y0,z0],[x0+EDGE,y0,z0],[x0+EDGE,y0+EDGE,z0],[x0,y0+EDGE,z0],
                     [x0,y0,z0+EDGE],[x0+EDGE,y0,z0+EDGE],[x0+EDGE,y0+EDGE,z0+EDGE],[x0,y0+EDGE,z0+EDGE]]
                o = rad.ObjHexahedron(v, [0,0,0]); rad.MatApl(o, mat); objs.append(o)
    return rad.ObjCnt(objs)


def loop_star_split(N, tol=1e-9):
    """SVD of N: right-singular vectors with ~zero singular value span ker N (loops);
    the rest span the star (non-loop) subspace U_s."""
    U, s, Vt = np.linalg.svd(N)
    smax = s[0]
    loop_mask = s < tol * smax
    V = Vt.T
    L = V[:, loop_mask]            # loops (ker N)
    U_s = V[:, ~loop_mask]         # star (complement)
    return L, U_s, int(loop_mask.sum())


def run(N, chi, b, U_s):
    A = np.diag(1.0/chi) - N
    sig_full = np.linalg.solve(A, b)
    Ar = U_s.T @ A @ U_s
    y = np.linalg.solve(Ar, U_s.T @ b)
    sig_r = U_s @ y
    # field-relevant part = N*sigma (demag field at collocation; loops are killed by N)
    f_full = N @ sig_full
    f_r = N @ sig_r
    fielddev = np.linalg.norm(f_r - f_full) / (np.linalg.norm(f_full) + 1e-30)
    return np.linalg.cond(A), np.linalg.cond(Ar), fielddev


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    grp = build_cube(n)
    Nmat, dof = rad.GetInteractMatrix(rad.BuildMatrix(grp))
    Nmat = np.asarray(Nmat, float)
    L, U_s, ndim = loop_star_split(Nmat)
    rng = np.random.default_rng(0)
    b = rng.standard_normal(dof)
    print(f"cube {n}^3: dof={dof}, null(loops)={ndim}, star={U_s.shape[1]}")
    print(f"{'case':<22} {'cond(A)':>11} {'cond(A_star)':>13} {'field dev':>11}")

    # LINEAR: uniform chi (mu_r = 1e5 -> chi = 1e5-1)
    chi_u = np.full(dof, 1e5 - 1.0)
    cA, cAr, fd = run(Nmat, chi_u, b, U_s)
    print(f"{'linear mu=1e5':<22} {cA:>11.2e} {cAr:>13.2e} {fd:>11.2e}")

    # NONLINEAR proxy: non-uniform chi (random per-DOF, spanning saturation)
    chi_n = 10.0 ** rng.uniform(0.0, 5.0, dof)   # chi in [1, 1e5], strongly varying
    cA, cAr, fd = run(Nmat, chi_n, b, U_s)
    print(f"{'nonlin proxy (varied)':<22} {cA:>11.2e} {cAr:>13.2e} {fd:>11.2e}")

    rad.UtiDelAll()
    print("\nExpect: LINEAR -> field dev ~ machine eps (field-exact) AND")
    print("        cond(A_star) << cond(A) (loops removed the ill-conditioning, alpha-free).")
    print("        NONLINEAR proxy -> field dev sizable (loops couple to star: the caveat).")


if __name__ == '__main__':
    main()
