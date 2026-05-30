#!/usr/bin/env python
"""Why the null (loop) modes appear in the LATTER HALF of BiCGSTAB, and how a
shifted / deflated BiCGSTAB avoids them.

Theory (Krylov):
  x_k lives in the Krylov space; the error is e_k = A^{-1} p_k(A) b with a
  residual polynomial p_k, p_k(0)=1. The component of the solution along an
  eigenvector v_i (eigenvalue lam_i) reaches its exact value (b_i/lam_i) as
  p_k(lam_i) -> 0.
  - Large lam_i (physical modes): p_k(lam_i) drops in a FEW iterations
    -> physical part converges first ("beautiful").
  - Near-zero lam_i = 1/(mu_r-1) (null/loop modes): p_k(lam_i) ~ p_k(0) = 1
    for low degree k, so the (HUGE, ~mu_r-amplified) coefficient b_i/lam_i
    only enters once the polynomial finally resolves the tiny eigenvalues
    -> the loop content appears LATE ("ugly").
  Also: the residual component of a null mode is lam_i * (solution error) ~
  (1/mu_r)*error, so the RESIDUAL is nearly blind to the loop error -- a
  relative-residual stop cannot tell beautiful from ugly (Ida's point).

This script measures null-fraction ||V^T x_k|| / ||x_k|| per BiCGSTAB
iteration for plain A vs an eigenvalue-shifted A_s, and overlays the residual.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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


def main():
    nx, ny, nz = 4, 4, 2
    mu_r = 1e5
    a = sys.argv[1:]
    if len(a) >= 3: nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4: mu_r = float(a[3])
    inv_chi = 1.0/(mu_r-1.0)

    blk = build(nx, ny, nz, mu_r)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(blk))
    N = np.asarray(N, float)
    A = np.diag(np.full(dof, inv_chi)) - N
    U, s, Vt = np.linalg.svd(N)
    null = s < 1e-9*s[0]
    V = Vt[null].T; W = U[:, null]
    print(f"{nx}x{ny}x{nz} dof={dof} mu_r={mu_r:.0e} null_dim={int(null.sum())} cond(A)={np.linalg.cond(A):.2e}")

    rng = np.random.default_rng(0)
    b = rng.standard_normal(dof); b /= np.linalg.norm(b)
    # physical-like RHS: a uniform field weakly excites the loop (left-null)
    # subspace, so reduce b's left-null projection to a small value -> the
    # amplified loop content then builds up over the LATTER iterations.
    b = b - 0.97 * (W @ (W.T @ b))
    b /= np.linalg.norm(b)
    bn = np.linalg.norm(b)

    def run(M):
        hist = {"nf": [], "res": []}
        def cb(xk):
            nx_ = np.linalg.norm(xk) + 1e-30
            hist["nf"].append(np.linalg.norm(V.T @ xk)/nx_)
            hist["res"].append(np.linalg.norm(b - A @ xk)/bn)
        bicgstab(LinearOperator((dof,dof), matvec=lambda x: M @ x), b,
                 rtol=1e-12, atol=0.0, maxiter=400, callback=cb)
        return hist

    plain = run(A)
    alpha = 1.0
    A_s = A + alpha * (V @ np.linalg.inv(W.T @ V) @ W.T)
    shifted = run(A_s)
    print(f"plain  : {len(plain['nf'])} iters, final null-frac={plain['nf'][-1]:.3f}")
    print(f"shifted: {len(shifted['nf'])} iters, final null-frac={shifted['nf'][-1]:.3e}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    kp = np.arange(1, len(plain["nf"])+1)
    ks = np.arange(1, len(shifted["nf"])+1)
    ax[0].plot(kp, plain["nf"], "o-", ms=3, label="plain BiCGSTAB")
    ax[0].plot(ks, shifted["nf"], "s-", ms=3, label="eigenvalue-shifted")
    ax[0].set_xlabel("BiCGSTAB iteration"); ax[0].set_ylabel("null (loop) fraction  ||V^T x|| / ||x||")
    ax[0].set_title("(a) loop content vs iteration"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)

    ax2 = ax[1]; ax3 = ax2.twinx()
    l1, = ax2.semilogy(kp, plain["res"], "b-", label="rel. residual (plain)")
    l2, = ax3.plot(kp, plain["nf"], "r-", label="null fraction (plain)")
    ax2.set_xlabel("BiCGSTAB iteration"); ax2.set_ylabel("rel. residual", color="b")
    ax3.set_ylabel("null (loop) fraction", color="r")
    ax2.set_title("(b) residual is blind to the loop error")
    ax2.legend(handles=[l1, l2], fontsize=8, loc="center right")
    fig.suptitle(f"Null-mode dynamics in BiCGSTAB: {nx}x{ny}x{nz}, mu_r={mu_r:.0e}")
    fig.tight_layout()
    import os
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       f"bicgstab_dynamics_{nx}x{ny}x{nz}_mu{int(mu_r)}.png")
    fig.savefig(out, dpi=140); print("saved:", out)
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
