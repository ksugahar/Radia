#!/usr/bin/env python
"""Visualize Yano's 'beautiful -> ugly' transition in magnetization space.

Phenomenon (Yano, 2026-05): BiCGSTAB from a zero initial guess passes a
'beautiful' (physically aligned) solution at intermediate iterations, then
twists into an 'ugly' (loop/vortex) solution = the LU / fully-converged
solution.  Tied to the near-null subspace of the MSC operator: at high mu_r the
system matrix A = diag(1/(mu_r-1)) - N has its smallest eigenvalues at ~1/(mu_r-1),
so a tiny RHS projection onto the near-null (loop) modes is amplified ~mu_r and
dominates the converged solution.  BiCGSTAB resolves large-eigenvalue (physical)
components first and the near-null (loop) components last.

This uses ONLY exposed Radia APIs (Solve + ObjM + SolverConfig): the BiCGSTAB
tolerance is the iteration knob (loose tol = few iterations = beautiful; tight
tol = converged = ugly).  No matrix extraction / face-normal mapping needed --
the magnetization is read directly per element via ObjM.

Usage: python beautiful_ugly_viz.py [nx ny nz] [mu_r]   default 6 6 1 1e5
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import radia as rad

EDGE = 0.01
D_HAT = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)   # diagonal applied field
B0 = 0.1


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
    cube = rad.ObjCnt(objs)
    ext = rad.ObjBckg(lambda p: [B0 * D_HAT[0], B0 * D_HAT[1], B0 * D_HAT[2]])
    return rad.ObjCnt([cube, ext]), cube


def solve_get(nx, ny, nz, mu_r, tol, method=1):
    grp, cube = build(nx, ny, nz, mu_r)
    rad.SolverConfig(bicgstab_tol=tol, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    rad.Solve(grp, 0.001, 2000, method)
    data = rad.ObjM(cube)
    C = np.array([d[0] for d in data], dtype=float)
    M = np.array([d[1] for d in data], dtype=float)
    return C, M


def alignment(M):
    nrm = np.linalg.norm(M, axis=1)
    g = nrm > 1e-12
    return float(np.mean((M[g] @ D_HAT) / nrm[g]))


def main():
    nx, ny, nz = 6, 6, 1
    mu_r = 1e5
    a = sys.argv[1:]
    if len(a) >= 3:
        nx, ny, nz = int(a[0]), int(a[1]), int(a[2])
    if len(a) >= 4:
        mu_r = float(a[3])

    tols = np.logspace(-1, -10, 19)
    aligns = []
    for t in tols:
        _, M = solve_get(nx, ny, nz, mu_r, t)
        aligns.append(alignment(M))
    aligns = np.array(aligns)

    C_b, M_b = solve_get(nx, ny, nz, mu_r, 1e-1)    # beautiful (few iter)
    C_u, M_u = solve_get(nx, ny, nz, mu_r, 1e-10)   # ugly (converged)
    print(f"block {nx}x{ny}x{nz} mu_r={mu_r:.0e}: beauty align={alignment(M_b):.3f}, "
          f"ugly align={alignment(M_u):.3f}")

    # pick first z-layer for quiver
    zmin = C_b[:, 2].min()
    sel = np.abs(C_b[:, 2] - zmin) < 0.5 * EDGE
    Xb, Yb = C_b[sel, 0] * 1e3, C_b[sel, 1] * 1e3

    def inplane(M):
        m = M[sel][:, :2]
        s = np.linalg.norm(M[sel], axis=1).mean()
        return m / (s + 1e-30)

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    ax[0].semilogx(tols, aligns, "o-")
    ax[0].invert_xaxis()
    ax[0].set_xlabel("BiCGSTAB tolerance (loose = few iter -> tight = converged)")
    ax[0].set_ylabel(r"mean alignment $\langle \hat M\cdot\hat d\rangle$")
    ax[0].set_title("(a) beautiful -> ugly")
    ax[0].axhline(1.0, color="g", ls=":", lw=1, label="perfect align")
    ax[0].grid(True, alpha=0.3); ax[0].legend(fontsize=8)

    mb = inplane(M_b)
    ax[1].quiver(Xb, Yb, mb[:, 0], mb[:, 1], np.hypot(mb[:, 0], mb[:, 1]),
                 cmap="viridis", scale=12, width=0.006)
    ax[1].set_title(f"(b) beautiful (tol=1e-1)\nalign={alignment(M_b):.3f}")
    ax[1].set_aspect("equal"); ax[1].set_xlabel("x [mm]"); ax[1].set_ylabel("y [mm]")

    mu = inplane(M_u)
    ax[2].quiver(Xb, Yb, mu[:, 0], mu[:, 1], np.hypot(mu[:, 0], mu[:, 1]),
                 cmap="viridis", scale=12, width=0.006)
    ax[2].set_title(f"(c) ugly / converged (tol=1e-10)\nalign={alignment(M_u):.3f}")
    ax[2].set_aspect("equal"); ax[2].set_xlabel("x [mm]"); ax[2].set_ylabel("y [mm]")

    fig.suptitle(f"MSC near-null loop contamination: {nx}x{ny}x{nz} block, "
                 f"mu_r={mu_r:.0e}, diagonal field (in-plane M, z={zmin*1e3:.0f}mm layer)")
    fig.tight_layout()
    import os
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       f"beautiful_ugly_{nx}x{ny}x{nz}_mu{int(mu_r)}.png")
    fig.savefig(out, dpi=140)
    print("saved:", out)
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
