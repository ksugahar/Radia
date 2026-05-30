#!/usr/bin/env python
"""(A) Local null basis + deflated BiCGSTAB -- scalable removal concept proof.

Goal: remove the MSC near-null (loop) modes with a LOCAL basis (so the
deflation is O(N), compatible with HACApK), and show a deflated/shifted
BiCGSTAB converges fast to a loop-free solution.

Local basis WITHOUT the face-DOF bridge: the null projector P = V V^T was
shown to be local (decays with element distance). SCDM (Selected Columns of
the Density Matrix; QR-with-pivoting on P) picks nd columns of P that are
(a) localized and (b) span the null space -- a principled local basis.

Checks:
  1. span     : ||(I - Pi_B) V||_F / ||V||_F  ~ 0   (B spans the null space)
  2. locality : mean fraction of each basis vector's mass within 1.5 cells
  3. removal  : projection deflation (exact, field-preserving) and
                eigenvalue shift A_s = A + alpha B (Bw^T B)^-1 Bw^T
  4. solver   : plain BiCGSTAB on A (ill-conditioned) vs on A_s (clean)

Usage: python nullmode_local_deflation.py [nx ny nz] [mu_r]   default 4 4 2 1e5
"""
import sys
import numpy as np
from scipy.linalg import qr
from scipy.sparse.linalg import bicgstab, LinearOperator
import radia as rad

EDGE = 0.01


def build(nx, ny, nz, mu_r):
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    objs, cent = [], []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                x0, y0, z0 = ix * EDGE, iy * EDGE, iz * EDGE
                v = [[x0, y0, z0], [x0+EDGE, y0, z0], [x0+EDGE, y0+EDGE, z0], [x0, y0+EDGE, z0],
                     [x0, y0, z0+EDGE], [x0+EDGE, y0, z0+EDGE], [x0+EDGE, y0+EDGE, z0+EDGE], [x0, y0+EDGE, z0+EDGE]]
                o = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(o, mat); objs.append(o)
                cent.append([(ix+0.5)*EDGE, (iy+0.5)*EDGE, (iz+0.5)*EDGE])
    return rad.ObjCnt(objs), np.asarray(cent)


def scdm(P, k):
    """Selected Columns of the Density Matrix: QR-pivot pick k local columns."""
    _, _, piv = qr(P, pivoting=True, mode="economic")
    cols = piv[:k]
    return P[:, cols].copy(), cols


def proj_onto(B):
    """orthogonal projector onto span(B) via normal equations (B kept local)."""
    G = B.T @ B
    Ginv = np.linalg.inv(G)
    return lambda x: B @ (Ginv @ (B.T @ x))


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

    blk, cent = build(nx, ny, nz, mu_r)
    N, dof = rad.GetInteractMatrix(rad.BuildMatrix(blk))
    N = np.asarray(N, float)
    A = np.diag(np.full(dof, inv_chi)) - N
    elem_of = np.arange(dof) // 6

    U, s, Vt = np.linalg.svd(N)
    null = s < 1e-9 * s[0]
    nd = int(null.sum())
    V = Vt[null].T.copy()        # right null space
    W = U[:, null].copy()        # left  null space
    print(f"{nx}x{ny}x{nz}={nx*ny*nz} elem, dof={dof}, mu_r={mu_r:.0e}, null-dim={nd}")

    # --- 1+2: local basis via SCDM on P = V V^T ---
    P = V @ V.T
    B, cols = scdm(P, nd)                 # local right-null basis (dof x nd)
    Pw = W @ W.T
    Bw, _ = scdm(Pw, nd)                  # local left-null basis

    PiB = proj_onto(B)
    span_err = np.linalg.norm(V - np.column_stack([PiB(V[:, j]) for j in range(nd)])) / np.linalg.norm(V)

    def locality(col):
        seed_elem = elem_of[cols[col]]
        m = np.array([np.sum(B[6*e:6*e+6, col]**2) for e in range(nx*ny*nz)])
        near = np.array([np.linalg.norm(cent[e]-cent[seed_elem]) <= 1.5*EDGE for e in range(nx*ny*nz)])
        return m[near].sum() / (m.sum() + 1e-30)
    loc = np.mean([locality(c) for c in range(nd)])
    print(f"local basis (SCDM): span error ||(I-Pi_B)V||/||V|| = {span_err:.2e}")
    print(f"                    mean mass within 1.5 cells of seed = {loc:.3f}  (1.0 = perfectly local)")

    # --- 3: removal ---
    rng = np.random.default_rng(0)
    b = rng.standard_normal(dof); b /= np.linalg.norm(b)
    sigma = np.linalg.solve(A, b)
    nf = lambda x: np.linalg.norm(V.T @ x) / (np.linalg.norm(x) + 1e-30)

    sigma_def = sigma - PiB(sigma)
    field_pres = np.linalg.norm(N @ sigma_def - N @ sigma) / (np.linalg.norm(N @ sigma) + 1e-30)
    print(f"\nconverged sigma:        null-frac = {nf(sigma):.3f}")
    print(f"projection (local B):   null-frac = {nf(sigma_def):.2e}, rel field change = {field_pres:.2e}")

    alpha = 1.0
    A_s = A + alpha * (B @ np.linalg.inv(Bw.T @ B) @ Bw.T)
    sigma_s = np.linalg.solve(A_s, b)
    print(f"shift A_s (local B,Bw): cond {np.linalg.cond(A):.2e} -> {np.linalg.cond(A_s):.2e}, "
          f"null-frac = {nf(sigma_s):.2e}")

    # --- 4: deflated BiCGSTAB iteration count ---
    def count(op, rhs, M=None):
        it = {"n": 0}
        bicgstab(op, rhs, rtol=1e-8, atol=0.0, maxiter=2000,
                 callback=lambda xk: it.__setitem__("n", it["n"]+1), M=M)
        return it["n"]
    Aop = LinearOperator((dof, dof), matvec=lambda x: A @ x)
    Asop = LinearOperator((dof, dof), matvec=lambda x: A_s @ x)
    it_plain = count(Aop, b)
    it_shift = count(Asop, b)
    print(f"\nBiCGSTAB iters: plain A = {it_plain},  shifted A_s = {it_shift}")
    print(f"(plain A converges to the loop-dominated 'ugly' solution; A_s to the clean one)")

    nzfrac = np.mean([np.sum(np.abs(B[:, c]) > 1e-6*np.max(np.abs(B[:, c]))) for c in range(nd)]) / dof
    print(f"\nSCALABILITY: nd={nd} ({100*nd/dof:.0f}% DOF); local basis avg nonzero support "
          f"~{100*nzfrac:.0f}% of DOF -> O(N) apply once localized via mesh cycles in C++.")
    rad.UtiDelAll()


if __name__ == "__main__":
    main()
