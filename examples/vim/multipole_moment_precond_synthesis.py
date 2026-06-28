"""COMPLETE synthesis: why no CHEAP preconditioner bounds multipole-moment MMM GMRES iters at high mu_r, and what the
loops actually do.  This reconciles the whole 2026-06-22 investigation (refines multipole_moment_precond_diagnosis.py
+ multipole_moment_block_size_test.py) and resolves the user's loop / block / center-charge hypotheses.

Chain of decisive facts (C-yoke, nxy=16, per-element block-Jacobi baseline = 174 GMRES iters at mu_r=1000):

1. LOOPS ARE FIELD-NULL IN THE DEMAG OPERATOR: ||D @ Lb|| / ||D|| = 5e-16 (D = the centroid field+gradient
   operator; Lb = rad.GetLoopBasis).  So A @ loop = L @ loop (local part only), chi-INDEPENDENT.

2. BLOCK-JACOBI OVER-SCALES THE LOOPS BY chi: the preconditioner M = block-diag(A) = block-diag(L - chi D)
   carries the chi-scaled self-field (D's diagonal blocks != 0), while A @ loop is chi-independent.  So
   M^-1 A @ loop ~ 1/chi -> the loop modes collapse to eigenvalue ~1/chi.  MEASURED: #(|eig(M^-1A)|<0.05)
   grows 0 -> 124 -> 205 for mu_r 2 -> 100 -> 1000 (205 = nLoop), min|eig| = 1.77e-3 ~ 1/999 = 1/chi, and the
   near-zero modes have loop-overlap = 1.000.  THIS is the chi-scaling of the iteration count (15->174).

3. BUT THE NEAR-ZERO LOOP CLUSTER IS NOT THE GMRES BOTTLENECK.  The operator is NON-NORMAL
   (||A A^T - A^T A|| / ||A||^2 = 7e-2), so GMRES convergence is governed by the pseudospectrum, NOT the
   eigenvalues.  TRUE operator-deflation of ALL 205 loop modes ((I - A Q) A, Q = Lb (Lb^T A Lb)^-1 Lb^T)
   improves only 174 -> 160.  Removing the entire near-zero cluster barely helps -> eigenvalues do not drive
   the count here.  The ~160 that remain are the NON-NORMAL DEMAG-COMPLEMENT -- the real bottleneck (this is
   what multipole_moment_precond_diagnosis.py correctly called "the demag mutual coupling").

VERDICT.  Every CHEAP lever fails to bound iters, each for the reason above:
  - per-element / bigger / flux-path-reordered BLOCKS: block_size_test -- loop-range HURTS (174->219), only
    near-global blocks help (non-local stiff structure);
  - net-charge / center-charge handling: precond_diagnosis -- slow modes have 0.000 net-charge overlap;
  - eigenvalue or loop-space DEFLATION: this script -- non-normal, so deflating eigenmodes barely helps.
The robust cure is an approximate FACTORIZATION that addresses the non-normal structure as a whole = H-LU on
the cluster tree (heavy A-build, but bounded iters -- the lab's standing conclusion, see HDiv-VIM memory).  If
H-LU is too heavy, the practical path is to ACCEPT the iter growth (it is ~dof^1.06, sub-linear; each iter is
a cheap HACApK matvec), exactly what the production surface-charge MSC / HDiv-VIM do at scale.  There is no free cheap
preconditioner here -- the obstacle is the non-normal high-mu_r demag operator, not loops / div(B) / blocks.

mdx-clean import (radia direct).
"""
import json
import os
import platform
from datetime import datetime

import numpy as np
import scipy.sparse.linalg as spla

import radia as rad

HERE = os.path.dirname(os.path.abspath(__file__))


def _inside_cyoke(cx, cy):
    return (-0.06 <= cx <= 0.06) and (-0.06 <= cy <= 0.06) and not (-0.035 < cx < 0.035 and -0.035 < cy < 0.035) and not (cx > 0.018)


def build_cyoke_hexes(nxy, nz):
    xs = np.linspace(-0.06, 0.06, nxy + 1); zs = np.linspace(-0.02, 0.02, nz + 1)
    hexes = []
    for k in range(nz):
        for j in range(nxy):
            for i in range(nxy):
                if not _inside_cyoke(0.5 * (xs[i] + xs[i + 1]), 0.5 * (xs[j] + xs[j + 1])):
                    continue
                x0, x1, y0, y1, z0, z1 = xs[i], xs[i + 1], xs[j], xs[j + 1], zs[k], zs[k + 1]
                hexes.append([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                              [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    return hexes


def _norm(row, rhs):
    nn = np.linalg.norm(row)
    return (row / nn, rhs / nn) if nn > 1e-300 else (row, rhs)


def assemble(hexes, Happ, mu_r):
    """Dense linear moment system A,b + per-element row blocks + loop basis Lb + the nonlocal operator D."""
    chi = mu_r - 1.0
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in hexes]
    for h in objs:
        rad.MatApl(h, rad.MatLin(mu_r))
    handle = rad.BuildMatrix(rad.ObjCnt(objs))
    G = np.asarray(rad.GetFaceGeom(handle), float); C = np.asarray(rad.GetCentroidFieldGrad(handle), float)
    Lb_raw, nLoop = rad.GetLoopBasis(handle); Lb = np.asarray(Lb_raw, float)
    dof = G.shape[0]; elem = G[:, 0].astype(int); area = G[:, 1]; fc = G[:, 2:5]; nrm = G[:, 5:8]; ecen = G[:, 8:11]
    n_el = int(elem.max()) + 1; dofs_of = [np.where(elem == e)[0] for e in range(n_el)]
    if Lb.size and Lb.shape[0] != dof and Lb.shape[-1] == dof:
        Lb = Lb.T
    A = np.zeros((dof, dof)); b = np.zeros(dof); r = 0; row_of = []
    for e in range(n_el):
        fs = dofs_of[e]; Ae = area[fs]; ne = nrm[fs]; ce = ecen[fs[0]]; d = fc[fs] - ce
        Ve = (1.0 / 3.0) * np.sum(Ae * np.sum(fc[fs] * ne, axis=1)); F0 = C[e, 0:3, :]; Ginv = C[e, 3:9, :]
        row_of.append((r, fs)); dip = np.zeros((3, dof))
        for k in range(3):
            dip[k, fs] = Ae * d[:, k]
        for k in range(3):
            row, rhs = _norm(dip[k, :] / Ve - chi * F0[k, :], chi * Happ[k]); A[r, :] = row; b[r] = rhs; r += 1
        mono = np.zeros(dof); mono[fs] = Ae
        row, rhs = _norm(mono, 0.0); A[r, :] = row; b[r] = rhs; r += 1
        for Bm in (d[:, 0]**2 - d[:, 1]**2, d[:, 1]**2 - d[:, 2]**2):
            row = np.zeros(dof); row[fs] = Ae * Bm
            cm = np.array([np.sum(Ae * ne[:, k] * Bm) for k in range(3)]); row -= (cm @ dip) / Ve
            Dm = np.array([[np.sum(Ae * d[:, jj] * ne[:, ii] * Bm) for jj in range(3)] for ii in range(3)])
            Dvec = np.array([Dm[0, 0], Dm[1, 1], Dm[2, 2], Dm[0, 1]+Dm[1, 0], Dm[0, 2]+Dm[2, 0], Dm[1, 2]+Dm[2, 1]])
            row -= chi * (Dvec @ Ginv)
            row, rhs = _norm(row, 0.0); A[r, :] = row; b[r] = rhs; r += 1
    Dop = np.vstack([C[:, 0:3, :].reshape(3 * n_el, dof), C[:, 3:9, :].reshape(6 * n_el, dof)])
    rad.UtiDelAll()
    return A, b, row_of, dof, n_el, Lb, Dop


def block_jacobi(A, row_of, dof):
    blocks = [(r0, fs, np.linalg.inv(A[r0:r0 + 6, fs])) for (r0, fs) in row_of]
    def apply(rvec):
        rvec = np.asarray(rvec, float); out = np.zeros(dof)
        for (r0, fs, Mi) in blocks:
            out[fs] += Mi @ rvec[r0:r0 + 6]
        return out
    return apply


def gmres_count(matvec, b, dof, M=None):
    A_op = spla.LinearOperator((dof, dof), matvec=matvec)
    M_op = spla.LinearOperator((dof, dof), matvec=M) if M is not None else None
    nit = [0]
    x, info = spla.gmres(A_op, b, M=M_op, rtol=1e-10, atol=0.0, restart=min(dof, 300), maxiter=3000,
                         callback=lambda rk: nit.__setitem__(0, nit[0] + 1), callback_type="pr_norm")
    return x, nit[0], info


def main():
    nxy = 16
    hexes = build_cyoke_hexes(nxy, 2); Happ = np.array([0.0, 1e3, 0.0])

    # --- mechanism: D@Lb=0, near-zero cluster grows with chi as 1/chi, modes are loops ---
    print(f"\nmultipole-moment MMM preconditioner SYNTHESIS (C-yoke nxy={nxy}).\n")
    rows_mech = []
    A0, b0, row_of0, dof, n_el, Lb0, Dop0 = assemble(hexes, Happ, 1000.0)
    Qloop, _ = np.linalg.qr(Lb0)
    dlb = float(np.linalg.norm(Dop0 @ Qloop) / np.linalg.norm(Dop0))
    nonnormal = float(np.linalg.norm(A0 @ A0.T - A0.T @ A0) / np.linalg.norm(A0) ** 2)
    print(f"  1. loops field-null in D:  ||D Lb||/||D|| = {dlb:.1e}   (=> A.loop = L.loop, chi-independent)")
    print(f"     operator non-normality: ||AA^T - A^T A||/||A||^2 = {nonnormal:.1e}   (0 = normal)\n")
    print(f"  2. block-Jacobi over-scales loops by chi -> near-zero cluster = loops:")
    print(f"     {'mu_r':>6} {'#|eig(M^-1A)|<0.05':>18} {'min|eig|':>10} {'~1/chi':>9} {'nearzero loop-overlap':>21}")
    for mu_r in (2.0, 100.0, 1000.0):
        A, b, row_of, dof, n_el, Lb, _ = assemble(hexes, Happ, mu_r)
        bj = block_jacobi(A, row_of, dof)
        MinvA = np.column_stack([bj(A[:, j]) for j in range(dof)])
        w, V = np.linalg.eig(MinvA); absA = np.abs(w)
        ncluster = int(np.sum(absA < 0.05)); mineig = float(absA.min())
        if ncluster > 0:
            Ql, _ = np.linalg.qr(Lb)
            near = np.argsort(absA)[:ncluster]
            ov = float(np.mean([np.linalg.norm(Ql.T @ (np.real(V[:, j]) / (np.linalg.norm(np.real(V[:, j])) + 1e-30))) for j in near]))
        else:
            ov = float("nan")
        print(f"     {mu_r:>6.0f} {ncluster:>18d} {mineig:>10.2e} {1.0/(mu_r-1.0):>9.2e} {ov:>21.3f}")
        rows_mech.append(dict(mu_r=mu_r, near_zero_count=ncluster, min_eig=mineig, inv_chi=1.0/(mu_r-1.0),
                              nearzero_loop_overlap=ov, nLoop=int(Lb.shape[1]) if Lb.size else 0))

    # --- the near-zero loop cluster is NOT the bottleneck: deflate all loops, barely helps ---
    bj0 = block_jacobi(A0, row_of0, dof)
    _, it_base, _ = gmres_count(lambda x: A0 @ x, b0, dof, M=bj0)
    Z, _ = np.linalg.qr(Lb0); AZ = A0 @ Z; Einv = np.linalg.inv(Z.T @ AZ)
    Pb = b0 - AZ @ (Einv @ (Z.T @ b0))
    def PA(x):
        Ax = A0 @ x; return Ax - AZ @ (Einv @ (Z.T @ Ax))
    xt, it_defl, _ = gmres_count(PA, Pb, dof, M=bj0)
    x = Z @ (Einv @ (Z.T @ b0)) + xt - Z @ (Einv @ (Z.T @ (A0 @ xt)))
    res = float(np.linalg.norm(A0 @ x - b0) / np.linalg.norm(b0))
    print(f"\n  3. near-zero LOOP cluster is NOT the bottleneck (non-normal -> eigenvalues don't drive GMRES):")
    print(f"     block-Jacobi GMRES          : {it_base} iters")
    print(f"     + TRUE loop-deflation (all {Z.shape[1]}): {it_defl} iters   (full-system res {res:.0e} -- correct solve)")
    print(f"     deflating the ENTIRE near-zero cluster barely helps -> the bottleneck is the non-normal")
    print(f"     demag-complement, not the loops.\n")

    out = dict(timestamp=datetime.now().isoformat(), hostname=platform.node(),
               benchmark="multipole_moment_precond_synthesis", nxy=nxy, n_el=n_el, dof=dof,
               D_dot_loop_rel=dlb, nonnormality=nonnormal, mechanism=rows_mech,
               iters_block_jacobi=it_base, iters_loop_deflated=it_defl, loop_deflation_residual=res,
               conclusion=(
                   "Loops are field-null in D (||D Lb||/||D||=5e-16), so block-Jacobi -- whose blocks carry the "
                   "chi-scaled self-field -- over-scales them, collapsing all nLoop loop modes to eigenvalue "
                   "~1/chi (near-zero count 0/124/205 for mu_r 2/100/1000, min|eig|~1/chi, loop-overlap 1.000). "
                   "That IS the chi-scaling of the iteration count. BUT the operator is non-normal (7e-2), so "
                   "GMRES is governed by the pseudospectrum, not eigenvalues: TRUE operator-deflation of all 205 "
                   f"loop modes only improves {it_base}->{it_defl}. The bottleneck is the non-normal demag-"
                   "complement (what precond_diagnosis correctly called the demag mutual coupling). Every cheap "
                   "lever fails -- blocks (loop-range hurts), net-charge handling (0 overlap), eigenvalue/loop "
                   "deflation (non-normal). Robust cure = H-LU on the cluster tree (heavy A-build, bounded iters, "
                   "the lab's standing conclusion); cheap alternative = ACCEPT the sub-linear iter growth "
                   "(~dof^1.06) with cheap HACApK matvecs, as the production surface-charge MSC / HDiv-VIM do. No free cheap "
                   "preconditioner exists here: the obstacle is the non-normal high-mu_r demag operator."))
    with open(os.path.join(HERE, "multipole_moment_precond_synthesis.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print("  VERDICT: cheap preconditioners (blocks / deflation / loop / center-charge) do NOT bound iters --")
    print("  the obstacle is the NON-NORMAL high-mu_r demag operator. H-LU (heavy, robust) or accept sub-linear")
    print("  iter growth with cheap HACApK matvecs (production surface-charge MSC / HDiv-VIM path).")
    print("  saved", os.path.join(HERE, "multipole_moment_precond_synthesis.json"))


if __name__ == "__main__":
    main()
