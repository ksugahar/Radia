"""Phase 3b — Krylov-Galerkin augmented CLN, replacing the broken Pade.

Paper 1 §III specifies:
    Y_CLN(s) = sigma * b_r^T K_r(s)^{-1} b_r
where K_r(s) = Q^T K(s) Q is the Galerkin-projected reduced matrix
with basis Q (an N x M matrix, N = FE DOFs, M = ROM size).

Construction of Q (paper 1 §II canonical CLN basis at s=0):
    q_0 = K_0^{-1} b                  (DC solution)
    q_k = K_0^{-1} (sigma * M * q_{k-1})  for k = 1..M-1
    Q = [q_0, q_1, ..., q_{M-1}], orthonormalised (modified Gram-Schmidt
        in the M-inner-product to preserve Cauer canonical form)

At runtime:
    K_r(s) = Q^T K_0 Q + s * sigma * Q^T M Q       (M x M complex matrix)
    b_r    = Q^T b                                  (M x 1)
    Y_CLN(s) = computed directly as the FE side dictates:
        Y_FE(s) = sigma*area - s*sigma^2 * b_const_r^T K_r(s)^{-1} b_const_r
        (volume-source split, exactly the same as Phase 3 attempt)

This avoids the Pade Hankel ill-conditioning entirely — all arithmetic
stays at the same numerical scale, and the M=4-5 ROM is well-conditioned
even for high-sigma problems.

Then the augmentation:
    Y_R(s) = Y_FE_via_ROM(s) + K_SIBC * sqrt(s) / (s + d)
"""

import math
import cmath
import numpy as np
import json
from pathlib import Path

from ngsolve import (
    Mesh, H1, BilinearForm, LinearForm, GridFunction,
    grad, dx, Integrate, x, y,
)
from netgen.geom2d import SplineGeometry
from scipy.special import iv


R     = 0.01
SIGMA = 1.0e6
MU0   = 4.0e-7 * math.pi
NU    = 1.0 / MU0
PERIM = 2 * math.pi * R
K_SIBC = PERIM * math.sqrt(SIGMA / MU0)
AREA  = math.pi * R * R


def make_disk_mesh(maxh):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), R, leftdomain=1, rightdomain=0, bc="outer")
    return Mesh(geo.GenerateMesh(maxh=maxh))


def bessel_Y(omega):
    """Closed-form Y(jw) for cylinder."""
    kappa = cmath.sqrt(1j * omega * SIGMA * MU0)
    kR    = kappa * R
    return SIGMA * 2 * math.pi * R * iv(1, kR) / (kappa * iv(0, kR))


def build_krylov_basis(K0_inv, M_mat, b_vec, fes, M_ROM):
    """Build CLN Krylov basis Q = [q_0, ..., q_{M_ROM-1}] (FE DOFs x M_ROM).

    q_0 = K0^{-1} b
    q_{k+1} = K0^{-1} (sigma * M * q_k)

    M-orthonormalisation (modified Gram-Schmidt with weight sigma*M)
    keeps the Cauer canonical form.
    """
    N = fes.ndof
    Q = np.zeros((N, M_ROM))
    q = b_vec.CreateVector()
    q.data = K0_inv * b_vec       # q_0
    tmp = b_vec.CreateVector()
    for k in range(M_ROM):
        v = np.array(q.FV().NumPy()).copy()
        # M-Gram-Schmidt against previous columns
        for j in range(k):
            # inner product <v, Q[:,j]>_{sigma*M}
            # We need v^T (sigma M) Q[:,j]. Use NGSolve mat-vec.
            qj_vec = b_vec.CreateVector()
            qj_vec.FV().NumPy()[:] = Q[:, j]
            tmp.data = M_mat * qj_vec
            inner = SIGMA * np.dot(v, np.array(tmp.FV().NumPy()))
            v -= inner * Q[:, j]
        # Normalise in sigma*M norm
        qv = b_vec.CreateVector()
        qv.FV().NumPy()[:] = v
        tmp.data = M_mat * qv
        nrm2 = SIGMA * np.dot(v, np.array(tmp.FV().NumPy()))
        if nrm2 <= 0:
            # Degenerate (e.g. M_ROM too large); break
            print(f"  Krylov breakdown at k={k}, nrm2={nrm2:.3e}", flush=True)
            return Q[:, :k]
        Q[:, k] = v / math.sqrt(nrm2)
        if k < M_ROM - 1:
            # Compute next q: K0^{-1} (sigma * M * q_k)
            qcurr = b_vec.CreateVector()
            qcurr.FV().NumPy()[:] = Q[:, k]
            tmp.data = M_mat * qcurr
            tmp.data *= SIGMA
            q.data = K0_inv * tmp
    return Q


def augmented_cln_krylov(mesh, M_ROM=5):
    """Galerkin-Krylov augmented CLN — replaces the broken Pade version."""
    fes = H1(mesh, order=1, dirichlet="outer")  # real-valued
    u, v_ = fes.TnT()
    a0 = BilinearForm(fes, symmetric=True)
    a0 += NU * grad(u) * grad(v_) * dx
    a0.Assemble()
    Mm = BilinearForm(fes, symmetric=True)
    Mm += u * v_ * dx
    Mm.Assemble()
    fL = LinearForm(fes)
    fL += 1.0 * v_ * dx
    fL.Assemble()

    K0_inv = a0.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    b_const = fL.vec
    b_const_np = np.array(b_const.FV().NumPy())

    # Build Krylov basis Q (FE DOFs x M_ROM)
    Q = build_krylov_basis(K0_inv, Mm.mat, b_const, fes, M_ROM)
    M_actual = Q.shape[1]
    print(f"  Krylov basis: {M_actual} columns built", flush=True)

    # Project K0, M, b
    # K0_r[i,j] = Q[:,i]^T K0 Q[:,j],  using NGSolve mat-vec
    K0_r = np.zeros((M_actual, M_actual))
    M_r  = np.zeros((M_actual, M_actual))
    tmp  = b_const.CreateVector()
    qj   = b_const.CreateVector()
    for j in range(M_actual):
        qj.FV().NumPy()[:] = Q[:, j]
        tmp.data = a0.mat * qj
        K0_r[:, j] = Q.T @ np.array(tmp.FV().NumPy())
        tmp.data = Mm.mat * qj
        M_r[:, j] = Q.T @ np.array(tmp.FV().NumPy())
    b_r = Q.T @ b_const_np

    # Find d: smallest generalized eigenvalue of (K0_r, sigma*M_r)
    # equivalent to smallest pole of Y_CLN
    from scipy.linalg import eigvals
    eigs = eigvals(K0_r, SIGMA * M_r)
    eigs_real_pos = [e.real for e in eigs if e.real > 0]
    d_from_eig = min(eigs_real_pos) if eigs_real_pos else 1.0e3

    # Compare with analytical Bessel first pole
    j01_sq = 2.40483 ** 2
    d_analytic = j01_sq / (SIGMA * MU0 * R * R)
    print(f"  smallest eig of (K0_r, sigma*M_r): {d_from_eig:.4e}",
          flush=True)
    print(f"  analytical 1st pole:               {d_analytic:.4e}",
          flush=True)

    # Use the eigenvalue-extracted d (intrinsic to the ROM)
    d = d_from_eig

    def Y_R(omega):
        s = 1j * omega
        # Reduced system: (K0_r + s*sigma*M_r) u_r = b_r
        K_r_s = K0_r + s * SIGMA * M_r
        u_r = np.linalg.solve(K_r_s, b_r)
        # Y_FE(s) = sigma*area - s*sigma^2 * b_const^T u_full
        #         = sigma*area - s*sigma^2 * b_r^T u_r   (since u_full = Q u_r)
        Y_cln = SIGMA * AREA - s * SIGMA * SIGMA * (b_r @ u_r)
        # sqrt(s) Schur-block augmentation
        Y_schur = K_SIBC * np.sqrt(s) / (s + d)
        return Y_cln + Y_schur

    return Y_R, {
        "M_ROM_actual": int(M_actual),
        "K0_r": K0_r.tolist(),
        "M_r":  M_r.tolist(),
        "b_r":  b_r.tolist(),
        "d_crossover": float(d),
        "d_analytical_wall_band": float(d_analytic),
        "K_SIBC": float(K_SIBC),
    }


# =====================================================================
# Main: sweep + compare with Bessel
# =====================================================================

def main():
    out_dir = Path(__file__).parent
    print(f"K_SIBC = {K_SIBC:.4e}, AREA = {AREA:.4e}", flush=True)

    mesh = make_disk_mesh(maxh=R / 3.5)
    print(f"coarse mesh: nv={mesh.nv}", flush=True)

    for M_ROM in [3, 4, 5, 6]:
        print()
        print(f"=== M_ROM = {M_ROM} ===", flush=True)
        Y_aug, diag = augmented_cln_krylov(mesh, M_ROM=M_ROM)

        # Sweep frequencies and compare with Bessel
        r_over_delta_list = [0.1, 0.3, 1.0, 3.0, 5.0, 10.0, 15.0]
        print(f"{'r/d':>5} {'|Y_bess|':>12} {'|Y_aug|':>12} {'err_aug':>10}",
              flush=True)
        for rd in r_over_delta_list:
            delta = R / rd
            omega = 2.0 / (delta * delta * SIGMA * MU0)
            Y_b = bessel_Y(omega)
            Y_a = Y_aug(omega)
            err = (abs(Y_a) - abs(Y_b)) / abs(Y_b) * 100
            print(f"{rd:5.2f} {abs(Y_b):12.4e} {abs(Y_a):12.4e} "
                  f"{err:+9.2f}%", flush=True)


if __name__ == "__main__":
    main()
