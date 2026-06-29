"""hex_vim_cupy_kassembly.py — GPU K matrix assembly via Spherical Duffy.

Vectorized CuPy port of compute_K_entry_spherical_duffy + assemble_K_bare.

Algorithm:
  1. Build outer (theta, phi, rho) quadrature: n_outer = n_th * n_ph * n_rh
  2. For each outer point: compute u = rho * (sin θ cos φ, sin θ sin φ, cos θ),
     c-box [|u_d|/2, 1 - |u_d|/2], jacobians, sphere weight
  3. Build inner (cx, cy, cz) Gauss-Legendre on [0, 1]^3 chunked per outer point
  4. Total quadrature points: n_q = n_outer * n_c^3
  5. Compute (rx, ry, rz) = c + u/2 and (rpx, rpy, rpz) = c - u/2 at all q
  6. Evaluate Phi_left = basis(r) and Phi_right = basis(r')  shape (n_q, n_dofs, 3)
  7. K[i,j] = sum_q sum_c diag_c * w_q * Phi_left[q,i,c] * Phi_right[q,j,c]
           = sum_c diag_c * (Phi_left[:,:,c] * sqrt(w))^T @ (Phi_right[:,:,c] * sqrt(w))

  Wait, w can be negative (it's not, all positive in our case), so use:
  K[i,j] = sum_c diag_c * (Phi_left[:,:,c]).T @ (w * Phi_right[:,:,c])
"""
from __future__ import annotations

import sys
import time
import math
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from hex_vim_cupy import evaluate_basis, integrated_legendre_with_deriv_np, legendre_np


def gauss_legendre_m1_1(n):
    """Gauss-Legendre nodes/weights on [-1, 1] (numpy, transferred to GPU later)."""
    nodes, weights = np.polynomial.legendre.leggauss(n)
    return nodes, weights


def assemble_K_cupy(p, a, b, c, n_th=12, n_ph=24, n_rh=12, n_c=6,
                    use_gpu=True, dtype=np.float64, verbose=True):
    """Assemble bare K matrix on GPU (no mu0/4pi prefactor).

    Args:
      p: HDiv basis order
      a, b, c: cuboid dimensions [m]
      n_th, n_ph, n_rh, n_c: quadrature orders
      use_gpu: if True use CuPy, else numpy
      dtype: float32 or float64
    """
    if use_gpu:
        import cupy as xp
        device = "GPU"
    else:
        xp = np
        device = "CPU"

    n_outer = n_th * n_ph * n_rh
    n_q = n_outer * n_c ** 3
    n_dofs = 2 * p ** 3 + 3 * p ** 2

    if verbose:
        print(f"  [{device}] order={p}, n_dofs={n_dofs}")
        print(f"  [{device}] quad: th={n_th}, ph={n_ph}, rh={n_rh}, c={n_c}")
        print(f"  [{device}] n_outer={n_outer}, n_q={n_q}")
        mem_gb = 2 * n_q * n_dofs * 3 * (8 if dtype == np.float64 else 4) / 1e9
        print(f"  [{device}] est mem (Phi_L + Phi_R): {mem_gb:.2f} GB")

    # ----- Build outer quadrature (theta, phi, rho) ------------------------
    th_n, th_w = gauss_legendre_m1_1(n_th)
    ph_n, ph_w = gauss_legendre_m1_1(n_ph)
    rh_n, rh_w = gauss_legendre_m1_1(n_rh)
    c_n, c_w   = gauss_legendre_m1_1(n_c)

    # Move to device
    th_n = xp.asarray(th_n, dtype=dtype)
    th_w = xp.asarray(th_w, dtype=dtype)
    ph_n = xp.asarray(ph_n, dtype=dtype)
    ph_w = xp.asarray(ph_w, dtype=dtype)
    rh_n = xp.asarray(rh_n, dtype=dtype)
    rh_w = xp.asarray(rh_w, dtype=dtype)
    c_n  = xp.asarray(c_n,  dtype=dtype)
    c_w  = xp.asarray(c_w,  dtype=dtype)

    pi = math.pi
    # theta in [0, pi]: theta = (th_n + 1) * pi/2
    theta = (th_n + 1.0) * (pi / 2.0)              # (n_th,)
    jac_th = pi / 2.0
    s_th = xp.sin(theta)                            # (n_th,)
    c_th = xp.cos(theta)                            # (n_th,)

    # phi in [0, 2 pi]: phi = (ph_n + 1) * pi
    phi = (ph_n + 1.0) * pi                        # (n_ph,)
    jac_ph = pi
    c_ph = xp.cos(phi)                              # (n_ph,)
    s_ph = xp.sin(phi)                              # (n_ph,)

    # omega = (s_th cos_ph, s_th sin_ph, c_th)  shape (n_th, n_ph, 3)
    omega_x = xp.einsum("t,p->tp", s_th, c_ph)       # (n_th, n_ph)
    omega_y = xp.einsum("t,p->tp", s_th, s_ph)       # (n_th, n_ph)
    omega_z = xp.broadcast_to(c_th[:, None], (n_th, n_ph)).copy()

    # omega_max = max |omega_d|
    omega_max = xp.maximum(xp.maximum(xp.abs(omega_x), xp.abs(omega_y)),
                           xp.abs(omega_z))                     # (n_th, n_ph)
    # Mask away tiny omega_max (set rho_max=0, contribution=0)
    eps = 1e-30
    rho_max = xp.where(omega_max > eps, 1.0 / omega_max, 0.0)   # (n_th, n_ph)
    valid_omega = omega_max > eps                                # (n_th, n_ph)

    # |Lambda omega| = sqrt(a² ωx² + b² ωy² + c² ωz²)
    Lomega = xp.sqrt(a*a*omega_x**2 + b*b*omega_y**2 + c*c*omega_z**2)
    Lomega = xp.where(valid_omega, Lomega, 1.0)  # avoid /0

    # sphere weight: th_w * ph_w * jac_th * jac_ph * sin_th / Lomega
    sphere_w = (xp.einsum("t,p->tp", th_w, ph_w)
                * jac_th * jac_ph
                * s_th[:, None] / Lomega)             # (n_th, n_ph)
    sphere_w = xp.where(valid_omega, sphere_w, 0.0)

    # ----- rho integration ---------------------------------------------------
    # rho_val[t, p, r] = (rh_n[r] + 1) * rho_max[t,p] / 2
    # jac_rh[t, p] = rho_max[t, p] / 2
    # rho_w[r] applied
    jac_rh = rho_max / 2.0                            # (n_th, n_ph)
    rho_val = (rh_n[None, None, :] + 1.0) * jac_rh[..., None]  # (n_th, n_ph, n_rh)

    # u[t,p,r,d] = rho * omega
    u_x = rho_val * omega_x[..., None]                 # (n_th, n_ph, n_rh)
    u_y = rho_val * omega_y[..., None]
    u_z = rho_val * omega_z[..., None]

    # c-box: cx_lo = |u_x|/2, cx_hi = 1 - |u_x|/2, cx_jac = (cx_hi - cx_lo)/2 = (1-|u_x|)/2
    abs_ux = xp.abs(u_x)
    abs_uy = xp.abs(u_y)
    abs_uz = xp.abs(u_z)
    cx_lo = abs_ux / 2.0
    cy_lo = abs_uy / 2.0
    cz_lo = abs_uz / 2.0
    cx_jac = (1.0 - abs_ux) / 2.0  # if negative, region empty
    cy_jac = (1.0 - abs_uy) / 2.0
    cz_jac = (1.0 - abs_uz) / 2.0

    # valid_c: all jacs > 0
    valid_c = (cx_jac > 0) & (cy_jac > 0) & (cz_jac > 0)
    c_volume_jac = xp.where(valid_c, cx_jac * cy_jac * cz_jac, 0.0)

    # ----- Outer weight (per (t, p, r)) -------------------------------------
    # outer_w[t, p, r] = sphere_w[t, p] * rh_w[r] * jac_rh[t, p] * rho_val * c_volume_jac
    outer_w = (sphere_w[..., None]
               * rh_w[None, None, :]
               * jac_rh[..., None]
               * rho_val
               * c_volume_jac)                       # (n_th, n_ph, n_rh)

    # ----- Inner c quadrature ------------------------------------------------
    # cx = cx_lo + (c_n + 1) * cx_jac      shape (..., n_c)
    # similarly cy, cz
    cx = cx_lo[..., None] + (c_n[None, None, None, :] + 1.0) * cx_jac[..., None]
    cy = cy_lo[..., None] + (c_n[None, None, None, :] + 1.0) * cy_jac[..., None]
    cz = cz_lo[..., None] + (c_n[None, None, None, :] + 1.0) * cz_jac[..., None]
    # Shapes: (n_th, n_ph, n_rh, n_c)

    # r and r' at each quad point (t, p, r, ix, iy, iz):
    # rx = cx[..., ix] + u_x/2,  ry = cy[..., iy] + u_y/2, rz = cz[..., iz] + u_z/2
    # rpx = cx[..., ix] - u_x/2, etc.

    # Shape outer = (n_th, n_ph, n_rh)
    outer_shape = (n_th, n_ph, n_rh)
    # cx shape (n_th, n_ph, n_rh, n_c) — broadcast for outer-product over (ix, iy, iz)

    # rx[t,p,r,ix,iy,iz] = cx[t,p,r,ix] + u_x[t,p,r]/2
    # Build (n_th, n_ph, n_rh, n_c, 1, 1) etc., then broadcast to full 6D shape.
    target_shape = (n_th, n_ph, n_rh, n_c, n_c, n_c)
    half_ux = (u_x / 2.0)[..., None, None, None]  # (n_th, n_ph, n_rh, 1, 1, 1)
    half_uy = (u_y / 2.0)[..., None, None, None]
    half_uz = (u_z / 2.0)[..., None, None, None]
    rx = xp.broadcast_to(cx[..., :, None, None] + half_ux, target_shape)
    ry = xp.broadcast_to(cy[..., None, :, None] + half_uy, target_shape)
    rz = xp.broadcast_to(cz[..., None, None, :] + half_uz, target_shape)
    rpx = xp.broadcast_to(cx[..., :, None, None] - half_ux, target_shape)
    rpy = xp.broadcast_to(cy[..., None, :, None] - half_uy, target_shape)
    rpz = xp.broadcast_to(cz[..., None, None, :] - half_uz, target_shape)

    # Inner weight: c_w[ix] * c_w[iy] * c_w[iz]
    cw_xyz = (c_w[None, None, None, :, None, None]
              * c_w[None, None, None, None, :, None]
              * c_w[None, None, None, None, None, :])
    # Combined weight per quadrature point (n_th, n_ph, n_rh, n_c, n_c, n_c):
    full_w = outer_w[..., None, None, None] * cw_xyz   # broadcasts to target_shape

    # Flatten to (n_q,) — must materialize broadcasted views
    rx_flat = xp.ascontiguousarray(rx).reshape(n_q)
    ry_flat = xp.ascontiguousarray(ry).reshape(n_q)
    rz_flat = xp.ascontiguousarray(rz).reshape(n_q)
    rpx_flat = xp.ascontiguousarray(rpx).reshape(n_q)
    rpy_flat = xp.ascontiguousarray(rpy).reshape(n_q)
    rpz_flat = xp.ascontiguousarray(rpz).reshape(n_q)
    w_flat = xp.ascontiguousarray(xp.broadcast_to(full_w, target_shape)).reshape(n_q)

    # Filter out invalid points (where outer_w=0). We can keep them (zero weight)
    # but they waste memory; for now keep all.

    if verbose:
        t0 = time.time()
    # ----- Evaluate basis at r and r' --------------------------------------
    Phi_L = evaluate_basis(rx_flat, ry_flat, rz_flat, p, xp_mod=xp)  # (n_q, n_dofs, 3)
    if verbose:
        if use_gpu:
            xp.cuda.Stream.null.synchronize()
        print(f"  [{device}] Phi_L assembled in {time.time()-t0:.2f}s")
        t0 = time.time()
    Phi_R = evaluate_basis(rpx_flat, rpy_flat, rpz_flat, p, xp_mod=xp)
    if verbose:
        if use_gpu:
            xp.cuda.Stream.null.synchronize()
        print(f"  [{device}] Phi_R assembled in {time.time()-t0:.2f}s")

    # ----- K[i,j] via einsum -----------------------------------------------
    # K[i,j] = sum_q sum_c (a²,b²,c²)[c] * w[q] * Phi_L[q,i,c] * Phi_R[q,j,c]
    # Compute per c-component matmul, sum.
    if verbose:
        t0 = time.time()
    diag_abc = xp.asarray([a*a, b*b, c*c], dtype=dtype)
    K = xp.zeros((n_dofs, n_dofs), dtype=dtype)
    for cc in range(3):
        # weighted Phi_R[q, j, cc] -> (n_q, n_dofs)
        Pr_w = Phi_R[..., cc] * w_flat[:, None]  # (n_q, n_dofs)
        # K_c[i, j] = Phi_L[:, i, cc].T @ Pr_w
        K_c = Phi_L[..., cc].T @ Pr_w
        K += diag_abc[cc] * K_c
    if verbose:
        if use_gpu:
            xp.cuda.Stream.null.synchronize()
        print(f"  [{device}] K einsum in {time.time()-t0:.2f}s")

    # Match C++ assemble_K_bare semantics: upper-triangle compute, mirror.
    # Spherical Duffy with GL quadrature on (theta, phi) ∈ [0,pi]x[0,2pi]
    # is NOT symmetric under omega↔-omega in the quadrature nodes (it
    # mirrors via (theta, phi) ↔ (pi-theta, 2pi-phi) instead of
    # (pi-theta, phi+pi)). Therefore compute_K_entry(i,j) ≠ compute_K_entry(j,i)
    # at finite quadrature. C++ uses upper triangle and mirrors. For
    # bit-exact agreement we do the same: take upper triangle, mirror.
    if use_gpu:
        K_np = xp.asnumpy(K)
    else:
        K_np = K
    iu = np.triu_indices_from(K_np, k=1)
    K_np[iu[1], iu[0]] = K_np[iu[0], iu[1]]  # mirror upper to lower
    return K_np


# =============================================================================
# Validation: GPU K vs C++ K for small case
# =============================================================================

def validate_K(p, a, b, c, n_th=8, n_ph=12, n_rh=8, n_c=4,
               use_gpu=True, atol_rel=1e-6):
    """Compare CuPy K vs C++ assemble_K_bare for given quadrature settings."""
    sys.path.insert(0,
        str(Path("S:/Radia/01_GitHub/src/ext/radia_vim/src")))
    import radia_vim
    basis = radia_vim.HDivDivFreeHexBasis(p)
    n_dofs = basis.n_dofs

    # CPU/GPU assembly
    print(f"\n=== Validate K, order={p}, n_dofs={n_dofs} ===")
    print(f"  Quad: th={n_th}, ph={n_ph}, rh={n_rh}, c={n_c}")
    t0 = time.time()
    K_gpu = assemble_K_cupy(p, a, b, c, n_th, n_ph, n_rh, n_c,
                             use_gpu=use_gpu, verbose=True)
    t_gpu = time.time() - t0
    print(f"  Total GPU/CPU time: {t_gpu:.2f}s")

    # C++ reference
    print(f"  Computing C++ reference (single-thread)...")
    t0 = time.time()
    K_cpp = radia_vim.assemble_K_bare(basis, a, b, c,
                                       n_omega_theta=n_th,
                                       n_omega_phi=n_ph,
                                       n_rho=n_rh,
                                       n_c=n_c,
                                       verbose=0)
    t_cpp = time.time() - t0
    print(f"  C++ time: {t_cpp:.2f}s")

    # Compare. Filter to entries with |K_cpp| > 1e-15 to ignore rounding noise.
    significant = np.abs(K_cpp) > 1e-15
    diff = np.abs(K_gpu - K_cpp)
    max_abs = float(np.max(diff))
    abs_K_cpp_sig = np.where(significant, np.abs(K_cpp), 1.0)
    rel_sig = np.where(significant, diff / abs_K_cpp_sig, 0.0)
    max_rel_sig = float(np.max(rel_sig))
    print(f"  max|K_gpu - K_cpp|             = {max_abs:.3e}")
    print(f"  max rel diff (|K_cpp|>1e-15)   = {max_rel_sig:.3e}")
    print(f"  K_cpp range: [{np.min(K_cpp):.3e}, {np.max(K_cpp):.3e}]")
    print(f"  K_cpp diag range: [{np.min(np.diag(K_cpp)):.3e}, {np.max(np.diag(K_cpp)):.3e}]")
    print(f"  Top 5 abs errors:")
    flat_diff = diff.flatten()
    top5 = np.argsort(-flat_diff)[:5]
    for idx in top5:
        i, j = np.unravel_index(idx, diff.shape)
        print(f"    K[{i},{j}]: gpu={K_gpu[i,j]:.6e}  cpp={K_cpp[i,j]:.6e}  "
              f"diff={diff[i,j]:.3e}  rel={diff[i,j]/abs(K_cpp[i,j]+1e-30):.2e}")
    print(f"  speedup vs C++:                = {t_cpp/t_gpu:.1f}x")

    if max_rel_sig < atol_rel:
        print(f"  PASS (rel tol = {atol_rel})")
    else:
        print(f"  FAIL (rel tol = {atol_rel}) — likely quadrature precision")
    return max_rel_sig, t_gpu, t_cpp


def main():
    print("=" * 72)
    print("CuPy Spherical Duffy K assembly validation (Phase 2)")
    print("=" * 72)

    # Cuboid 5x2x1 mm
    a, b, c = 5e-3, 2e-3, 1e-3

    # Production settings (n_th=12, n_ph=24, n_rh=12, n_c=6) + order=4
    print("\n### GPU @ production quad (12,24,12,6) order=4 ###")
    print("    (C++ reference takes ~4075 s; expected GPU ~3-5 s)")
    validate_K(p=4, a=a, b=b, c=c,
               n_th=12, n_ph=24, n_rh=12, n_c=6,
               use_gpu=True, atol_rel=1e-6)


if __name__ == "__main__":
    main()
