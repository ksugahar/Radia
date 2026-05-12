"""hex_vim_cupy.py — GPU-accelerated K matrix assembly for radia-vim hex VIM.

Vectorized CuPy port of:
  - HDivDivFreeHexBasis (Schöberl-Zaglmayr 5-block: A, B, C, D, E)
  - Spherical Duffy quadrature (sauter_schwab_hex.hpp)

Computes K[i,j] = ∫∫_{[0,1]^6} φ_i(r)·diag(a²,b²,c²)·φ_j(r')
                                / |diag(a,b,c)·(r-r')| dr dr'
on GPU via einsum, validated against C++ single-thread reference.

Memory budget for order=4, n_dofs=176, n_quad ≈ 750k:
  Phi_left + Phi_right ≈ 6 GB FP64 / 3 GB FP32. Fits in 16 GB Quadro RTX 5000.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np


def integrated_legendre_with_deriv_np(xi, p):
    """Vectorized IntegratedLegendre L_{j+2}(xi) and dL/dxi for j=0..p-1.

    Mirrors the C++ recurrence in hdiv_divfree_hex_basis.hpp::
      integrated_legendre_with_deriv. xi can be array of any shape.

    Returns:
      L:   shape xi.shape + (p,) with L[..., j] = L_{j+2}(xi)
      Lp:  shape xi.shape + (p,) with Lp[..., j] = dL_{j+2}/dxi
    """
    xp = type(xi).__module__  # 'numpy' or 'cupy'
    is_cupy = xp == "cupy"
    if is_cupy:
        import cupy as cp
        zeros = cp.zeros
        stack = cp.stack
    else:
        zeros = np.zeros
        stack = np.stack

    out_L = zeros(xi.shape + (p,), dtype=xi.dtype)
    out_Lp = zeros(xi.shape + (p,), dtype=xi.dtype)
    if p < 1:
        return out_L, out_Lp

    p3 = zeros(xi.shape, dtype=xi.dtype)
    p2 = zeros(xi.shape, dtype=xi.dtype) - 1.0
    p1 = xi.copy()
    p3p = zeros(xi.shape, dtype=xi.dtype)
    p2p = zeros(xi.shape, dtype=xi.dtype)
    p1p = zeros(xi.shape, dtype=xi.dtype) + 1.0

    for j in range(2, p + 2):
        a_coef = (2.0 * j - 3.0) / j
        b_coef = (j - 3.0) / j

        p3_new = p2.copy()
        p2_new = p1.copy()
        p3p_new = p2p.copy()
        p2p_new = p1p.copy()

        new_p1 = a_coef * xi * p2_new - b_coef * p3_new
        new_p1p = a_coef * (p2_new + xi * p2p_new) - b_coef * p3p_new

        out_L[..., j - 2] = new_p1
        out_Lp[..., j - 2] = new_p1p

        p3, p2, p1 = p3_new, p2_new, new_p1
        p3p, p2p, p1p = p3p_new, p2p_new, new_p1p

    return out_L, out_Lp


def legendre_np(xi, p):
    """Vectorized classical Legendre P_n(xi) for n=0..p.
    Returns shape xi.shape + (p+1,)."""
    xp = type(xi).__module__
    is_cupy = xp == "cupy"
    if is_cupy:
        import cupy as cp
        zeros = cp.zeros
    else:
        zeros = np.zeros

    P = zeros(xi.shape + (p + 1,), dtype=xi.dtype)
    P[..., 0] = 1.0
    if p >= 1:
        P[..., 1] = xi
    for k in range(2, p + 1):
        inv_k = 1.0 / k
        P[..., k] = ((2.0 - inv_k) * xi * P[..., k - 1]
                     + (inv_k - 1.0) * P[..., k - 2])
    return P


def evaluate_basis(x, y, z, p, xp_mod=None):
    """Evaluate ALL HDiv div-free basis functions at points (x, y, z).

    Args:
      x, y, z: arrays of same shape (n_quad,) in [0, 1]^3
      p:       order
      xp_mod:  numpy or cupy module (auto-detected if None)

    Returns:
      Phi: shape (n_quad, n_dofs, 3) where n_dofs = 2 p^3 + 3 p^2
           Phi[q, i, c] = c-component of basis function i at point q.
    """
    if xp_mod is None:
        xp_mod_name = type(x).__module__
        if xp_mod_name == "cupy":
            import cupy as xp_mod
        else:
            xp_mod = np

    xi = 2.0 * x - 1.0
    eta = 2.0 * y - 1.0
    zeta = 2.0 * z - 1.0

    LXi, LXi_d = integrated_legendre_with_deriv_np(xi, p)        # (n_quad, p)
    LEta, LEta_d = integrated_legendre_with_deriv_np(eta, p)
    LZeta, LZeta_d = integrated_legendre_with_deriv_np(zeta, p)
    PXi = legendre_np(xi, p + 1)                                  # (n_quad, p+2)
    PEta = legendre_np(eta, p + 1)
    PZeta = legendre_np(zeta, p + 1)

    n_quad = x.shape[0]
    n_dofs = 2 * p ** 3 + 3 * p ** 2
    Phi = xp_mod.zeros((n_quad, n_dofs, 3), dtype=x.dtype)

    eight = 8.0
    offset = 0

    # Block A: (i, j, k) ∈ [0..p-1]^3, p^3 dofs
    # S.x = +8 LXi[i] PEta[j+1] LZeta_d[k]
    # S.y = 0
    # S.z = -8 LXi_d[i] PEta[j+1] LZeta[k]
    # Outer-product flatten: idx = i*p^2 + j*p + k
    # We do this with einsum: outputs (n_quad, p, p, p)
    A_x = eight * xp_mod.einsum("qi,qj,qk->qijk", LXi, PEta[..., 1:p + 1], LZeta_d)
    A_z = -eight * xp_mod.einsum("qi,qj,qk->qijk", LXi_d, PEta[..., 1:p + 1], LZeta)
    Phi[:, offset:offset + p ** 3, 0] = A_x.reshape(n_quad, p ** 3)
    Phi[:, offset:offset + p ** 3, 2] = A_z.reshape(n_quad, p ** 3)
    offset += p ** 3

    # Block B: (i, j, k) ∈ [0..p-1]^3, p^3 dofs
    # S.x = 0
    # S.y = +8 PXi[i+1] LEta[j] LZeta_d[k]
    # S.z = -8 PXi[i+1] LEta_d[j] LZeta[k]
    B_y = eight * xp_mod.einsum("qi,qj,qk->qijk", PXi[..., 1:p + 1], LEta, LZeta_d)
    B_z = -eight * xp_mod.einsum("qi,qj,qk->qijk", PXi[..., 1:p + 1], LEta_d, LZeta)
    Phi[:, offset:offset + p ** 3, 1] = B_y.reshape(n_quad, p ** 3)
    Phi[:, offset:offset + p ** 3, 2] = B_z.reshape(n_quad, p ** 3)
    offset += p ** 3

    # Block C: (j, k) ∈ [0..p-1]^2, p^2 dofs (i unused/0)
    # S.x = 0
    # S.y = +8 LEta[j] LZeta_d[k]
    # S.z = -8 LEta_d[j] LZeta[k]
    C_y = eight * xp_mod.einsum("qj,qk->qjk", LEta, LZeta_d)
    C_z = -eight * xp_mod.einsum("qj,qk->qjk", LEta_d, LZeta)
    Phi[:, offset:offset + p ** 2, 1] = C_y.reshape(n_quad, p ** 2)
    Phi[:, offset:offset + p ** 2, 2] = C_z.reshape(n_quad, p ** 2)
    offset += p ** 2

    # Block D: (i, k) ∈ [0..p-1]^2, p^2 dofs (j unused/0)
    # S.x = -8 LXi[i] LZeta_d[k]
    # S.y = 0
    # S.z = +8 LXi_d[i] LZeta[k]
    D_x = -eight * xp_mod.einsum("qi,qk->qik", LXi, LZeta_d)
    D_z = eight * xp_mod.einsum("qi,qk->qik", LXi_d, LZeta)
    Phi[:, offset:offset + p ** 2, 0] = D_x.reshape(n_quad, p ** 2)
    Phi[:, offset:offset + p ** 2, 2] = D_z.reshape(n_quad, p ** 2)
    offset += p ** 2

    # Block E: (i, j) ∈ [0..p-1]^2, p^2 dofs (k unused/0)
    # S.x = +8 LXi[i] LEta_d[j]
    # S.y = -8 LXi_d[i] LEta[j]
    # S.z = 0
    E_x = eight * xp_mod.einsum("qi,qj->qij", LXi, LEta_d)
    E_y = -eight * xp_mod.einsum("qi,qj->qij", LXi_d, LEta)
    Phi[:, offset:offset + p ** 2, 0] = E_x.reshape(n_quad, p ** 2)
    Phi[:, offset:offset + p ** 2, 1] = E_y.reshape(n_quad, p ** 2)
    offset += p ** 2

    assert offset == n_dofs, f"DOF count mismatch: {offset} vs {n_dofs}"
    return Phi


# =============================================================================
# Validation: compare CuPy basis vs C++ basis at sample points
# =============================================================================

def validate_basis(p, n_test=20, atol=1e-12, use_cupy=False):
    """Compare evaluate_basis against radia_vim C++ at n_test random points."""
    sys.path.insert(0,
        str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))
    import radia_vim
    basis = radia_vim.HDivDivFreeHexBasis(p)
    n_dofs = basis.n_dofs
    print(f"  Validating order={p}, n_dofs={n_dofs}")

    rng = np.random.default_rng(42)
    pts = rng.uniform(0.01, 0.99, size=(n_test, 3))

    # C++ reference
    Phi_cpp = np.zeros((n_test, n_dofs, 3), dtype=np.float64)
    for q in range(n_test):
        x, y, z = pts[q]
        for i in range(n_dofs):
            v = basis.evaluate(i, float(x), float(y), float(z))
            Phi_cpp[q, i] = v

    # Python CPU
    if use_cupy:
        import cupy as cp
        x_arr = cp.asarray(pts[:, 0])
        y_arr = cp.asarray(pts[:, 1])
        z_arr = cp.asarray(pts[:, 2])
        Phi_py = evaluate_basis(x_arr, y_arr, z_arr, p)
        Phi_py = cp.asnumpy(Phi_py)
        device = "GPU"
    else:
        x_arr = pts[:, 0]
        y_arr = pts[:, 1]
        z_arr = pts[:, 2]
        Phi_py = evaluate_basis(x_arr, y_arr, z_arr, p)
        device = "CPU"

    diff = np.abs(Phi_py - Phi_cpp)
    max_abs = float(np.max(diff))
    max_rel = float(np.max(diff / (np.abs(Phi_cpp) + 1e-30)))
    print(f"  [{device}] max|Phi_py - Phi_cpp|     = {max_abs:.3e}")
    print(f"  [{device}] max relative diff         = {max_rel:.3e}")
    if max_abs < atol:
        print(f"  [{device}] PASS (tol={atol:.0e})")
    else:
        print(f"  [{device}] FAIL (tol={atol:.0e})")
        # Show worst case
        worst = np.argmax(diff)
        q, i, c = np.unravel_index(worst, diff.shape)
        print(f"    worst at (q,i,c)=({q},{i},{c}): "
              f"py={Phi_py[q,i,c]:.6e} cpp={Phi_cpp[q,i,c]:.6e}")
    return max_abs


def main():
    print("=" * 72)
    print("CuPy HDiv basis validation (Phase 1)")
    print("=" * 72)

    # Test order = 1, 2, 3, 4 on CPU first
    for p in [1, 2, 3, 4]:
        print(f"\n--- order = {p} ---")
        validate_basis(p, n_test=10, use_cupy=False)

    # Test on GPU
    print(f"\n--- order = 4 on GPU ---")
    validate_basis(4, n_test=10, use_cupy=True)


if __name__ == "__main__":
    main()
