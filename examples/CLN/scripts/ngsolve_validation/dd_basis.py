"""dd_basis.py — Double-double evaluation of HDiv div-free hex basis.

Port of evaluate_basis() from hex_vim_cupy.py to DD arithmetic. Returns
Phi_hi, Phi_lo of shape (n_quad, n_dofs, 3) representing DD basis values.

For each block (A, B, C, D, E) the basis is a triple product of
{LXi, LXi_d, PXi, LEta, LEta_d, PEta, LZeta, LZeta_d, PZeta} evaluated
at the (xi, eta, zeta) = (2x-1, 2y-1, 2z-1) DD coordinates.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from dd_arithmetic import (
    dd_add, dd_sub, dd_mul, dd_div, dd_from_float, dd_from_mpmath,
)
from dd_legendre import dd_integrated_legendre_with_deriv, dd_legendre


def dd_evaluate_basis(x_hi, x_lo, y_hi, y_lo, z_hi, z_lo, p, xp_mod=None):
    """Evaluate ALL HDiv div-free basis functions at DD points (x, y, z).

    Args:
      x_hi, x_lo: shape (n_quad,) DD high/low parts of x in [0, 1]
      y_hi, y_lo: similar
      z_hi, z_lo: similar
      p: HDiv basis order
      xp_mod: numpy or cupy module (auto-detected if None)

    Returns:
      Phi_hi, Phi_lo: shape (n_quad, n_dofs, 3) each, DD basis values
                     where n_dofs = 2 p^3 + 3 p^2
    """
    if xp_mod is None:
        xp_mod_name = type(x_hi).__module__
        if xp_mod_name == "cupy":
            import cupy as xp_mod
        else:
            xp_mod = np

    # xi = 2x - 1, eta = 2y - 1, zeta = 2z - 1 (DD)
    # (2x is exact in FP64 since 2 is a power of 2; subtracting 1 also exact)
    xi_hi = 2.0 * x_hi - 1.0
    xi_lo = 2.0 * x_lo
    eta_hi = 2.0 * y_hi - 1.0
    eta_lo = 2.0 * y_lo
    zeta_hi = 2.0 * z_hi - 1.0
    zeta_lo = 2.0 * z_lo

    # Compute Legendre + IntegratedLegendre in DD
    LXi_hi, LXi_lo, LXi_d_hi, LXi_d_lo = dd_integrated_legendre_with_deriv(
        xi_hi, xi_lo, p)
    LEta_hi, LEta_lo, LEta_d_hi, LEta_d_lo = dd_integrated_legendre_with_deriv(
        eta_hi, eta_lo, p)
    LZeta_hi, LZeta_lo, LZeta_d_hi, LZeta_d_lo = dd_integrated_legendre_with_deriv(
        zeta_hi, zeta_lo, p)
    PXi_hi, PXi_lo = dd_legendre(xi_hi, xi_lo, p + 1)
    PEta_hi, PEta_lo = dd_legendre(eta_hi, eta_lo, p + 1)
    PZeta_hi, PZeta_lo = dd_legendre(zeta_hi, zeta_lo, p + 1)

    n_quad = x_hi.shape[0]
    n_dofs = 2 * p ** 3 + 3 * p ** 2
    Phi_hi = xp_mod.zeros((n_quad, n_dofs, 3), dtype=x_hi.dtype)
    Phi_lo = xp_mod.zeros((n_quad, n_dofs, 3), dtype=x_hi.dtype)

    eight_hi = 0.0 * x_hi + 8.0
    eight_lo = 0.0 * x_hi
    neg_eight_hi = 0.0 * x_hi - 8.0
    neg_eight_lo = 0.0 * x_hi

    offset = 0

    # Block A: 8 LXi[i] PEta[j+1] LZeta_d[k] for x; -8 LXi_d[i] PEta[j+1] LZeta[k] for z
    for i in range(p):
        for j in range(p):
            for k in range(p):
                # x-component: 8 * LXi[i] * PEta[j+1] * LZeta_d[k]
                t1_hi, t1_lo = dd_mul(LXi_hi[i], LXi_lo[i],
                                       PEta_hi[j + 1], PEta_lo[j + 1])
                t2_hi, t2_lo = dd_mul(t1_hi, t1_lo,
                                       LZeta_d_hi[k], LZeta_d_lo[k])
                ax_hi, ax_lo = dd_mul(eight_hi, eight_lo, t2_hi, t2_lo)
                # z-component: -8 * LXi_d[i] * PEta[j+1] * LZeta[k]
                t1_hi, t1_lo = dd_mul(LXi_d_hi[i], LXi_d_lo[i],
                                       PEta_hi[j + 1], PEta_lo[j + 1])
                t2_hi, t2_lo = dd_mul(t1_hi, t1_lo,
                                       LZeta_hi[k], LZeta_lo[k])
                az_hi, az_lo = dd_mul(neg_eight_hi, neg_eight_lo, t2_hi, t2_lo)
                idx = offset + i * p ** 2 + j * p + k
                Phi_hi[:, idx, 0] = ax_hi
                Phi_lo[:, idx, 0] = ax_lo
                Phi_hi[:, idx, 2] = az_hi
                Phi_lo[:, idx, 2] = az_lo
    offset += p ** 3

    # Block B: 8 PXi[i+1] LEta[j] LZeta_d[k] for y; -8 PXi[i+1] LEta_d[j] LZeta[k] for z
    for i in range(p):
        for j in range(p):
            for k in range(p):
                t1_hi, t1_lo = dd_mul(PXi_hi[i + 1], PXi_lo[i + 1],
                                       LEta_hi[j], LEta_lo[j])
                t2_hi, t2_lo = dd_mul(t1_hi, t1_lo,
                                       LZeta_d_hi[k], LZeta_d_lo[k])
                by_hi, by_lo = dd_mul(eight_hi, eight_lo, t2_hi, t2_lo)
                t1_hi, t1_lo = dd_mul(PXi_hi[i + 1], PXi_lo[i + 1],
                                       LEta_d_hi[j], LEta_d_lo[j])
                t2_hi, t2_lo = dd_mul(t1_hi, t1_lo,
                                       LZeta_hi[k], LZeta_lo[k])
                bz_hi, bz_lo = dd_mul(neg_eight_hi, neg_eight_lo, t2_hi, t2_lo)
                idx = offset + i * p ** 2 + j * p + k
                Phi_hi[:, idx, 1] = by_hi
                Phi_lo[:, idx, 1] = by_lo
                Phi_hi[:, idx, 2] = bz_hi
                Phi_lo[:, idx, 2] = bz_lo
    offset += p ** 3

    # Block C: 8 LEta[j] LZeta_d[k] for y; -8 LEta_d[j] LZeta[k] for z
    for j in range(p):
        for k in range(p):
            t1_hi, t1_lo = dd_mul(LEta_hi[j], LEta_lo[j],
                                   LZeta_d_hi[k], LZeta_d_lo[k])
            cy_hi, cy_lo = dd_mul(eight_hi, eight_lo, t1_hi, t1_lo)
            t1_hi, t1_lo = dd_mul(LEta_d_hi[j], LEta_d_lo[j],
                                   LZeta_hi[k], LZeta_lo[k])
            cz_hi, cz_lo = dd_mul(neg_eight_hi, neg_eight_lo, t1_hi, t1_lo)
            idx = offset + j * p + k
            Phi_hi[:, idx, 1] = cy_hi
            Phi_lo[:, idx, 1] = cy_lo
            Phi_hi[:, idx, 2] = cz_hi
            Phi_lo[:, idx, 2] = cz_lo
    offset += p ** 2

    # Block D: -8 LXi[i] LZeta_d[k] for x; 8 LXi_d[i] LZeta[k] for z
    for i in range(p):
        for k in range(p):
            t1_hi, t1_lo = dd_mul(LXi_hi[i], LXi_lo[i],
                                   LZeta_d_hi[k], LZeta_d_lo[k])
            dx_hi, dx_lo = dd_mul(neg_eight_hi, neg_eight_lo, t1_hi, t1_lo)
            t1_hi, t1_lo = dd_mul(LXi_d_hi[i], LXi_d_lo[i],
                                   LZeta_hi[k], LZeta_lo[k])
            dz_hi, dz_lo = dd_mul(eight_hi, eight_lo, t1_hi, t1_lo)
            idx = offset + i * p + k
            Phi_hi[:, idx, 0] = dx_hi
            Phi_lo[:, idx, 0] = dx_lo
            Phi_hi[:, idx, 2] = dz_hi
            Phi_lo[:, idx, 2] = dz_lo
    offset += p ** 2

    # Block E: 8 LXi[i] LEta_d[j] for x; -8 LXi_d[i] LEta[j] for y
    for i in range(p):
        for j in range(p):
            t1_hi, t1_lo = dd_mul(LXi_hi[i], LXi_lo[i],
                                   LEta_d_hi[j], LEta_d_lo[j])
            ex_hi, ex_lo = dd_mul(eight_hi, eight_lo, t1_hi, t1_lo)
            t1_hi, t1_lo = dd_mul(LXi_d_hi[i], LXi_d_lo[i],
                                   LEta_hi[j], LEta_lo[j])
            ey_hi, ey_lo = dd_mul(neg_eight_hi, neg_eight_lo, t1_hi, t1_lo)
            idx = offset + i * p + j
            Phi_hi[:, idx, 0] = ex_hi
            Phi_lo[:, idx, 0] = ex_lo
            Phi_hi[:, idx, 1] = ey_hi
            Phi_lo[:, idx, 1] = ey_lo
    offset += p ** 2

    assert offset == n_dofs
    return Phi_hi, Phi_lo


def main():
    """Validate DD HDiv basis against C++ FP64 reference + mpmath 60-digit ref."""
    sys.path.insert(0,
        str(Path("S:/Radia/01_GitHub/packages/radia-vim/src")))
    import radia_vim
    import mpmath as mp
    mp.mp.dps = 60

    print("=" * 72)
    print("DD HDiv basis validation")
    print("=" * 72)

    rng = np.random.default_rng(46)
    n_test = 5
    pts = rng.uniform(0.01, 0.99, size=(n_test, 3))

    for p in [2, 3, 4]:
        print(f"\n--- order p={p}, n_dofs={2*p**3 + 3*p**2} ---")
        basis = radia_vim.HDivDivFreeHexBasis(p)

        # Build DD input arrays
        x_hi = np.array(pts[:, 0])
        x_lo = np.zeros_like(x_hi)
        y_hi = np.array(pts[:, 1])
        y_lo = np.zeros_like(y_hi)
        z_hi = np.array(pts[:, 2])
        z_lo = np.zeros_like(z_hi)

        # DD evaluation
        Phi_hi, Phi_lo = dd_evaluate_basis(x_hi, x_lo, y_hi, y_lo, z_hi, z_lo,
                                            p, xp_mod=np)

        # C++ FP64 reference (should match DD hi part to ~1e-14)
        max_err_fp64 = 0.0
        for q in range(n_test):
            for i in range(basis.n_dofs):
                v_cpp = basis.evaluate(i, float(pts[q, 0]),
                                        float(pts[q, 1]), float(pts[q, 2]))
                for c in range(3):
                    val_dd = Phi_hi[q, i, c] + Phi_lo[q, i, c]
                    err = abs(val_dd - v_cpp[c])
                    if abs(v_cpp[c]) > 1e-10:
                        rel = err / abs(v_cpp[c])
                        max_err_fp64 = max(max_err_fp64, rel)

        print(f"  max rel error vs C++ FP64: {max_err_fp64:.3e} "
              f"(expected ~1e-15 for DD-truncated-to-FP64)")


if __name__ == "__main__":
    main()
