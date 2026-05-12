"""dd_axisym_kernel.py — DD axisymmetric Green's function evaluation.

For axisym VIM, the kernel involves complete elliptic integrals K(m), E(m).
For DD precision, use mpmath ellipk/ellipe at 40-digit precision.

Per pair cost: ~1 ms (mpmath calls). For O(N²) pairs, this is the bottleneck.
Use multiprocessing for CPU parallelism (8 cores → ~8× speedup).
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import mpmath as mp

sys.path.insert(0, str(Path(__file__).parent))
from dd_arithmetic import (
    dd_add, dd_sub, dd_mul, dd_div, dd_from_mpmath, dd_sqrt,
)


def axisym_kernel_dd_mpmath(r_hi, r_lo, z_hi, z_lo,
                              r_p_hi, r_p_lo, z_p_hi, z_p_lo, dps=40):
    """Compute axisym Green's function G(r,z; r',z') in DD via mpmath.

    G = (1/π) sqrt(r'/r) [(2/k - k) K(m) - (2/k) E(m)]
    m = 4 r r' / ((r+r')² + (z-z')²),  k = sqrt(m)
    """
    mp.mp.dps = dps
    # Reconstruct full precision from DD
    r_mp = mp.mpf(float(r_hi)) + mp.mpf(float(r_lo))
    z_mp = mp.mpf(float(z_hi)) + mp.mpf(float(z_lo))
    r_p_mp = mp.mpf(float(r_p_hi)) + mp.mpf(float(r_p_lo))
    z_p_mp = mp.mpf(float(z_p_hi)) + mp.mpf(float(z_p_lo))

    if r_mp < mp.mpf("1e-30") or r_p_mp < mp.mpf("1e-30"):
        return 0.0, 0.0

    d2 = (r_mp + r_p_mp) ** 2 + (z_mp - z_p_mp) ** 2
    m = 4 * r_mp * r_p_mp / d2
    if m >= mp.mpf("1") - mp.mpf("1e-12"):
        m = mp.mpf("1") - mp.mpf("1e-12")

    k = mp.sqrt(m)
    K = mp.ellipk(m)
    E = mp.ellipe(m)
    coef = mp.sqrt(r_p_mp / r_mp) / mp.pi
    val = coef * ((2 / k - k) * K - (2 / k) * E)

    return dd_from_mpmath(val)


def axisym_kernel_dd_batch_cpu(r_list_hi, r_list_lo, z_list_hi, z_list_lo,
                                 r_p_list_hi, r_p_list_lo,
                                 z_p_list_hi, z_p_list_lo, dps=40):
    """Batch compute axisym G for many pairs via mpmath. Returns DD arrays.

    Optionally parallelize via multiprocessing.
    """
    n = r_list_hi.size
    G_hi = np.zeros(n)
    G_lo = np.zeros(n)
    for i in range(n):
        gh, gl = axisym_kernel_dd_mpmath(
            r_list_hi[i], r_list_lo[i], z_list_hi[i], z_list_lo[i],
            r_p_list_hi[i], r_p_list_lo[i], z_p_list_hi[i], z_p_list_lo[i],
            dps=dps)
        G_hi[i] = gh
        G_lo[i] = gl
    return G_hi, G_lo


def main():
    """Quick test of DD axisym kernel."""
    print("=" * 72)
    print("DD axisym Green's function (mpmath elliptic) test")
    print("=" * 72)

    # Test: G(r=1, z=0; r'=1.1, z'=0.2) at FP64 vs DD vs mpmath
    r, z = 1.0, 0.0
    r_p, z_p = 1.1, 0.2

    # DD evaluation
    g_hi, g_lo = axisym_kernel_dd_mpmath(np.array(r), np.array(0.0),
                                            np.array(z), np.array(0.0),
                                            np.array(r_p), np.array(0.0),
                                            np.array(z_p), np.array(0.0))
    print(f"\n  DD G({r}, {z}; {r_p}, {z_p}) = "
          f"{float(g_hi):.16e} + {float(g_lo):.6e}")

    # FP64 reference (scipy.special)
    import scipy.special as sp
    d2 = (r + r_p) ** 2 + (z - z_p) ** 2
    m_fp = 4.0 * r * r_p / d2
    k_fp = math.sqrt(m_fp)
    K_fp = sp.ellipk(m_fp)
    E_fp = sp.ellipe(m_fp)
    coef_fp = math.sqrt(r_p / r) / math.pi
    g_fp = coef_fp * ((2.0 / k_fp - k_fp) * K_fp - (2.0 / k_fp) * E_fp)
    print(f"  FP64 scipy.special G = {g_fp:.16e}")

    # Compare
    diff = abs(g_hi - g_fp)
    rel = diff / abs(g_fp)
    print(f"\n  DD_hi vs FP64: abs diff = {diff:.3e}, rel = {rel:.3e}")
    print(f"  DD_lo (sub-FP64 precision) = {g_lo:.3e}")


if __name__ == "__main__":
    main()
