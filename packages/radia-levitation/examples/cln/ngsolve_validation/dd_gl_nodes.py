"""dd_gl_nodes.py — High-precision Gauss-Legendre nodes/weights for DD pipeline.

GL nodes from np.polynomial.legendre.leggauss are FP64 (~1e-15 precision).
For DD precision in subsequent integrations, refine via mpmath Newton.
"""
from __future__ import annotations

import numpy as np
import mpmath as mp

from dd_arithmetic import dd_from_mpmath


def gauss_legendre_mpf(n, dps=40):
    """Compute Gauss-Legendre nodes/weights at high precision via mpmath Newton.

    Returns lists of mp.mpf for nodes and weights, sorted ascending.
    """
    mp.mp.dps = dps
    # Initial guess from FP64
    n_f, w_f = np.polynomial.legendre.leggauss(n)
    nodes_mp = []
    weights_mp = []
    for i in range(n):
        x = mp.mpf(float(n_f[i]))
        # Newton iteration on P_n(x) = 0
        for _ in range(30):
            Pnm1 = mp.mpf(1)
            Pn = x
            for k in range(1, n):
                Pnp1 = ((2 * k + 1) * x * Pn - k * Pnm1) / mp.mpf(k + 1)
                Pnm1 = Pn
                Pn = Pnp1
            Pn_deriv = n * (x * Pn - Pnm1) / (x * x - 1)
            x_new = x - Pn / Pn_deriv
            if abs(x_new - x) < mp.mpf(10) ** (-dps + 2):
                x = x_new
                break
            x = x_new
        nodes_mp.append(x)
        # Weight: 2 / ((1 - x²) (P_n'(x))²)
        Pnm1 = mp.mpf(1)
        Pn = x
        for k in range(1, n):
            Pnp1 = ((2 * k + 1) * x * Pn - k * Pnm1) / mp.mpf(k + 1)
            Pnm1 = Pn
            Pn = Pnp1
        Pn_deriv = n * (x * Pn - Pnm1) / (x * x - 1)
        w = mp.mpf(2) / ((1 - x * x) * Pn_deriv * Pn_deriv)
        weights_mp.append(w)
    # Sort
    pairs = sorted(zip(nodes_mp, weights_mp), key=lambda p: float(p[0]))
    nodes_mp = [p[0] for p in pairs]
    weights_mp = [p[1] for p in pairs]
    return nodes_mp, weights_mp


def gauss_legendre_dd_arrays(n, dps=40):
    """GL nodes/weights as DD numpy arrays."""
    nodes_mp, weights_mp = gauss_legendre_mpf(n, dps)
    nodes_hi = np.zeros(n)
    nodes_lo = np.zeros(n)
    weights_hi = np.zeros(n)
    weights_lo = np.zeros(n)
    for i in range(n):
        nodes_hi[i], nodes_lo[i] = dd_from_mpmath(nodes_mp[i])
        weights_hi[i], weights_lo[i] = dd_from_mpmath(weights_mp[i])
    return nodes_hi, nodes_lo, weights_hi, weights_lo, nodes_mp, weights_mp


def main():
    """Verify high-precision GL nodes vs np.leggauss."""
    print("=" * 72)
    print("DD GL nodes/weights verification")
    print("=" * 72)
    for n in [4, 8, 12, 24]:
        n_f, w_f = np.polynomial.legendre.leggauss(n)
        nodes_mp, weights_mp = gauss_legendre_mpf(n, dps=40)
        max_node_err = max(abs(mp.mpf(float(n_f[i])) - nodes_mp[i])
                            for i in range(n))
        max_weight_err = max(abs(mp.mpf(float(w_f[i])) - weights_mp[i])
                              for i in range(n))
        print(f"  n={n}: max node err vs FP64 = {float(max_node_err):.3e}, "
              f"max weight err = {float(max_weight_err):.3e}")
        # Sanity: weights sum to 2
        wsum = sum(weights_mp)
        print(f"        sum(weights) = {float(wsum):.20f}, "
              f"err = {float(abs(wsum - 2)):.3e}")


if __name__ == "__main__":
    main()
