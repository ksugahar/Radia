"""Gauss-Legendre quadrature: nodes / weights table + convergence.

Reference: Part 3 §4 (SA-03-52 / RM-03-54, 2003), Table 3.

Prints the ``n = 1, ..., 24`` Gauss-Legendre nodes and weights and
shows convergence of a smooth-function integral as ``n`` increases.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import (
    gauss_legendre_integrate,
    gauss_legendre_nodes_weights,
)
from radia.analytical_formulas.gauss_legendre import N_MAX


def main() -> None:
    print("Gauss-Legendre nodes and weights (positive half, n = 1..6):")
    for n in (1, 2, 3, 4, 5, 6):
        x, w = gauss_legendre_nodes_weights(n)
        # Print only the non-negative half (the table is symmetric).
        half = x[x >= -1e-15]
        wh = w[x >= -1e-15]
        print(f"  n = {n}:")
        for xi, wi in zip(half, wh):
            print(f"    x = {xi:+.15f}    w = {wi:.15f}")
    print(f"  ...   table extends to n = {N_MAX}.")
    print()

    # Convergence of a smooth integral.
    # int_0^pi sin(x) dx = 2 ; the n-point rule is exact for polynomials of
    # degree 2 n - 1, and converges spectrally on smooth integrands.
    truth = 2.0
    print("Convergence of int_0^pi sin(x) dx  (truth = 2.0):")
    print(f"  {'n':>3}  {'estimate':>22}  {'abs error':>12}")
    for n in (1, 2, 3, 4, 5, 6, 8, 10, 12, 16):
        val = gauss_legendre_integrate(np.sin, 0.0, math.pi, n=n)
        err = abs(val - truth)
        print(f"  {n:>3}   {val:>22.15f}   {err:.2e}")

    # Plot: error vs n on a log scale.
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return
    ns = np.arange(1, N_MAX + 1)
    errs_sin = []
    errs_smooth = []
    for n in ns:
        val_sin = gauss_legendre_integrate(np.sin, 0.0, math.pi, n=n)
        val_smooth = gauss_legendre_integrate(
            lambda x: np.exp(-x ** 2), -3.0, 3.0, n=n
        )
        errs_sin.append(abs(val_sin - 2.0))
        # int_-3^3 exp(-x**2) dx = sqrt(pi) * erf(3)  (very close to sqrt(pi)).
        from scipy.special import erf
        truth_smooth = math.sqrt(math.pi) * erf(3.0)
        errs_smooth.append(abs(val_smooth - truth_smooth))

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.semilogy(ns, errs_sin, "o-", lw=1.5, label=r"$\int_0^\pi \sin(x)\, dx$")
    ax.semilogy(ns, errs_smooth, "s-", lw=1.5,
                label=r"$\int_{-3}^{3} e^{-x^2}\, dx$")
    ax.set_xlabel("n (number of quadrature points)")
    ax.set_ylabel("absolute error")
    ax.set_title("Gauss-Legendre spectral convergence on smooth integrands")
    ax.set_yscale("log")
    ax.grid(True, ls=":", alpha=0.5, which="both")
    ax.legend()
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "gauss_legendre_demo.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
