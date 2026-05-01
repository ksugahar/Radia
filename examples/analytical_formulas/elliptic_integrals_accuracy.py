"""Hastings polynomial approximations for K(k), E(k): accuracy table.

Reference: Part 3 §3 (SA-03-52 / RM-03-54, 2003), eq 13 + Tables 1-2.

Compares the degree-2 and degree-4 Hastings polynomial approximations
against ``scipy.special.ellipk`` / ``ellipe`` (machine-precision
reference) and reproduces Fig. 2 of the PDF.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import (
    E_hastings_2,
    E_hastings_4,
    K_hastings_2,
    K_hastings_4,
)


def main() -> None:
    from scipy.special import ellipe as scipy_E, ellipk as scipy_K

    ks = np.linspace(0.0, 0.99, 200)
    K_ref = scipy_K(ks * ks)
    E_ref = scipy_E(ks * ks)
    K2 = K_hastings_2(ks)
    K4 = K_hastings_4(ks)
    E2 = E_hastings_2(ks)
    E4 = E_hastings_4(ks)

    print(f'{"k":>6}  {"K err deg2":>12}  {"K err deg4":>12}  '
          f'{"E err deg2":>12}  {"E err deg4":>12}')
    for k in (0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99):
        Kr = float(scipy_K(k * k))
        Er = float(scipy_E(k * k))
        print(f"{k:>6.2f}  {abs(K_hastings_2(k) - Kr):>12.2e}  "
              f"{abs(K_hastings_4(k) - Kr):>12.2e}  "
              f"{abs(E_hastings_2(k) - Er):>12.2e}  "
              f"{abs(E_hastings_4(k) - Er):>12.2e}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return

    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].semilogy(ks, np.abs(K2 - K_ref), label="K, degree 2")
    axes[0].semilogy(ks, np.abs(K4 - K_ref), label="K, degree 4")
    axes[0].set_ylabel("absolute error in K")
    axes[0].grid(True, ls=":", alpha=0.5, which="both")
    axes[0].legend()
    axes[0].set_title("Hastings polynomial approximation accuracy "
                      "(ref: scipy.special)")

    axes[1].semilogy(ks, np.abs(E2 - E_ref), label="E, degree 2")
    axes[1].semilogy(ks, np.abs(E4 - E_ref), label="E, degree 4")
    axes[1].set_ylabel("absolute error in E")
    axes[1].set_xlabel("k")
    axes[1].grid(True, ls=":", alpha=0.5, which="both")
    axes[1].legend()

    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__),
                       "elliptic_integrals_accuracy.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
