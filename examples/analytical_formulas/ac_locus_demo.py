"""Major / minor axis of an AC vector quantity's time-locus ellipse.

Reference: Part 5 (SA-04-44 / RM-04-68, 2004), eq 29-37.

The closed form is verified against a brute-force time sweep for a
randomly-chosen 3-component complex phasor. Cross-section: a noisy
"flux-density" sample point in an AC eddy-current run might give
B_x = 0.6 + 0.1j, B_y = -0.4 + 0.3j, B_z = 0.1 - 0.05j (in T).
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from radia.analytical_formulas import ac_locus_axes


def brute_force_axes(B: np.ndarray, n_t: int = 20000):
    """Sample |b(t)| on a fine time grid and return its (max, min)."""
    omega_t = np.linspace(0.0, 2.0 * np.pi, n_t, endpoint=False)
    bt = np.real(B[None, :] * np.exp(1j * omega_t)[:, None])  # (n_t, 3)
    mag = np.linalg.norm(bt, axis=1)
    return float(mag.max()), float(mag.min())


def main() -> None:
    rng = np.random.default_rng(seed=2025)

    cases = {
        "linear (along x)": np.array([1.0, 0.0, 0.0], dtype=complex),
        "circular (xy)":     np.array([1.0, 1j, 0.0]),
        "elliptical 2:1":    np.array([2.0, 1j, 0.0]),
        "tilted 3D ellipse": np.array([0.6 + 0.1j, -0.4 + 0.3j, 0.1 - 0.05j]),
        "random sample":     rng.standard_normal(3) + 1j * rng.standard_normal(3),
    }

    hdr = f'{"case":<22} {"|B|_max (analytic)":>20} {"|B|_max (sweep)":>18}  {"|B|_min (analytic)":>20} {"|B|_min (sweep)":>18}'
    print(hdr)
    print("-" * len(hdr))
    for name, B in cases.items():
        mx_an, mn_an = ac_locus_axes(B)
        mx_bf, mn_bf = brute_force_axes(B)
        print(
            f"{name:<22} {mx_an:>20.6f} {mx_bf:>18.6f}  "
            f"{mn_an:>20.6f} {mn_bf:>18.6f}"
        )

    # Optional plot of the locus for the tilted case.
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    except ImportError:
        print("\n(matplotlib not installed -- skipping plot)")
        return
    B = cases["tilted 3D ellipse"]
    omega_t = np.linspace(0.0, 2.0 * np.pi, 200)
    bt = np.real(B[None, :] * np.exp(1j * omega_t)[:, None])
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(bt[:, 0], bt[:, 1], bt[:, 2])
    ax.set_xlabel("B_x")
    ax.set_ylabel("B_y")
    ax.set_zlabel("B_z")
    ax.set_title("3D AC vector locus (tilted ellipse case)")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "ac_locus_demo.png")
    fig.savefig(out, dpi=150)
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
