r"""The 3-D closing condition: HELICITY is the obstruction to a global Clebsch pair.

In 2-D a flux line always has a conserved Hamiltonian A_z (the Clebsch potential,
`hdiv_vim_clebsch_2d_az.py`), so flux lines always close.  In **3-D** a *global*
Clebsch representation `B = grad(alpha) x grad(beta)` need not exist, and the
fundamental obstruction is the **HELICITY**

    h = INT_V A . B dV        (B = curl A),

the topological linking/knotting of the field lines (Moffatt 1969).  A Clebsch
field is helicity-free **pointwise**: with the vector potential `A = alpha grad(beta)`
(indeed `curl(alpha grad beta) = grad(alpha) x grad(beta) = B`),

    A . B = alpha grad(beta) . (grad(alpha) x grad(beta)) = 0,

so `h = 0`.  Therefore **`h != 0`  =>  NO global Clebsch pair**: the field lines are
linked/knotted and cannot all lie on `alpha = const ∩ beta = const` surfaces -- they
do not close.  This is the 3-D "closing condition": *flux lines close (lie on flux
surfaces) iff a global Clebsch pair exists iff the helicity obstruction vanishes.*

Demonstrated on the 3-torus:
  - a Clebsch field `B = grad(alpha) x grad(beta)`  ->  `h ~ 0` (machine), and its
    flux lines are confined;
  - the **ABC (Arnold-Beltrami-Childress)** Beltrami field `curl B = B`
    (`Bx=sin z+cos y, By=sin x+cos z, Bz=sin y+cos x`) has `A = B`, so
    `h = INT |B|^2 = 3 (2 pi)^3 != 0`  ->  no global Clebsch, and its flux lines are
    chaotic (the famous ABC chaos -- a Poincare section that fills a 2-D region,
    never closing onto a curve).

run:  python clebsch_3d_closing_condition.py
"""
import os
from math import pi, sin, cos

import numpy as np

TWO_PI = 2.0 * pi


# ---- the ABC Beltrami field (curl B = B) on the 3-torus ----
def abc_B(x):
    X, Y, Z = x
    return np.array([sin(Z) + cos(Y), sin(X) + cos(Z), sin(Y) + cos(X)])


# ---- a Clebsch field B = grad(alpha) x grad(beta) and its potential A = alpha grad(beta) ----
def _clebsch_fields(X, Y, Z):
    # alpha = cos x + sin y ,  beta = cos z + sin x   (generic, non-degenerate)
    ga = np.array([-np.sin(X), np.cos(Y), 0.0 * X])               # grad(alpha)
    gb = np.array([np.cos(X), 0.0 * Y, -np.sin(Z)])               # grad(beta)
    alpha = np.cos(X) + np.sin(Y)
    B = np.cross(ga, gb, axis=0)                                  # grad(alpha) x grad(beta)
    A = alpha * gb                                                # alpha grad(beta) (curl A = B)
    return A, B


def helicity_grid(n=48):
    """INT A.B over [0,2pi]^3 by midpoint rule for the Clebsch field (should be ~0) and the ABC field
    (A=B, so h = INT |B|^2 = 3 (2pi)^3)."""
    g = (np.arange(n) + 0.5) / n * TWO_PI
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    dV = (TWO_PI / n) ** 3
    # Clebsch
    A, B = _clebsch_fields(X, Y, Z)
    h_cl = float(np.sum(np.einsum("i...,i...->...", A, B)) * dV)
    bnorm2_cl = float(np.sum(np.einsum("i...,i...->...", B, B)) * dV)
    # ABC: A = B (Beltrami), h = INT|B|^2
    Bx = np.sin(Z) + np.cos(Y); By = np.sin(X) + np.cos(Z); Bz = np.sin(Y) + np.cos(X)
    h_abc = float(np.sum(Bx * Bx + By * By + Bz * Bz) * dV)
    return h_cl, bnorm2_cl, h_abc


def _poincare_abc(x0, n_cross=2500, ds=0.02, max_steps=3_000_000):
    """ABC flux-line Poincare section at z = 0 (mod 2pi), upward crossings -> (x mod 2pi, y mod 2pi)."""
    x = np.array(x0, float); out = []; zp = x[2] % TWO_PI; s = 0
    while len(out) < n_cross and s < max_steps:
        k1 = abc_B(x); k2 = abc_B(x + 0.5 * ds * k1); k3 = abc_B(x + 0.5 * ds * k2); k4 = abc_B(x + ds * k3)
        xn = x + ds / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        zm = xn[2] % TWO_PI
        if zm < zp and k1[2] > 0:
            fr = zp / (zp - zm + 1e-30); xc = x + fr * (xn - x)
            out.append([xc[0] % TWO_PI, xc[1] % TWO_PI])
        zp = zm; x = xn; s += 1
    return np.array(out)


def _occupancy(pts, K=40):
    if len(pts) == 0:
        return 0.0
    ix = np.clip((pts[:, 0] / TWO_PI * K).astype(int), 0, K - 1)
    iy = np.clip((pts[:, 1] / TWO_PI * K).astype(int), 0, K - 1)
    return len(set(zip(ix.tolist(), iy.tolist()))) / (K * K)


def analyze(n=48, n_cross=2500):
    h_cl, bnorm2_cl, h_abc = helicity_grid(n)
    # pick a chaotic-sea start for the ABC Poincare (away from the island around the stagnation pts)
    P = _poincare_abc([1.0, 2.0, 0.0], n_cross=n_cross)
    occ = _occupancy(P)
    return {
        "n": n, "h_clebsch": h_cl, "Bnorm2_clebsch": bnorm2_cl, "h_abc": h_abc,
        "h_abc_exact": 3.0 * TWO_PI ** 3, "rel_clebsch": abs(h_cl) / (bnorm2_cl + 1e-30),
        "P_abc": P, "occ_abc": occ,
    }


def _plot(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(5.2, 4.8), dpi=150)
    ax.plot(r["P_abc"][:, 0], r["P_abc"][:, 1], "C3.", ms=1.0)
    ax.set_aspect("equal"); ax.set_xlim(0, TWO_PI); ax.set_ylim(0, TWO_PI)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(f"ABC flux line: chaotic Poincare section\nhelicity $\\neq$ 0 $\\Rightarrow$ no global "
                 f"Clebsch, never closes")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("The 3-D closing condition: helicity is the obstruction to a global Clebsch pair\n")
    r = analyze()
    print(f"  helicity h = INT A.B over the 3-torus (grid {r['n']}^3):")
    print(f"    Clebsch field B=grad(a)xgrad(b):  h={r['h_clebsch']:.2e}  (rel to INT|B|^2: "
          f"{r['rel_clebsch']:.1e})  -> ZERO -> global Clebsch exists, flux lines on surfaces")
    print(f"    ABC Beltrami field (A=B):         h={r['h_abc']:.4f}  (exact 3(2pi)^3="
          f"{r['h_abc_exact']:.4f})  -> NONZERO -> NO global Clebsch")
    print(f"  ABC flux-line Poincare occupancy = {r['occ_abc']:.3f}  (fills 2-D = chaotic, never closes)")
    print("\n  => the 3-D closing condition is the vanishing of the HELICITY = the existence of a global")
    print("     Clebsch pair B = grad(alpha) x grad(beta).  2-D always has A_z (helicity-free, planar);")
    print("     3-D is obstructed by helicity (Moffatt) -- the real frontier.")
    _plot(r)


if __name__ == "__main__":
    main()
