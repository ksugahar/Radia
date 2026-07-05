"""Cylinder Gx fingerprint: (homogeneity, peak-J) Pareto front via folded
Tikhonov + ACA+TSVD, with the cylinder LENGTH as the geometry lever.

The same folded Tikhonov + ACA routine (RegularizedTSVD) as the planar demo,
but on the cylindrical SF basis-loop kernel (the standard MRI / shim geometry).
Reuses the committed cylinder primitives from demo_sf_to_peec_gx.py:
loop_corners / _loop_Hz / make_dsv / _cyl_laplacian.

Target = Gx transverse-gradient fingerprint (Bz = Gx * x over the DSV).
Geometry lever = cylinder length L: a longer cylinder gives the azimuthal
current more axial room -> lower peak.  Unlike the planar former (monotone),
the cylinder length has an OPTIMUM: once L covers the DSV (loops beyond
~+/- 2*dsv barely reach the target) a longer cylinder stops helping.  nz is
scaled with L so the grid spacing dz stays ~constant (fair peak comparison).

Outputs JSON + PNG next to this script.
"""
import os
import sys
import json
import argparse

import numpy as np

from radia_mcp.figure import lab_figure, save_lab_figure

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "..", "..", "src"))

from demo_sf_to_peec_gx import loop_corners, _loop_Hz, make_dsv, _cyl_laplacian
from radia.stream_function import aca_tsvd, RegularizedTSVD


def cyl_front(L, a, Gx, dsv, nphi, dz_target, ndsv, lams):
    """Folded-Tikhonov (homogeneity, peak) front for one cylinder length."""
    nz = max(10, int(round(L / dz_target)))
    phi_grid = np.linspace(0.0, 2 * np.pi, nphi, endpoint=False)
    z_grid = np.linspace(-L / 2, L / 2, nz)
    dphi = 2 * np.pi / nphi
    dz = L / (nz - 1)
    corners = [loop_corners(a, pp, zq, dphi, dz)
               for zq in z_grid for pp in phi_grid]
    N = len(corners)
    obs = make_dsv(dsv, ndsv)
    M = len(obs)
    B = Gx * obs[:, 0]
    A = np.array([[_loop_Hz(obs[i], corners[j]) for j in range(N)]
                  for i in range(M)])
    S = _cyl_laplacian(nz, nphi)
    md = float(np.mean(np.diag(A.T @ A)))
    base = aca_tsvd(M, N, lambda i, j: float(A[i, j]),
                    modes=M, kmax=min(M, N), aca_eps=1.0e-10)
    reg = RegularizedTSVD.from_stiffness(base, S)

    def peak(psi):
        P = psi.reshape(nz, nphi)
        dpz = np.gradient(P, dz, axis=0)                          # z free
        dpp = (np.roll(P, -1, 1) - np.roll(P, 1, 1)) / (2 * dphi)  # phi periodic
        return float(np.sqrt(dpz ** 2 + (dpp / a) ** 2).max())

    pts = []
    for lam in lams:
        psi = reg.solve(B, alpha=lam * md)
        misfit = float(np.linalg.norm(A @ psi - B)
                       / (np.linalg.norm(B) + 1e-30))
        pts.append(dict(lam=float(lam), misfit=misfit, peak=peak(psi)))
    return pts, M, N


def nondominated(pts):
    s = sorted(pts, key=lambda p: (p["misfit"], p["peak"]))
    out, best = [], float("inf")
    for p in s:
        if p["peak"] < best - 1e-15:
            out.append(p); best = p["peak"]
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--radius", type=float, default=0.15)
    ap.add_argument("--Gx", type=float, default=1.0)
    ap.add_argument("--dsv", type=float, default=0.08)
    ap.add_argument("--nphi", type=int, default=24)
    ap.add_argument("--ndsv", type=int, default=5)
    ap.add_argument("--dz-target", type=float, default=0.03)
    ap.add_argument("--lengths", type=float, nargs="+",
                    default=[0.32, 0.40, 0.50, 0.62, 0.78])
    ap.add_argument("--out", default=os.path.join(_HERE, "pareto_cylinder"))
    args = ap.parse_args()

    lams = np.concatenate([[0.0], np.logspace(-4, 1.5, 11)])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = lab_figure(embed_width_cm=8.0)
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(args.lengths)))

    per_L, allpts = {}, []
    for L, col in zip(args.lengths, cmap):
        pts, M, N = cyl_front(L, args.radius, args.Gx, args.dsv, args.nphi,
                              args.dz_target, args.ndsv, lams)
        per_L[L] = pts
        allpts += pts
        nd = nondominated(pts)
        ax.plot([p["misfit"] * 100 for p in nd], [p["peak"] for p in nd],
                "o-", color=col, ms=4, lw=1.2, label=f"L={L*100:.0f} cm")
        print(f"L={L:.2f} m  M={M} N={N}  exact-homog peak={nd[0]['peak']:.3e} "
              f"(misfit {nd[0]['misfit']*100:.1e}%)", flush=True)

    env = nondominated(allpts)
    ax.plot([p["misfit"] * 100 for p in env], [p["peak"] for p in env],
            "k--", lw=2.0, label="envelope (length-pushed)")
    ax.set_xscale("log")
    ax.set_xlabel("inhomogeneity ||A psi - B|| / ||B|| (%)")
    ax.set_ylabel("peak surface current density max|grad psi|")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(alpha=0.3)
    save_lab_figure(fig, args.out, embed_width_cm=8.0)

    peaks0 = {L: nondominated(per_L[L])[0]["peak"] for L in args.lengths}
    Lbest = min(peaks0, key=peaks0.get)
    json.dump({"radius": args.radius, "Gx": args.Gx, "dsv": args.dsv,
               "lengths": args.lengths, "per_L": {str(k): v
                                                  for k, v in per_L.items()},
               "envelope": env, "best_length": Lbest},
              open(args.out + ".json", "w"), indent=2)
    print(f"\nlength lever (exact-homog peak): "
          f"{peaks0[args.lengths[0]]:.3e} -> {peaks0[Lbest]:.3e} "
          f"(optimum at L={Lbest*100:.0f} cm)", flush=True)
    print(f"saved {args.out}.json / .png", flush=True)


if __name__ == "__main__":
    main()
